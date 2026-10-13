"""Tabular MVP models: histogram GBM + logistic baseline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from critical_minerals_aster.prospectivity.metrics import (
    FEATURE_COLS_ALL,
    FEATURE_COLS_HSI,
    FEATURE_COLS_TIR,
    FEATURE_COLS_CTX,
    evaluate_site,
    site_mean_metrics,
)
from critical_minerals_aster.prospectivity.split import era_block_split, site_group_kfold


@dataclass
class FitResult:
    model: Any
    feature_names: list[str]
    kind: str


def _make_gbm(**kwargs):
    from sklearn.ensemble import HistGradientBoostingClassifier

    params = dict(
        max_depth=4,
        learning_rate=0.08,
        max_iter=120,
        l2_regularization=0.1,
        random_state=0,
    )
    params.update(kwargs)
    return HistGradientBoostingClassifier(**params)


def _make_logistic(**kwargs):
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    params = dict(max_iter=500, class_weight="balanced", random_state=0)
    params.update(kwargs)
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(**params)),
    ])


def fit_model(
    X: np.ndarray,
    y: np.ndarray,
    *,
    kind: str = "gbm",
    feature_names: list[str] | None = None,
) -> FitResult:
    if kind == "logistic":
        model = _make_logistic()
    else:
        model = _make_gbm()
    model.fit(X, y)
    return FitResult(model=model, feature_names=list(feature_names or []), kind=kind)


def predict_proba(fit: FitResult, X: np.ndarray) -> np.ndarray:
    if hasattr(fit.model, "predict_proba"):
        return fit.model.predict_proba(X)[:, 1]
    # decision_function fallback
    s = fit.model.decision_function(X)
    return 1.0 / (1.0 + np.exp(-s))


def feature_importance(fit: FitResult) -> pd.DataFrame | None:
    m = fit.model
    names = fit.feature_names
    if hasattr(m, "feature_importances_"):
        imp = m.feature_importances_
    elif hasattr(m, "named_steps") and hasattr(m.named_steps.get("clf"), "coef_"):
        imp = np.abs(m.named_steps["clf"].coef_).ravel()
    else:
        return None
    return pd.DataFrame({"feature": names, "importance": imp}).sort_values(
        "importance", ascending=False
    )


def ablation_columns(mode: str) -> list[str]:
    mode = mode.lower()
    if mode == "tir":
        return FEATURE_COLS_TIR + FEATURE_COLS_CTX
    if mode == "hsi":
        return FEATURE_COLS_HSI + ["hsi_missing"]
    return FEATURE_COLS_ALL


def _xy(df: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    from critical_minerals_aster.prospectivity.features import model_matrix

    X, y, _ = model_matrix(df, cols=cols)
    return X, y


def run_site_cv(
    df: pd.DataFrame,
    *,
    n_splits: int = 5,
    kind: str = "gbm",
    ablation: str = "all",
    random_state: int = 0,
) -> tuple[pd.DataFrame, dict]:
    """Site-blocked GroupKFold; returns per-test-site metrics + site-mean summary."""
    cols = ablation_columns(ablation)
    rows = []
    for fold, (tr, te, test_sites) in enumerate(
        site_group_kfold(df["site_id"], n_splits=n_splits, random_state=random_state)
    ):
        train_df = df.iloc[tr]
        test_df = df.iloc[te]
        Xtr, ytr = _xy(train_df, cols)
        if len(np.unique(ytr)) < 2:
            continue
        fit = fit_model(Xtr, ytr, kind=kind, feature_names=cols)
        for sid in test_sites:
            sub = test_df[test_df["site_id"] == sid]
            if sub.empty:
                continue
            Xte, yte = _xy(sub, cols)
            scores = predict_proba(fit, Xte)
            rows.append(evaluate_site(yte, scores, site_id=sid) | {
                "fold": fold,
                "ablation": ablation,
                "model": kind,
            })
    per_site = pd.DataFrame(rows)
    summary = site_mean_metrics(rows)
    summary["ablation"] = ablation
    summary["model"] = kind
    summary["splitter"] = "site_group_kfold"
    return per_site, summary


def run_era_block(
    df: pd.DataFrame,
    *,
    kind: str = "gbm",
    ablation: str = "all",
    cutoff: int = 1950,
) -> dict:
    """Era-blocked stress test; report separately from site-CV (never average)."""
    cols = ablation_columns(ablation)
    split = era_block_split(df, cutoff=cutoff)
    train_df = df.iloc[split.train_idx]
    test_df = df.iloc[split.test_idx]
    out: dict[str, Any] = {
        "splitter": "era_block",
        "cutoff": cutoff,
        "ablation": ablation,
        "model": kind,
        "n_train": int(len(train_df)),
        "n_test": int(len(test_df)),
        "n_train_pos": split.n_train_pos,
        "n_test_pos": split.n_test_pos,
        "n_undated_pos": split.n_undated_pos,
        "n_held_out": int(len(split.held_out_idx)),
    }
    if split.n_train_pos < 5 or split.n_test_pos < 3:
        out["status"] = "underpowered"
        out["note"] = "Too few dated positives for a stable era-block estimate"
        return out

    Xtr, ytr = _xy(train_df, cols)
    if len(np.unique(ytr)) < 2:
        out["status"] = "single_class_train"
        return out
    fit = fit_model(Xtr, ytr, kind=kind, feature_names=cols)

    # Evaluate per site on post-1950 positives plus backgrounds at those sites
    # (test_df alone is positives-only under the default era_role rule).
    test_sites = set(test_df["site_id"].astype(str))
    eval_df = pd.concat(
        [
            test_df,
            df[(df["label"].astype(int) == 0) & df["site_id"].astype(str).isin(test_sites)],
        ],
        ignore_index=True,
    ).drop_duplicates(subset=["site_id", "row", "col"], keep="first")

    site_rows = []
    for sid, sub in eval_df.groupby("site_id"):
        Xte, yte = _xy(sub, cols)
        if yte.sum() == 0 or len(np.unique(yte)) < 2:
            continue
        scores = predict_proba(fit, Xte)
        site_rows.append(evaluate_site(yte, scores, site_id=str(sid)))
    out["per_site"] = site_rows
    out["site_mean"] = site_mean_metrics(site_rows)
    out["status"] = "ok"
    # Also pooled test AP for diagnostics only (not headline)
    Xte_all, yte_all = _xy(test_df, cols)
    scores_all = predict_proba(fit, Xte_all)
    from critical_minerals_aster.prospectivity.metrics import average_precision_safe

    out["pooled_pr_auc_diagnostic_only"] = average_precision_safe(yte_all, scores_all)
    out["importance"] = (
        feature_importance(fit).to_dict(orient="records")
        if feature_importance(fit) is not None
        else []
    )
    return out
