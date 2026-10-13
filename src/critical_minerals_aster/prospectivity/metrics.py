"""Site-mean evaluation metrics (never pool pixels across sites as headline)."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd


FEATURE_COLS_TIR = [
    "qi",
    "ci",
    "mi",
    "qi_rank",
    "ci_rank",
    "mi_rank",
    "score",
    "dist_strong_m",
]
FEATURE_COLS_HSI = [
    "aloh",
    "kaolinite",
    "alunite",
    "aloh_rank",
    "kaolinite_rank",
    "alunite_rank",
    "clay_mask",
]
FEATURE_COLS_CTX = [
    "stratum_code",
    "hsi_missing",
]
FEATURE_COLS_ALL = FEATURE_COLS_TIR + FEATURE_COLS_HSI + FEATURE_COLS_CTX


def average_precision_safe(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """PR-AUC; returns NaN when only one class is present."""
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    if y_true.size == 0 or len(np.unique(y_true)) < 2:
        return float("nan")
    try:
        from sklearn.metrics import average_precision_score

        return float(average_precision_score(y_true, y_score))
    except ImportError:
        # Rank-based AP fallback
        order = np.argsort(-y_score)
        yt = y_true[order]
        tp = np.cumsum(yt)
        fp = np.cumsum(1 - yt)
        recall = tp / max(yt.sum(), 1)
        precision = tp / np.maximum(tp + fp, 1)
        return float(np.trapz(precision, recall)) if recall.size > 1 else float("nan")


def site_hit_rate(y_true: np.ndarray, y_score: np.ndarray, *, top_frac: float = 0.1) -> float:
    """Fraction of positives falling in the top ``top_frac`` of scores."""
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    n_pos = int(y_true.sum())
    if n_pos == 0 or y_true.size == 0:
        return float("nan")
    k = max(1, int(np.ceil(top_frac * len(y_score))))
    top = np.argsort(-y_score)[:k]
    return float(y_true[top].sum() / n_pos)


def enrichment_vs_base(y_true: np.ndarray, y_score: np.ndarray, *, top_frac: float = 0.1) -> float:
    """(hit rate in top frac) / (base positive rate)."""
    base = float(np.mean(y_true)) if len(y_true) else float("nan")
    if not np.isfinite(base) or base <= 0:
        return float("nan")
    hr = site_hit_rate(y_true, y_score, top_frac=top_frac)
    # Convert hit-among-positives to precision-style enrichment of top band
    k = max(1, int(np.ceil(top_frac * len(y_score))))
    top = np.argsort(-y_score)[:k]
    prec = float(y_true[top].mean())
    return prec / base


def median_score_distance(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Median percentile-rank gap: positives vs all cells (higher = better separation)."""
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, dtype=float)
    if y_true.sum() == 0:
        return float("nan")
    order = y_score.argsort()
    ranks = np.empty_like(y_score, dtype=float)
    ranks[order] = np.linspace(0, 1, len(y_score), endpoint=False)
    return float(np.median(ranks[y_true == 1]) - np.median(ranks))


def evaluate_site(
    y_true: np.ndarray,
    y_score: np.ndarray,
    *,
    site_id: str,
    top_frac: float = 0.1,
) -> dict:
    return {
        "site_id": site_id,
        "n": int(len(y_true)),
        "n_pos": int(np.sum(y_true)),
        "pr_auc": average_precision_safe(y_true, y_score),
        "top10_hit_rate": site_hit_rate(y_true, y_score, top_frac=top_frac),
        "enrichment_top10": enrichment_vs_base(y_true, y_score, top_frac=top_frac),
        "median_rank_gap": median_score_distance(y_true, y_score),
    }


def site_mean_metrics(rows: Iterable[dict]) -> dict:
    """Average metrics across test sites (headline reporting rule)."""
    df = pd.DataFrame(list(rows))
    if df.empty:
        return {"n_sites": 0}
    out = {"n_sites": int(df["site_id"].nunique()), "n_rows": int(df["n"].sum())}
    for col in ("pr_auc", "top10_hit_rate", "enrichment_top10", "median_rank_gap"):
        if col in df.columns:
            out[f"site_mean_{col}"] = float(np.nanmean(df[col].to_numpy(dtype=float)))
            out[f"site_median_{col}"] = float(np.nanmedian(df[col].to_numpy(dtype=float)))
    return out
