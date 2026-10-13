"""Blocked validation splitters: site GroupKFold and era pre/post-1950."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

import numpy as np
import pandas as pd

ERA_CUTOFF = 1950


@dataclass(frozen=True)
class EraBlockSplit:
    """Train on pre-cutoff dated positives; test on post-cutoff dated positives.

    Undated rows are held out of both train and test positive sets. Background
    rows may follow the train/test membership of their paired deposit site, or
    be assigned via ``bg_mask`` columns already in the frame.
    """

    cutoff: int
    train_idx: np.ndarray
    test_idx: np.ndarray
    held_out_idx: np.ndarray
    n_train_pos: int
    n_test_pos: int
    n_undated_pos: int


def site_group_kfold(
    site_ids: pd.Series | np.ndarray,
    *,
    n_splits: int = 5,
    shuffle: bool = True,
    random_state: int = 0,
) -> Iterator[tuple[np.ndarray, np.ndarray, list[str]]]:
    """Yield (train_idx, test_idx, test_site_ids) with sites never mixed.

    Uses sklearn GroupKFold when available; falls back to a deterministic
    leave-group-out cycle when n_groups < n_splits.
    """
    groups = pd.Series(site_ids).astype(str).reset_index(drop=True)
    unique = sorted(groups.unique())
    n_groups = len(unique)
    if n_groups < 2:
        raise ValueError("site_group_kfold requires at least 2 distinct site_id values")

    n_splits = min(n_splits, n_groups)
    try:
        from sklearn.model_selection import GroupKFold

        gkf = GroupKFold(n_splits=n_splits)
        X_dummy = np.zeros((len(groups), 1))
        y_dummy = np.zeros(len(groups))
        for train_idx, test_idx in gkf.split(X_dummy, y_dummy, groups=groups.to_numpy()):
            test_sites = sorted(groups.iloc[test_idx].unique().tolist())
            # Sanity: no site appears in both folds
            train_sites = set(groups.iloc[train_idx].unique())
            assert train_sites.isdisjoint(test_sites), "site leakage in GroupKFold"
            yield train_idx, test_idx, test_sites
        return
    except ImportError:
        pass

    # Deterministic fallback without sklearn
    if shuffle:
        rng = np.random.default_rng(random_state)
        order = list(unique)
        rng.shuffle(order)
    else:
        order = list(unique)
    folds = [order[i::n_splits] for i in range(n_splits)]
    for i in range(n_splits):
        test_sites = folds[i]
        train_sites = [s for j, fold in enumerate(folds) if j != i for s in fold]
        train_idx = np.where(groups.isin(train_sites))[0]
        test_idx = np.where(groups.isin(test_sites))[0]
        yield train_idx, test_idx, list(test_sites)


def era_block_split(
    df: pd.DataFrame,
    *,
    year_col: str = "disc_yr",
    label_col: str = "label",
    cutoff: int = ERA_CUTOFF,
) -> EraBlockSplit:
    """Split positives by discovery year; backgrounds stay with train by default.

    Rules
    -----
    - Positive + year < cutoff → train
    - Positive + year >= cutoff → test
    - Positive + year missing → held out (neither)
    - Background (label==0): included in train (stratum-matched hard negatives
      for the pre-cutoff regime). Backgrounds are not used as era-test
      positives; optional ``era_role`` column can override.
    """
    if year_col not in df.columns:
        raise KeyError(f"era_block_split requires column {year_col!r}")
    if label_col not in df.columns:
        raise KeyError(f"era_block_split requires column {label_col!r}")

    idx = np.arange(len(df))
    y = df[year_col]
    lab = df[label_col].astype(int)
    is_pos = lab == 1
    dated = y.notna()

    train_pos = is_pos & dated & (y < cutoff)
    test_pos = is_pos & dated & (y >= cutoff)
    undated_pos = is_pos & ~dated
    bg = ~is_pos

    if "era_role" in df.columns:
        train_mask = (df["era_role"] == "train") | train_pos
        test_mask = (df["era_role"] == "test") | test_pos
        held = ~(train_mask | test_mask)
    else:
        train_mask = train_pos | bg
        test_mask = test_pos
        held = undated_pos

    return EraBlockSplit(
        cutoff=cutoff,
        train_idx=idx[train_mask.to_numpy()],
        test_idx=idx[test_mask.to_numpy()],
        held_out_idx=idx[held.to_numpy()],
        n_train_pos=int(train_pos.sum()),
        n_test_pos=int(test_pos.sum()),
        n_undated_pos=int(undated_pos.sum()),
    )


def assert_no_site_leakage(train_sites: set[str], test_sites: set[str]) -> None:
    overlap = train_sites & test_sites
    if overlap:
        raise AssertionError(f"site leakage: {sorted(overlap)}")
