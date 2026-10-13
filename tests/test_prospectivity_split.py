"""Blocked validation splitters: no site leakage; era respects disc_yr."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from critical_minerals_aster.prospectivity.split import (
    ERA_CUTOFF,
    era_block_split,
    site_group_kfold,
)


def _synthetic_multi_site(n_sites: int = 6, cells_per: int = 20) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(0)
    for i in range(n_sites):
        sid = f"site_{i}"
        for j in range(cells_per):
            year = rng.choice([1900, 1940, 1960, 1980, np.nan])
            rows.append({
                "site_id": sid,
                "label": int(j < 4),
                "disc_yr": year,
                "x": float(j),
            })
    return pd.DataFrame(rows)


def test_group_kfold_never_mixes_sites():
    df = _synthetic_multi_site(6, 15)
    seen_test = set()
    for train_idx, test_idx, test_sites in site_group_kfold(df["site_id"], n_splits=3):
        train_sites = set(df.iloc[train_idx]["site_id"])
        test_set = set(test_sites)
        assert train_sites.isdisjoint(test_set)
        assert set(df.iloc[test_idx]["site_id"]) == test_set
        # All rows for a test site are in test_idx
        for sid in test_sites:
            site_rows = set(df.index[df["site_id"] == sid])
            assert site_rows.issubset(set(test_idx))
        seen_test |= test_set
    assert seen_test == set(df["site_id"].unique())


def test_group_kfold_requires_two_sites():
    df = pd.DataFrame({"site_id": ["a"] * 10})
    with pytest.raises(ValueError):
        list(site_group_kfold(df["site_id"], n_splits=2))


def test_era_split_respects_disc_yr():
    df = pd.DataFrame({
        "site_id": ["a", "a", "b", "b", "c", "c"],
        "label": [1, 0, 1, 0, 1, 0],
        "disc_yr": [1920, np.nan, 1975, np.nan, np.nan, np.nan],
    })
    split = era_block_split(df, cutoff=ERA_CUTOFF)
    train = df.iloc[split.train_idx]
    test = df.iloc[split.test_idx]
    held = df.iloc[split.held_out_idx]

    # Pre-1950 positive in train; backgrounds may accompany them
    pre = train[(train["label"] == 1) & train["disc_yr"].notna()]
    assert len(pre) == 1
    assert (pre["disc_yr"] < ERA_CUTOFF).all()

    # Post-1950 positive in test
    assert split.n_test_pos == 1
    assert (test.loc[test["label"] == 1, "disc_yr"] >= ERA_CUTOFF).all()

    # Undated positive held out
    assert split.n_undated_pos == 1
    assert (held["label"] == 1).sum() == 1
    assert held.loc[held["label"] == 1, "disc_yr"].isna().all()


def test_era_split_backgrounds_in_train_by_default():
    df = pd.DataFrame({
        "label": [1, 0, 0, 1],
        "disc_yr": [1910, np.nan, np.nan, 1960],
    })
    split = era_block_split(df)
    assert set(split.train_idx.tolist()) >= {0, 1, 2}  # pre pos + both bg
    assert 3 in set(split.test_idx.tolist())
