"""Tests for paired ASTER TIR vs GEMx SWIR deposit proximity metrics."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gemx_deposit_comparison import (
    paired_sensor_tests,
    site_direction_test,
    summarize_mask_at_points,
)


def test_summarize_mask_at_points_uses_only_shared_footprint():
    """Catches counting deposits outside the common TIR–SWIR footprint."""
    overlap = np.ones((3, 4), dtype=bool)
    overlap[2, 3] = False
    mask = np.zeros_like(overlap)
    mask[:, 0] = True

    result = summarize_mask_at_points(
        mask,
        overlap,
        rows=np.array([0, 1, 2, 2]),
        cols=np.array([0, 1, 2, 3]),
        pixel_size_m=90.0,
    )

    assert result["n_deposits"] == 3
    assert result["n_hits"] == 1
    assert result["hit_rate"] == pytest.approx(1 / 3)
    assert result["coverage"] == pytest.approx(3 / 11)
    assert result["enrichment"] == pytest.approx((1 / 3) / (3 / 11))
    assert result["median_distance_m"] == pytest.approx(90.0)


def test_paired_sensor_tests_count_directional_discordance():
    """Catches swapping TIR-only and SWIR-only deposit hits."""
    result = paired_sensor_tests(
        tir_hits=np.array([True, True, True, False, False]),
        swir_hits=np.array([False, False, True, True, False]),
        tir_distances_m=np.array([0.0, 0.0, 0.0, 90.0, 180.0]),
        swir_distances_m=np.array([90.0, 180.0, 0.0, 0.0, 270.0]),
    )

    assert result["tir_only_hits"] == 2
    assert result["swir_only_hits"] == 1
    assert result["mcnemar_exact_p"] == pytest.approx(1.0)
    assert result["median_distance_delta_m_swir_minus_tir"] == pytest.approx(90.0)


def test_site_direction_test_ignores_ties_and_respects_lower_is_better():
    """Catches treating larger distances as better or ties as sensor wins."""
    result = site_direction_test(
        tir=np.array([90.0, 100.0, 50.0, 20.0]),
        swir=np.array([180.0, 100.0, 70.0, 10.0]),
        lower_is_better=True,
    )

    assert result["tir_site_wins"] == 2
    assert result["swir_site_wins"] == 1
    assert result["site_ties"] == 1
    assert result["site_sign_test_p"] == pytest.approx(1.0)
