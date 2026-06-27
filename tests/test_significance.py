"""Tests for significance.py — Wilson CI, uncertainty annotation, buffer sensitivity."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.significance import (  # noqa: E402
    wilson_ci,
    add_uncertainty_columns,
    buffer_sensitivity,
)


# ---------------------------------------------------------------------------
# wilson_ci
# ---------------------------------------------------------------------------


def test_wilson_ci_midrange():
    # 50 hits / 100 deposits: interval should straddle 0.5
    lo, hi = wilson_ci(50, 100)
    assert lo < 0.5 < hi
    assert 0.0 <= lo <= hi <= 1.0


def test_wilson_ci_zero_hits():
    # 0 hits: lower bound = 0, upper bound > 0
    lo, hi = wilson_ci(0, 20)
    assert lo == pytest.approx(0.0)
    assert hi > 0.0


def test_wilson_ci_all_hits():
    # All hits: upper bound = 1, lower bound < 1
    lo, hi = wilson_ci(20, 20)
    assert hi == pytest.approx(1.0)
    assert lo < 1.0


def test_wilson_ci_zero_n():
    # No deposits: return (0, 1) — maximum uncertainty
    lo, hi = wilson_ci(0, 0)
    assert lo == 0.0 and hi == 1.0


def test_wilson_ci_narrower_with_more_n():
    # Larger sample at same rate → narrower interval
    lo_small, hi_small = wilson_ci(5, 10)
    lo_large, hi_large = wilson_ci(50, 100)
    assert (hi_small - lo_small) > (hi_large - lo_large)


# ---------------------------------------------------------------------------
# add_uncertainty_columns
# ---------------------------------------------------------------------------


def _make_summary(n_bbox, n_in_zones, site_id="test") -> pd.DataFrame:
    """Minimal summary DataFrame matching the shape compute_site_summary produces."""
    return pd.DataFrame([
        {
            "site_id": site_id,
            "site_name": "Test Site",
            "row_type": "site",
            "n_deposits_bbox": n_bbox,
            "n_deposits_in_zones": n_in_zones,
            "hit_rate_pct": round(n_in_zones / n_bbox * 100, 1) if n_bbox else 0.0,
        }
    ])


def _make_zones_and_footprint():
    """Synthetic zones covering exactly 20% of a 100×100 km footprint.

    Both zone and footprint are in EPSG:32611 (metres).  Passing the footprint
    explicitly avoids the WGS84-to-UTM reprojection that would make a bbox-
    based coverage calculation inconsistent with the synthetic geometry.
    """
    import geopandas as gpd
    from shapely.geometry import box as shapely_box

    # Footprint: 100×100 units → area = 10 000
    footprint = gpd.GeoDataFrame(
        geometry=[shapely_box(0, 0, 100, 100)],
        crs="EPSG:32611",
    )
    # Zone covers 20×100 = 2 000 → 20% of footprint
    zone = gpd.GeoDataFrame(
        {"area_km2": [0.2]},
        geometry=[shapely_box(0, 0, 20, 100)],
        crs="EPSG:32611",
    )
    return zone, footprint


def test_add_uncertainty_columns_shape():
    zones, footprint = _make_zones_and_footprint()
    bbox = (-117.0, 37.0, -116.0, 38.0)  # not used when footprint is None
    summary = _make_summary(100, 30)

    result = add_uncertainty_columns(summary, zones, footprint, bbox)

    # All new columns must be present
    for col in ["null_hit_rate_pct", "hit_rate_ci_low", "hit_rate_ci_high", "p_binomial"]:
        assert col in result.columns, f"Missing column: {col}"


def test_add_uncertainty_columns_ci_ordering():
    zones, footprint = _make_zones_and_footprint()
    bbox = (-117.0, 37.0, -116.0, 38.0)
    summary = _make_summary(50, 20)

    result = add_uncertainty_columns(summary, zones, footprint, bbox)
    row = result.iloc[0]

    # CI must bracket the observed rate
    assert row["hit_rate_ci_low"] <= row["hit_rate_pct"] <= row["hit_rate_ci_high"]
    # Bounds must be percentages in [0, 100]
    assert 0.0 <= row["hit_rate_ci_low"] <= 100.0
    assert 0.0 <= row["hit_rate_ci_high"] <= 100.0


def test_add_uncertainty_columns_p_value_significant():
    """Observed rate well above null → p_binomial should be small."""
    zones, footprint = _make_zones_and_footprint()
    bbox = (-117.0, 37.0, -116.0, 38.0)
    # Zone covers ~20%; 40 of 50 deposits are in zones (80% hit rate)
    summary = _make_summary(50, 40)

    result = add_uncertainty_columns(summary, zones, footprint, bbox)
    assert result.iloc[0]["p_binomial"] < 0.001


def test_add_uncertainty_columns_p_value_not_significant():
    """Observed rate at null → p_binomial should not be significant."""
    zones, footprint = _make_zones_and_footprint()
    bbox = (-117.0, 37.0, -116.0, 38.0)
    # Zone covers ~20%; 10 of 50 deposits (20% = null rate exactly)
    summary = _make_summary(50, 10)

    result = add_uncertainty_columns(summary, zones, footprint, bbox)
    assert result.iloc[0]["p_binomial"] > 0.4


def test_add_uncertainty_columns_zero_deposits():
    """Zero deposits → no crash, p_binomial = 1.0."""
    zones, footprint = _make_zones_and_footprint()
    bbox = (-117.0, 37.0, -116.0, 38.0)
    summary = _make_summary(0, 0)

    result = add_uncertainty_columns(summary, zones, footprint, bbox)
    assert result.iloc[0]["p_binomial"] == 1.0


# ---------------------------------------------------------------------------
# buffer_sensitivity
# ---------------------------------------------------------------------------


def _make_point_deposits_and_zones():
    """Deposits just outside a zone at radius 0; inside at radius >= 300 m."""
    import geopandas as gpd
    from shapely.geometry import Point, box as shapely_box

    # Zone: square from x=0–100, y=0–100 in EPSG:32611
    zones = gpd.GeoDataFrame(
        {"area_km2": [0.01]},
        geometry=[shapely_box(0, 0, 100, 100)],
        crs="EPSG:32611",
    )
    # 4 deposits: 2 clearly inside (50,50 and 20,20), 2 just outside (at x=200,50 — 100m gap)
    deposits = gpd.GeoDataFrame(
        {"commod1": ["Gold"] * 4},
        geometry=[Point(50, 50), Point(20, 20), Point(200, 50), Point(200, 80)],
        crs="EPSG:32611",
    )
    return deposits, zones


def test_buffer_sensitivity_returns_correct_columns():
    deposits, zones = _make_point_deposits_and_zones()
    result = buffer_sensitivity(deposits, zones, radii_m=(0, 500))
    for col in ["radius_m", "n_deposits", "n_hits", "hit_rate_pct",
                "hit_rate_ci_low", "hit_rate_ci_high", "delta_pct"]:
        assert col in result.columns


def test_buffer_sensitivity_hits_increase_with_radius():
    """Points just outside zone at r=0 should become hits at larger radius."""
    deposits, zones = _make_point_deposits_and_zones()
    result = buffer_sensitivity(deposits, zones, radii_m=(0, 500))
    hits_0 = result.loc[result["radius_m"] == 0, "n_hits"].iloc[0]
    hits_500 = result.loc[result["radius_m"] == 500, "n_hits"].iloc[0]
    assert hits_500 >= hits_0


def test_buffer_sensitivity_baseline_delta_zero():
    """Delta at radius=0 is always 0 (it's the baseline)."""
    deposits, zones = _make_point_deposits_and_zones()
    result = buffer_sensitivity(deposits, zones, radii_m=(0, 250, 500))
    assert result.loc[result["radius_m"] == 0, "delta_pct"].iloc[0] == 0.0


def test_buffer_sensitivity_geographic_crs_raises():
    """Geographic CRS (degrees) should raise ValueError — buffer is meaningless in degrees."""
    import geopandas as gpd
    from shapely.geometry import Point, box as shapely_box

    zones = gpd.GeoDataFrame(geometry=[shapely_box(0, 0, 1, 1)], crs="EPSG:4326")
    deposits = gpd.GeoDataFrame(geometry=[Point(0.5, 0.5)], crs="EPSG:4326")
    try:
        buffer_sensitivity(deposits, zones, radii_m=(0, 500))
        assert False, "Should have raised ValueError"
    except ValueError:
        pass
