"""Tests for lithology.py — stratified expectation, host-restricted test (offline, synthetic)."""

import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from scipy.stats import binomtest
from shapely.geometry import Point, box
from shapely.ops import unary_union

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.lithology import (  # noqa: E402
    UNMAPPED,
    assign_strata,
    build_strata,
    generalize_lithology,
    host_restricted_test,
    poisson_binomial_sf,
    stratified_test,
    stratum_table,
)

CRS = "EPSG:32612"


def _synthetic_site():
    """100×100 footprint: west half carbonate, east half clastic.

    Zones: 50×20 strip inside carbonate (20% of carbonate, 10% of footprint)
    plus nothing in clastic. Deposits: 10 in carbonate (4 inside zones),
    5 in clastic (none inside zones).
    """
    footprint = box(0, 0, 100, 100)
    units = gpd.GeoDataFrame(
        {"lith_class": ["carbonate_chemical", "clastic"]},
        geometry=[box(0, 0, 50, 100), box(50, 0, 100, 100)],
        crs=CRS,
    )
    zone_union = box(0, 0, 50, 20)
    pts = (
        [Point(10 + i * 10, 10) for i in range(4)]          # carbonate, in zone
        + [Point(5 + i * 7, 60) for i in range(6)]          # carbonate, out of zone
        + [Point(60 + i * 7, 50) for i in range(5)]         # clastic, out of zone
    )
    deposits = gpd.GeoDataFrame(geometry=pts, crs=CRS)
    deposits["hit"] = deposits.geometry.within(zone_union)
    return footprint, units, zone_union, deposits


def test_build_strata_partitions_footprint():
    footprint, units, _, _ = _synthetic_site()
    strata = build_strata(units, footprint)
    assert set(strata["lith_class"]) == {"carbonate_chemical", "clastic"}
    assert strata.geometry.area.sum() == pytest.approx(footprint.area)


def test_build_strata_unmapped_and_overlap():
    footprint = box(0, 0, 100, 100)
    # Overlapping units (state-line mismatch) and a gap on the east edge.
    units = gpd.GeoDataFrame(
        {"lith_class": ["a", "b"]},
        geometry=[box(0, 0, 60, 100), box(40, 0, 80, 100)],
        crs=CRS,
    )
    strata = build_strata(units, footprint).set_index("lith_class")
    assert strata.loc["a"].geometry.area == pytest.approx(6000)
    assert strata.loc["b"].geometry.area == pytest.approx(2000)  # overlap assigned to "a"
    assert strata.loc[UNMAPPED].geometry.area == pytest.approx(2000)
    assert strata.geometry.area.sum() == pytest.approx(footprint.area)


def test_stratum_table_and_stratified_expectation():
    footprint, units, zone_union, deposits = _synthetic_site()
    strata = build_strata(units, footprint)
    labels = assign_strata(deposits, strata)
    tab = stratum_table(strata, zone_union, labels, deposits["hit"]).set_index("lith_class")

    assert tab.loc["carbonate_chemical", "p0"] == pytest.approx(0.2)
    assert tab.loc["clastic", "p0"] == pytest.approx(0.0)
    assert tab.loc["carbonate_chemical", "n"] == 10
    assert tab.loc["carbonate_chemical", "hits"] == 4
    assert tab.loc["clastic", "n"] == 5

    res = stratified_test(tab.reset_index())
    # Expected = 10·0.2 + 5·0.0 = 2.0 vs unconditioned 15·0.1 = 1.5
    assert res["expected"] == pytest.approx(2.0)
    assert res["enrichment"] == pytest.approx(2.0)
    assert res["n"] == 15 and res["hits"] == 4
    # Only the carbonate stratum is random → equals Binomial(10, 0.2) tail.
    assert res["p"] == pytest.approx(binomtest(4, 10, 0.2, alternative="greater").pvalue)


def test_stratified_equals_binomial_when_p0_uniform():
    tab = pd.DataFrame({
        "lith_class": ["a", "b", "c"],
        "area": [1.0, 2.0, 3.0],
        "zone_area": [0.1, 0.2, 0.3],
        "n": [20, 30, 50],
        "hits": [5, 6, 9],
    })
    tab["p0"] = tab["zone_area"] / tab["area"]
    res = stratified_test(tab)
    assert res["expected"] == pytest.approx(10.0)
    assert res["p"] == pytest.approx(binomtest(20, 100, 0.1, alternative="greater").pvalue)


def test_poisson_binomial_sf_matches_monte_carlo():
    n_per = np.array([12, 7, 30])
    p_per = np.array([0.3, 0.05, 0.12])
    rng = np.random.default_rng(42)
    sims = sum(rng.binomial(n, p, size=200_000) for n, p in zip(n_per, p_per))
    for k in (0, 5, 9, 14):
        assert poisson_binomial_sf(n_per, p_per, k) == pytest.approx((sims >= k).mean(), abs=4e-3)
    assert poisson_binomial_sf(n_per, p_per, 0) == 1.0
    assert poisson_binomial_sf(n_per, p_per, 50) == 0.0


def test_host_restricted_drops_barren_strata():
    footprint = box(0, 0, 100, 100)
    units = gpd.GeoDataFrame(
        {"lith_class": ["host", "barren"]},
        geometry=[box(0, 0, 25, 100), box(25, 0, 100, 100)],
        crs=CRS,
    )
    # Zones: 500 m² inside host (20% of host), 500 m² inside barren.
    zone_union = unary_union([box(0, 0, 25, 20), box(50, 0, 75, 20)])
    deposits = gpd.GeoDataFrame(geometry=[Point(5 + i * 2, 10) for i in range(3)]
                                + [Point(5 + i * 2, 50) for i in range(7)], crs=CRS)
    deposits["hit"] = deposits.geometry.within(zone_union)
    strata = build_strata(units, footprint)
    tab = stratum_table(strata, zone_union, assign_strata(deposits, strata), deposits["hit"])

    res = host_restricted_test(tab)
    assert res["host_classes"] == ["host"]
    assert res["p0"] == pytest.approx(0.2)            # vs unconditioned 1000/10000 = 0.1
    assert res["host_area_frac"] == pytest.approx(0.25)
    assert res["n"] == 10 and res["hits"] == 3
    assert res["enrichment"] == pytest.approx(1.5)
    assert res["p"] == pytest.approx(binomtest(3, 10, 0.2, alternative="greater").pvalue)


def test_generalize_lithology_mapping():
    assert generalize_lithology("Sedimentary, carbonate") == "carbonate_chemical"
    assert generalize_lithology("Metamorphic, carbonate") == "carbonate_chemical"
    assert generalize_lithology("Igneous, intrusive") == "intrusive"
    assert generalize_lithology("Metamorphic, sedimentary clastic") == "clastic"
    assert generalize_lithology("Sedimentary, undifferentiated") == "sedimentary_mixed"
    assert generalize_lithology(" Metamorphic, gneiss ") == "metamorphic"
    assert generalize_lithology("Unconsolidated, undifferentiated") == "unconsolidated"
    assert generalize_lithology(None) == "other"
