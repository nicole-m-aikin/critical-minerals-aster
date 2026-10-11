"""Tests for the bedrock / cover / unmapped stratified null (offline, synthetic geometries)."""

import importlib.util
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from scipy.stats import binom, binomtest, mannwhitneyu
from shapely.geometry import Point, box
from shapely.ops import unary_union

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster import lithology  # noqa: E402
from critical_minerals_aster.lithology import (  # noqa: E402
    BEDROCK,
    COVER,
    STRATUM_CODES,
    UNMAPPED,
    assign_cover_stratum,
    cover_strata,
    cover_strata_from_units,
    rasterize_cover_strata,
    stratified_null,
    stratum_zone_fractions,
)
from critical_minerals_aster.significance import (  # noqa: E402
    _bbox_in_crs,
    run_permutation_stratified,
    stratified_mannwhitney,
)

CRS = "EPSG:32612"
SEED = 42


def _units(polys_and_liths, crs=CRS):
    return gpd.GeoDataFrame(
        {"GENERALIZED_LITH": [l for _, l in polys_and_liths]},
        geometry=[p for p, _ in polys_and_liths],
        crs=crs,
    )


def _synthetic_site():
    """100×100 footprint: x<40 bedrock, 40≤x<80 cover (alluvium), x≥80 unmapped.

    Zones: 40×50 block in bedrock (50% of bedrock, 20% of footprint).
    Deposits: 6 in bedrock (3 in zone), 4 in cover (0 in zone).
    """
    footprint = box(0, 0, 100, 100)
    units = _units([
        (box(0, 0, 40, 100), "Igneous, intrusive"),
        (box(40, 0, 80, 100), "Unconsolidated, undifferentiated"),
    ])
    zone_union = box(0, 0, 40, 50)
    pts = (
        [Point(5 + i * 10, 25) for i in range(3)]      # bedrock, in zone
        + [Point(5 + i * 10, 75) for i in range(3)]    # bedrock, out of zone
        + [Point(45 + i * 8, 50) for i in range(4)]    # cover, out of zone
    )
    deposits = gpd.GeoDataFrame(geometry=pts, crs=CRS)
    deposits["hit"] = deposits.geometry.within(zone_union)
    return footprint, units, zone_union, deposits


# ---------------------------------------------------------------------------
# Strata construction
# ---------------------------------------------------------------------------


def test_cover_strata_partition_footprint():
    footprint, units, _, _ = _synthetic_site()
    s = cover_strata_from_units(units, footprint)
    assert s[BEDROCK].area == pytest.approx(4000)
    assert s[COVER].area == pytest.approx(4000)
    assert s[UNMAPPED].area == pytest.approx(2000)
    assert sum(g.area for g in s.values()) == pytest.approx(footprint.area)
    assert s[BEDROCK].intersection(s[COVER]).area == pytest.approx(0)


def test_cover_strata_overlap_cover_wins_and_other_is_bedrock():
    footprint = box(0, 0, 100, 100)
    units = _units([
        (box(0, 0, 60, 100), "Sedimentary, carbonate"),
        (box(50, 0, 90, 100), "Unconsolidated and Sedimentary, undifferentiated"),
        (box(90, 0, 100, 100), "Water"),
    ])
    s = cover_strata_from_units(units, footprint)
    assert s[COVER].area == pytest.approx(4000)      # overlap 50–60 assigned to cover
    assert s[BEDROCK].area == pytest.approx(6000)    # carbonate 0–50 + "other" (water) 90–100
    assert s[UNMAPPED].area == pytest.approx(0)


def test_cover_strata_no_sgmc_is_all_unmapped():
    footprint = box(0, 0, 10, 10)
    empty = gpd.GeoDataFrame({"GENERALIZED_LITH": []}, geometry=[], crs=CRS)
    s = cover_strata_from_units(empty, footprint)
    assert s[UNMAPPED].area == pytest.approx(100)
    assert s[BEDROCK].is_empty and s[COVER].is_empty


def test_cover_strata_reads_cache_without_network(tmp_path, monkeypatch):
    footprint, units, _, _ = _synthetic_site()
    cache = tmp_path / "data" / "lithology" / "toy_sgmc_units.geojson"
    cache.parent.mkdir(parents=True)
    units.to_crs("EPSG:4326").to_file(cache, driver="GeoJSON")

    def _no_network(*_a, **_k):
        raise AssertionError("network fetch attempted despite cache")

    monkeypatch.setattr(lithology, "fetch_sgmc_units", _no_network)
    s = cover_strata("toy", footprint, CRS, repo_root=tmp_path)
    assert s[BEDROCK].area == pytest.approx(4000, rel=1e-3)
    assert s[COVER].area == pytest.approx(4000, rel=1e-3)
    assert s[UNMAPPED].area == pytest.approx(2000, rel=1e-3)


def test_cover_strata_fetches_and_caches_when_missing(tmp_path, monkeypatch):
    footprint = box(500_000, 3_500_000, 501_000, 3_501_000)
    cover_poly = gpd.GeoSeries([box(500_000, 3_500_000, 500_500, 3_501_000)], crs=CRS).to_crs("EPSG:4326")
    calls = []

    def _fake_fetch(bbox_wgs84, **_k):
        calls.append(bbox_wgs84)
        return {"type": "FeatureCollection", "features": [{
            "type": "Feature",
            "properties": {"STATE": "AZ", "GENERALIZED_LITH": "Unconsolidated, undifferentiated"},
            "geometry": cover_poly.iloc[0].__geo_interface__,
        }]}

    monkeypatch.setattr(lithology, "fetch_sgmc_units", _fake_fetch)
    s = cover_strata("fresh", footprint, CRS, repo_root=tmp_path)
    assert len(calls) == 1
    assert (tmp_path / "data" / "lithology" / "fresh_sgmc_units.geojson").is_file()
    assert s[COVER].area == pytest.approx(500_000, rel=1e-3)
    assert s[UNMAPPED].area == pytest.approx(500_000, rel=1e-3)
    cover_strata("fresh", footprint, CRS, repo_root=tmp_path)
    assert len(calls) == 1  # second call served from cache


def test_assign_cover_stratum_snaps_off_footprint_points_to_nearest():
    footprint, units, _, _ = _synthetic_site()
    s = cover_strata_from_units(units, footprint)
    pts = gpd.GeoDataFrame(geometry=[Point(10, 10), Point(60, 10), Point(90, 10),
                                     Point(-3, 50), Point(60, 101)], crs=CRS)
    labels = assign_cover_stratum(pts, s)
    assert labels.tolist() == [BEDROCK, COVER, UNMAPPED, BEDROCK, COVER]


# ---------------------------------------------------------------------------
# Stratified expectation and p-value
# ---------------------------------------------------------------------------


def test_stratified_expectation_and_p():
    footprint, units, zone_union, deposits = _synthetic_site()
    s = cover_strata_from_units(units, footprint)
    frac = stratum_zone_fractions(s, zone_union)
    assert frac.loc[BEDROCK, "p"] == pytest.approx(0.5)
    assert frac.loc[COVER, "p"] == pytest.approx(0.0)
    assert frac.loc[UNMAPPED, "p"] == pytest.approx(0.0)

    p_i = assign_cover_stratum(deposits, s).map(frac["p"]).to_numpy(dtype=float)
    hits = int(deposits["hit"].sum())
    res = stratified_null(p_i, hits)
    # E = 6·0.5 + 4·0 = 3 vs geometric 10·0.2 = 2 — cover deposits cannot be hit.
    assert res["expected"] == pytest.approx(3.0)
    assert res["p0_eff"] == pytest.approx(0.3)
    assert res["enrichment"] == pytest.approx(1.0)
    assert res["p"] == pytest.approx(binomtest(3, 6, 0.5, alternative="greater").pvalue)
    geo = binomtest(3, 10, 0.2, alternative="greater").pvalue
    assert res["p"] > geo  # the geometric null overstated significance here


def test_stratified_null_reduces_to_binomial_for_uniform_p():
    p_i = np.full(40, 0.15)
    res = stratified_null(p_i, 11)
    assert res["expected"] == pytest.approx(6.0)
    assert res["p"] == pytest.approx(binomtest(11, 40, 0.15, alternative="greater").pvalue)


def test_stratified_null_empty_and_zero_p():
    assert stratified_null(np.array([]), 0)["p"] == 1.0
    res = stratified_null(np.zeros(5), 0)
    assert res["p"] == 1.0 and np.isnan(res["enrichment"])
    assert stratified_null(np.zeros(5), 1)["p"] == 0.0


def test_rasterize_cover_strata_codes():
    from rasterio.transform import from_bounds

    footprint, units, _, _ = _synthetic_site()
    s = cover_strata_from_units(units, footprint)
    codes = rasterize_cover_strata(s, from_bounds(0, 0, 100, 100, 10, 10), (10, 10))
    assert (codes[:, :4] == STRATUM_CODES[BEDROCK]).all()
    assert (codes[:, 4:8] == STRATUM_CODES[COVER]).all()
    assert (codes[:, 8:] == STRATUM_CODES[UNMAPPED]).all()


# ---------------------------------------------------------------------------
# Stratified Monte Carlo and Mann–Whitney
# ---------------------------------------------------------------------------


def test_stratified_permutation_matches_poisson_binomial():
    bbox = (-111.0, 32.0, -110.9, 32.1)
    minx, miny, maxx, maxy = _bbox_in_crs(bbox, CRS)
    midx, midy = (minx + maxx) / 2, (miny + maxy) / 2
    strata = {
        BEDROCK: box(minx, miny, midx, maxy),
        COVER: box(midx, miny, maxx, maxy),
        UNMAPPED: box(0, 0, 0, 0).buffer(0),
    }
    # Zone covers half of bedrock and none of cover.
    zones = gpd.GeoDataFrame(geometry=[box(minx, miny, midx, midy)], crs=CRS)
    p_mc = run_permutation_stratified(zones, bbox, strata, {BEDROCK: 10, COVER: 5}, 7,
                                      n_iter=10_000, seed=SEED)
    exact = float(binom.sf(6, 10, 0.5))  # cover deposits never hit
    assert p_mc == pytest.approx(exact, abs=0.02)

    # Geometric null would have used p = 0.25 for all 15 → much smaller p.
    assert float(binom.sf(6, 15, 0.25)) < exact / 2


def test_stratified_mwu_single_stratum_matches_scipy():
    rng = np.random.default_rng(SEED)
    dep = rng.integers(0, 7, 40).astype(float) + 0.5
    bg = rng.integers(0, 7, 2000).astype(float)
    res = stratified_mannwhitney(dep, np.ones(40), bg, np.ones(2000))
    ref = mannwhitneyu(dep, bg, alternative="greater", method="asymptotic", use_continuity=False)
    assert res["auc"] == pytest.approx(ref.statistic / (40 * 2000))
    assert res["p"] == pytest.approx(ref.pvalue, rel=1e-6)


def test_stratified_mwu_removes_bedrock_vs_cover_confound():
    """Deposits all on bedrock; bedrock scores high, cover scores low.

    Within bedrock deposits are indistinguishable from background, so the
    pooled test is spuriously significant and the stratified one is not.
    """
    rng = np.random.default_rng(SEED)
    bedrock_bg = rng.integers(2, 7, 5000).astype(float)
    cover_bg = rng.integers(0, 3, 5000).astype(float)
    dep = rng.integers(2, 7, 60).astype(float)

    pooled = mannwhitneyu(dep, np.concatenate([bedrock_bg, cover_bg]), alternative="greater")
    assert pooled.pvalue < 1e-6

    res = stratified_mannwhitney(
        dep, np.full(60, 1),
        np.concatenate([bedrock_bg, cover_bg]),
        np.concatenate([np.full(5000, 1), np.full(5000, 2)]),
    )
    assert res["n_used"] == 60
    assert abs(res["auc"] - 0.5) < 0.1
    assert res["p"] > 0.05


# ---------------------------------------------------------------------------
# Phase 5: declustered E uses the medoid's own p_i
# ---------------------------------------------------------------------------


def _load_phase5():
    path = Path(__file__).resolve().parents[1] / "scripts" / "phase5_clustering_sensitivity.py"
    spec = importlib.util.spec_from_file_location("phase5_clustering_sensitivity", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_phase5_declustered_expectation_uses_unit_p_i():
    phase5 = _load_phase5()
    # A tight 5-point bedrock camp (p=0.5) and 3 isolated cover points (p=0.05).
    camp = [Point(1000 + i, 1000) for i in range(5)]
    cover = [Point(5000 + 3000 * i, 9000) for i in range(3)]
    deposits = gpd.GeoDataFrame(
        {"in_zone": [True, True, True, False, False, False, False, False],
         "p_null": [0.5] * 5 + [0.05] * 3},
        geometry=camp + cover, crs=CRS,
    )
    raw = phase5.effective_units(deposits, 0)
    assert len(raw) == 8
    assert stratified_null(raw["p_null"].to_numpy(), int(raw["in_zone"].sum()))["expected"] == pytest.approx(2.65)

    units = phase5.effective_units(deposits, 250)
    assert len(units) == 4                      # camp → 1 medoid unit + 3 noise singletons
    assert int(units["in_zone"].sum()) == 1     # medoid is the middle camp point (in zone)
    assert sorted(units["p_null"]) == [0.05, 0.05, 0.05, 0.5]
    res = stratified_null(units["p_null"].to_numpy(), int(units["in_zone"].sum()))
    assert res["expected"] == pytest.approx(0.5 + 3 * 0.05)
    # The site-level p0 (2.65/8) would have given E = 4 · 0.33 = 1.33 instead.
    assert res["p"] == pytest.approx(
        lithology.poisson_binomial_sf(np.array([1, 3]), np.array([0.5, 0.05]), 1)
    )
    assert phase5.effective_n_hits(deposits, 250) == (4, 1)
