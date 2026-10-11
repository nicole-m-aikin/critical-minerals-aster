"""GEE parity module — no live Earth Engine, no local ASTER rasters required."""

from __future__ import annotations

import sys
from datetime import timezone
from pathlib import Path

import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.config import (  # noqa: E402
    ClassificationParams,
    SiteConfig,
    load_site_by_id,
)
from critical_minerals_aster.gee_parity import (  # noqa: E402
    DEFAULT_PARITY_SITES,
    GEE_ASTER_L1T_COLLECTION,
    SITE_ROLES,
    TIR_UCC,
    TIRWorkflowProduct,
    GeeRecipe,
    cohen_kappa,
    compare_classified,
    dn_to_radiance,
    gee_image_id_from_granule,
    gee_recipe_javascript,
    histogram_percentile,
    jaccard_binary,
    overall_accuracy,
    parse_aster_l1t_granule_id,
    probe_clip_order,
    probe_histogram_percentile,
    recipe_for_site,
    summarize_probe_effect,
)
from critical_minerals_aster.spectral import alteration_ratios  # noqa: E402


# ---------------------------------------------------------------------------
# Granule ID mapping
# ---------------------------------------------------------------------------


def test_parse_goldfield_granule_id():
    """AST_L1T_00402072010054840 → 2010-02-07 05:48:40 UTC."""
    gid = "AST_L1T_00402072010054840_20250702175723"
    parsed = parse_aster_l1t_granule_id(gid)
    assert parsed is not None
    assert parsed.version == "004"
    assert parsed.acquired.year == 2010
    assert parsed.acquired.month == 2
    assert parsed.acquired.day == 7
    assert parsed.acquired.hour == 5
    assert parsed.acquired.minute == 48
    assert parsed.acquired.second == 40
    assert parsed.acquired.tzinfo == timezone.utc
    assert parsed.gee_image_id == f"{GEE_ASTER_L1T_COLLECTION}/20100207T054840"


def test_parse_mcdermitt_granule_id():
    gid = "AST_L1T_00407232010184946_20250705074029"
    parsed = parse_aster_l1t_granule_id(gid)
    assert parsed is not None
    assert parsed.acquired.month == 7
    assert parsed.acquired.day == 23
    assert gee_image_id_from_granule(gid) == f"{GEE_ASTER_L1T_COLLECTION}/20100723T184946"


def test_mosaic_id_is_not_a_granule():
    assert parse_aster_l1t_granule_id("goldfield_mosaic") is None
    assert gee_image_id_from_granule("goldfield_mosaic") is None


# ---------------------------------------------------------------------------
# Recipe / JavaScript (the GEE port, inspectable without ee)
# ---------------------------------------------------------------------------


def test_recipe_for_site_pins_goldfield_acquisition():
    site = load_site_by_id("goldfield", Path(__file__).resolve().parents[1] / "sites")
    recipe = recipe_for_site(site)
    assert recipe.site_id == "goldfield"
    assert recipe.collection == GEE_ASTER_L1T_COLLECTION
    assert recipe.bands == ("B10", "B11", "B12", "B13", "B14")
    assert recipe.low_pct == 70
    assert recipe.high_pct == 90
    assert recipe.pin_image_id == f"{GEE_ASTER_L1T_COLLECTION}/20150714T184607"
    assert recipe.bbox_wgs84 == tuple(site.bbox_wgs84)


def test_javascript_recipe_contains_the_six_steps():
    recipe = GeeRecipe(
        site_id="goldfield",
        bbox_wgs84=(-117.65, 37.55, -117.05, 38.05),
        start="2010-01-01",
        end="2023-12-31",
        pin_image_id=f"{GEE_ASTER_L1T_COLLECTION}/20100207T054840",
        crs="EPSG:32611",
    )
    js = gee_recipe_javascript(recipe)
    assert "ASTER/AST_L1T_003" in js
    assert "B10" in js and "B14" in js
    assert "clip(region)" in js
    assert "B13').divide(image.select('B14')" in js
    assert "B13').divide(image.select('B12')" in js
    assert "B12').divide(image.select('B13')" in js
    assert "Reducer.percentile([70.0, 90.0])" in js
    assert "crs: 'EPSG:32611'" in js
    assert recipe.pin_image_id in js


def test_javascript_flags_default_crs_resampling():
    recipe = GeeRecipe(
        site_id="bisbee",
        bbox_wgs84=(-110.15, 31.25, -109.65, 31.70),
        start="2010-01-01",
        end="2023-12-31",
        crs=None,
    )
    js = gee_recipe_javascript(recipe)
    assert "EPSG:4326" in js
    assert "resampling" in js.lower()
    assert "crs:" not in js.split("reduceRegion")[1].split("maxPixels")[0] or "omitted" in js


def test_default_panel_sites_exist_and_are_documented():
    sites_dir = Path(__file__).resolve().parents[1] / "sites"
    assert 5 <= len(DEFAULT_PARITY_SITES) <= 10
    for site_id in DEFAULT_PARITY_SITES:
        site = load_site_by_id(site_id, sites_dir)
        assert site.id == site_id
        assert len(site.bbox_wgs84) == 4
        lon0, lat0, lon1, lat1 = site.bbox_wgs84
        assert lon0 < lon1 and lat0 < lat1
        assert site_id in SITE_ROLES


# ---------------------------------------------------------------------------
# Agreement metrics
# ---------------------------------------------------------------------------


def test_identical_maps_are_perfect_agreement():
    score = np.array([[0, 1, 3], [2, 4, 6]], dtype=np.uint8)
    valid = np.ones_like(score, dtype=bool)
    stats = compare_classified(score, score.copy(), valid)
    assert stats.n_overlap == 6
    assert stats.overall_accuracy == pytest.approx(1.0)
    assert stats.kappa == pytest.approx(1.0)
    assert stats.jaccard_strong == pytest.approx(1.0)
    assert stats.score_mae == pytest.approx(0.0)
    assert stats.strong_frac_delta == pytest.approx(0.0)


def test_complementary_strong_zones_have_zero_jaccard():
    ref = np.array([[3, 3, 0], [0, 0, 0]], dtype=np.uint8)
    other = np.array([[0, 0, 0], [3, 3, 0]], dtype=np.uint8)
    valid = np.ones_like(ref, dtype=bool)
    stats = compare_classified(ref, other, valid)
    assert stats.jaccard_strong == pytest.approx(0.0)
    assert stats.strong_frac_ref == pytest.approx(2 / 6)
    assert stats.strong_frac_other == pytest.approx(2 / 6)


def test_both_empty_strong_zones_jaccard_is_one():
    a = np.zeros((3, 3), dtype=bool)
    assert jaccard_binary(a, a) == 1.0


def test_kappa_chance_agreement_is_near_zero():
    # 50/50 labels, independent — κ should sit near 0, not 0.5.
    rng = np.random.default_rng(0)
    a = rng.integers(0, 2, size=10_000)
    b = rng.integers(0, 2, size=10_000)
    k = cohen_kappa(a, b, labels=(0, 1))
    assert abs(k) < 0.05


def test_compare_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="shape"):
        compare_classified(np.zeros((2, 2)), np.zeros((3, 3)), np.ones((2, 2), dtype=bool))


def test_empty_overlap_returns_nans():
    score = np.ones((2, 2), dtype=np.uint8)
    valid = np.zeros((2, 2), dtype=bool)
    stats = compare_classified(score, score, valid)
    assert stats.n_overlap == 0
    assert np.isnan(stats.jaccard_strong)


def test_threshold_deltas_are_other_minus_ref():
    ref = np.array([[0, 3], [3, 0]], dtype=np.uint8)
    other = np.array([[0, 3], [3, 0]], dtype=np.uint8)
    valid = np.ones((2, 2), dtype=bool)
    stats = compare_classified(
        ref,
        other,
        valid,
        ref_thresholds={"silica": (1.00, 1.10)},
        other_thresholds={"silica": (1.02, 1.15)},
    )
    assert stats.p70_delta["silica"] == pytest.approx(0.02)
    assert stats.p90_delta["silica"] == pytest.approx(0.05)


# ---------------------------------------------------------------------------
# Radiance UCC — DN ratios ≠ radiance ratios
# ---------------------------------------------------------------------------


def test_ucc_changes_silica_ratio_because_coefficients_differ():
    # Constant DN → silica_DN = 1, silica_rad = UCC13/UCC14 ≠ 1.
    b12 = np.full((4, 4), 500.0)
    b13 = np.full((4, 4), 500.0)
    b14 = np.full((4, 4), 500.0)
    silica_dn, _, _ = alteration_ratios(b12, b13, b14)
    silica_rad, _, _ = alteration_ratios(
        dn_to_radiance(b12, "B12"),
        dn_to_radiance(b13, "B13"),
        dn_to_radiance(b14, "B14"),
    )
    assert np.allclose(silica_dn, 1.0)
    expected = TIR_UCC["B13"] / TIR_UCC["B14"]
    assert np.allclose(silica_rad, expected)
    assert not np.allclose(silica_dn, silica_rad)


def test_dn_zero_stays_nodata_after_ucc():
    dn = np.array([[0.0, 500.0], [np.nan, 400.0]])
    rad = dn_to_radiance(dn, "B13")
    assert np.isnan(rad[0, 0])
    assert np.isnan(rad[1, 0])
    assert rad[0, 1] == pytest.approx((500.0 - 1.0) * TIR_UCC["B13"])


# ---------------------------------------------------------------------------
# Histogram percentile vs exact
# ---------------------------------------------------------------------------


def test_histogram_percentile_matches_exact_on_uniform_sample():
    rng = np.random.default_rng(1)
    x = rng.uniform(0.8, 1.2, size=50_000)
    exact = float(np.nanpercentile(x, 70))
    hist = histogram_percentile(x, 70, n_bins=10_000)
    assert hist == pytest.approx(exact, abs=5e-4)


def test_histogram_percentile_empty_is_nan():
    assert np.isnan(histogram_percentile(np.array([np.nan, np.nan]), 70))


# ---------------------------------------------------------------------------
# Clip-order probe — the project's own known bug class
# ---------------------------------------------------------------------------


def _toy_local_product() -> tuple[SiteConfig, TIRWorkflowProduct]:
    """Clipped 2×2 box sitting inside a 2×6 unclipped scene.

    In-box pixels are all ≤ 0.4.  The extra four columns are 10.0, so
    eight of twelve scene pixels sit at the spike.  The 70th and 90th
    of the full scene are both ~10; applying those cuts to the box
    classifies every in-box pixel as background.
    """
    site = SiteConfig(
        id="toy",
        name="Toy",
        bbox_wgs84=(-117.0, 37.0, -116.5, 37.5),
        classification=ClassificationParams(low_pct=70, high_pct=90),
    )
    # Clipped grid (what the real pipeline analyses).
    silica = np.array([[0.10, 0.20], [0.30, 0.40]])
    carbonate = silica.copy()
    mafic = silica.copy()
    # Unclipped scene: same 2×2 plus four extra columns at 10.0.
    unclipped = np.array(
        [
            [0.10, 0.20, 10.0, 10.0, 10.0, 10.0],
            [0.30, 0.40, 10.0, 10.0, 10.0, 10.0],
        ]
    )
    from critical_minerals_aster.classification import classify_percentiles, combined_score

    s_cls, s_lo, s_hi = classify_percentiles(silica, 70, 90)
    c_cls, c_lo, c_hi = classify_percentiles(carbonate, 70, 90)
    m_cls, m_lo, m_hi = classify_percentiles(mafic, 70, 90)
    local = TIRWorkflowProduct(
        site_id="toy",
        source="local",
        granule_id=None,
        gee_image_id=None,
        crs="EPSG:32611",
        silica=silica,
        carbonate=carbonate,
        mafic=mafic,
        silica_cls=s_cls,
        carbonate_cls=c_cls,
        mafic_cls=m_cls,
        combined=combined_score(s_cls, c_cls, m_cls),
        valid=np.ones((2, 2), dtype=bool),
        thresholds={
            "silica": (s_lo, s_hi),
            "carbonate": (c_lo, c_hi),
            "mafic": (m_lo, m_hi),
        },
        unclipped_ratios={
            "silica": unclipped,
            "carbonate": unclipped.copy(),
            "mafic": unclipped.copy(),
        },
    )
    return site, local


def test_clip_order_probe_disagrees_when_scene_has_out_of_box_spike():
    site, local = _toy_local_product()
    probed = probe_clip_order(local, site)
    assert probed is not None
    # Local (clip first) must flag at least one strong pixel — the 70/90
    # of {0.10, 0.20, 0.30, 0.40} sit inside that set.
    assert int((local.combined >= 3).sum()) >= 1
    # Unclipped 90th is near 10, so every in-box pixel is background.
    assert int((probed.combined >= 3).sum()) == 0
    stats = compare_classified(
        local.combined,
        probed.combined,
        local.valid,
        ref_thresholds=local.thresholds,
        other_thresholds=probed.thresholds,
    )
    assert stats.jaccard_strong == pytest.approx(0.0)
    assert stats.p90_delta["silica"] > 1.0


def test_histogram_probe_same_grid_as_local():
    site, local = _toy_local_product()
    probed = probe_histogram_percentile(local, site)
    assert probed.combined.shape == local.combined.shape
    assert probed.source == "probe:histogram_percentile"


def test_summarize_probe_effect_averages_by_probe():
    from critical_minerals_aster.gee_parity import ParityRow

    rows = [
        ParityRow(
            site_id="a",
            site_name="A",
            role="",
            probe="probe:clip_order",
            n_overlap=10,
            overall_accuracy=0.5,
            kappa=0.2,
            jaccard_strong=0.1,
            score_mae=1.0,
            strong_frac_local=0.2,
            strong_frac_other=0.4,
            strong_frac_delta=0.2,
            p70_delta_silica=None,
            p90_delta_silica=None,
            local_granule_id=None,
            gee_image_id=None,
            local_crs="",
            notes="",
        ),
        ParityRow(
            site_id="b",
            site_name="B",
            role="",
            probe="probe:clip_order",
            n_overlap=10,
            overall_accuracy=0.6,
            kappa=0.4,
            jaccard_strong=0.3,
            score_mae=0.5,
            strong_frac_local=0.2,
            strong_frac_other=0.1,
            strong_frac_delta=-0.1,
            p70_delta_silica=None,
            p90_delta_silica=None,
            local_granule_id=None,
            gee_image_id=None,
            local_crs="",
            notes="",
        ),
    ]
    summary = summarize_probe_effect(rows)
    assert summary["probe:clip_order"]["n_sites"] == 2
    assert summary["probe:clip_order"]["mean_jaccard_strong"] == pytest.approx(0.2)
    assert summary["probe:clip_order"]["mean_abs_strong_frac_delta"] == pytest.approx(0.15)


def test_gee_extra_is_not_imported_by_the_module():
    """CI and the offline benchmark must not require earthengine-api."""
    import critical_minerals_aster.gee_parity as mod

    assert "ee" not in sys.modules or not hasattr(mod, "ee")
    # The live entry points exist but do not run at import.
    assert hasattr(mod, "run_gee_workflow")
    assert hasattr(mod, "build_gee_image")
    assert hasattr(mod, "initialize_earth_engine")


def test_overall_accuracy_helper():
    a = np.array([0, 1, 2, 2])
    b = np.array([0, 1, 1, 2])
    assert overall_accuracy(a, b) == pytest.approx(0.75)
