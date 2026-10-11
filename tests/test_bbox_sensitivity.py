"""Tests for bbox_sensitivity.py — box variants, reclassification, critical hit test (offline)."""

import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from shapely.geometry import box

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.bbox_sensitivity import (  # noqa: E402
    bbox_overlap_fraction,
    centred_box,
    classify_strong_zones,
    critical_hit_test,
    footprint_centroid_wgs84,
    perturb_bbox,
    perturbation_variants,
    raster_extent_polygon,
    restrict_to_region,
    summarize_hits,
)
from critical_minerals_aster.config import ClassificationParams  # noqa: E402

UTM = rasterio.crs.CRS.from_epsg(32612)
# 100 x 100 grid of 90 m pixels near (-111.0, 32.0)
TRANSFORM = rasterio.Affine(90.0, 0.0, 500_000.0, 0.0, -90.0, 3_540_000.0)


# ---------------------------------------------------------------------------
# Box geometry
# ---------------------------------------------------------------------------


def test_perturb_bbox_shift_preserves_size():
    b = (-110.0, 31.0, -109.5, 31.4)
    out = perturb_bbox(b, dx_frac=0.2, dy_frac=-0.1)
    assert out[2] - out[0] == pytest.approx(0.5)
    assert out[3] - out[1] == pytest.approx(0.4)
    assert out[0] == pytest.approx(-110.0 + 0.1)
    assert out[1] == pytest.approx(31.0 - 0.04)


def test_perturb_bbox_scale_about_centre():
    b = (0.0, 0.0, 1.0, 2.0)
    out = perturb_bbox(b, scale=0.5)
    assert out == pytest.approx((0.25, 0.5, 0.75, 1.5))


def test_perturbation_variants_eight_distinct_boxes():
    b = (-110.0, 31.0, -109.5, 31.4)
    v = perturbation_variants(b, frac=0.15)
    assert len(v) == 8
    assert len({tuple(round(x, 9) for x in bb) for bb in v.values()}) == 8
    assert "shrink_15" in v and "expand_15" in v
    shrink = v["shrink_15"]
    assert shrink[2] - shrink[0] == pytest.approx(0.5 * 0.85)


def test_centred_box():
    assert centred_box(-110.0, 32.0, 0.46) == pytest.approx((-110.23, 31.77, -109.77, 32.23))


def test_bbox_overlap_fraction():
    a = (0.0, 0.0, 2.0, 2.0)
    assert bbox_overlap_fraction(a, (1.0, 1.0, 3.0, 3.0)) == pytest.approx(0.25)
    assert bbox_overlap_fraction(a, (5.0, 5.0, 6.0, 6.0)) == 0.0
    assert bbox_overlap_fraction(a, (-1.0, -1.0, 3.0, 3.0)) == pytest.approx(1.0)


def test_raster_extent_polygon_matches_grid():
    poly = raster_extent_polygon(TRANSFORM, (100, 50))
    assert poly.bounds == pytest.approx((500_000.0, 3_531_000.0, 504_500.0, 3_540_000.0))


def test_footprint_centroid_wgs84_roundtrip():
    fp = gpd.GeoDataFrame(geometry=[box(500_000, 3_530_000, 510_000, 3_540_000)], crs=UTM)
    lon, lat = footprint_centroid_wgs84(fp)
    # UTM 12N central meridian is -111; centre easting 505 km is just east of it
    assert -111.0 < lon < -110.9
    assert 31.8 < lat < 32.0


# ---------------------------------------------------------------------------
# Classification — extent is the experimental variable
# ---------------------------------------------------------------------------


def _synthetic_scene(seed=0):
    rng = np.random.default_rng(seed)
    shape = (100, 100)
    silica = rng.normal(1.0, 0.05, shape)
    carbonate = rng.normal(1.0, 0.05, shape)
    mafic = rng.normal(1.0, 0.05, shape)
    # A coherent anomaly block elevated in all three ratios
    for arr in (silica, carbonate, mafic):
        arr[40:50, 40:50] += 0.5
    b10 = np.ones(shape)
    b10[:, :10] = np.nan  # no-data strip on the west edge
    for arr in (silica, carbonate, mafic):
        arr[np.isnan(b10)] = np.nan
    return silica, carbonate, mafic, b10


def test_classify_strong_zones_finds_anomaly_block():
    silica, carbonate, mafic, b10 = _synthetic_scene()
    zones, fp, rb = classify_strong_zones(
        silica, carbonate, mafic, b10, TRANSFORM, UTM, ClassificationParams()
    )
    assert len(zones) > 0
    block = box(500_000 + 40 * 90, 3_540_000 - 50 * 90, 500_000 + 50 * 90, 3_540_000 - 40 * 90)
    assert zones.geometry.union_all().contains(block.buffer(-1))
    # Footprint excludes the NaN strip: 90 of 100 columns valid
    assert fp.geometry.area.sum() == pytest.approx(90 * 100 * 90 * 90, rel=0.02)
    assert rb[0] < rb[2] and rb[1] < rb[3]


def test_percentiles_depend_on_extent():
    silica, carbonate, mafic, b10 = _synthetic_scene()
    cp = ClassificationParams()
    z_full, fp_full, _ = classify_strong_zones(silica, carbonate, mafic, b10, TRANSFORM, UTM, cp)
    # Sub-window dominated by the anomaly block: thresholds rise, so strong area shrinks there
    sub = (slice(35, 55), slice(35, 55))
    t_sub = TRANSFORM * rasterio.Affine.translation(35, 35)
    z_sub, fp_sub, _ = classify_strong_zones(
        silica[sub], carbonate[sub], mafic[sub], b10[sub], t_sub, UTM, cp
    )
    frac_full = z_full.geometry.area.sum() / fp_full.geometry.area.sum()
    frac_sub = z_sub.geometry.area.sum() / fp_sub.geometry.area.sum()
    assert frac_sub != pytest.approx(frac_full, abs=1e-3)


# ---------------------------------------------------------------------------
# Region restriction + hit test
# ---------------------------------------------------------------------------


def test_restrict_to_region_clips_zones_and_footprint():
    zones = gpd.GeoDataFrame(geometry=[box(0, 0, 10, 10), box(20, 20, 30, 30)], crs=UTM)
    fp = gpd.GeoDataFrame(geometry=[box(0, 0, 40, 40)], crs=UTM)
    region = box(0, 0, 40, 40).difference(box(0, 0, 15, 15))
    z, f = restrict_to_region(zones, fp, region)
    assert len(z) == 1
    assert z["area_km2"].sum() == pytest.approx(100 / 1e6)
    assert f.geometry.area.sum() == pytest.approx(40 * 40 - 15 * 15)


def test_summarize_hits_enrichment_and_p():
    r = summarize_hits(100, 30, 0.1)
    assert r["enrichment"] == pytest.approx(3.0)
    assert r["expected_hits_null"] == pytest.approx(10.0)
    assert r["p_binomial"] < 1e-6
    r0 = summarize_hits(0, 0, 0.1)
    assert r0["p_binomial"] == 1.0
    assert np.isnan(r0["enrichment"])


def _mrds_frame(points_utm):
    pts = gpd.GeoSeries([p for p, _ in points_utm], crs=UTM).to_crs("EPSG:4326")
    return pd.DataFrame({
        "longitude": pts.x.values,
        "latitude": pts.y.values,
        "commod1": [c for _, c in points_utm],
    })


def test_critical_hit_test_counts_only_critical_inside_footprint():
    from shapely.geometry import Point

    zones = gpd.GeoDataFrame(geometry=[box(500_000, 3_530_000, 502_000, 3_532_000)], crs=UTM)
    fp = gpd.GeoDataFrame(geometry=[box(500_000, 3_530_000, 510_000, 3_540_000)], crs=UTM)
    mrds = _mrds_frame([
        (Point(501_000, 3_531_000), "Copper"),       # critical, in zone
        (Point(501_500, 3_531_500), "Gold"),         # critical, in zone
        (Point(505_000, 3_535_000), "Zinc"),         # critical, outside zone
        (Point(501_200, 3_531_200), "Sand and Gravel, Construction"),  # non-critical, in zone
        (Point(520_000, 3_535_000), "Copper"),       # outside footprint
    ])
    rb = tuple(gpd.GeoSeries([fp.geometry.iloc[0]], crs=UTM).to_crs("EPSG:4326").total_bounds)
    rb = (rb[0] - 0.5, rb[1] - 0.5, rb[2] + 0.5, rb[3] + 0.5)
    r = critical_hit_test(mrds, zones, fp, rb)
    assert r["n_crit"] == 3
    assert r["hits_crit"] == 2
    assert r["p0_null"] == pytest.approx(4 / 100, abs=1e-4)
    assert r["enrichment"] == pytest.approx((2 / 3) / 0.04, rel=1e-2)


def test_critical_hit_test_empty_zones():
    from shapely.geometry import Point

    zones = gpd.GeoDataFrame(geometry=[], crs=UTM)
    fp = gpd.GeoDataFrame(geometry=[box(500_000, 3_530_000, 510_000, 3_540_000)], crs=UTM)
    mrds = _mrds_frame([(Point(505_000, 3_535_000), "Copper")])
    rb = (-112.0, 31.0, -110.0, 33.0)
    r = critical_hit_test(mrds, zones, fp, rb)
    assert r["n_crit"] == 1 and r["hits_crit"] == 0
    assert r["p_binomial"] == 1.0
