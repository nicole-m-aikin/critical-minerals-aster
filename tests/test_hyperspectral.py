import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Polygon, box

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.hyperspectral import (  # noqa: E402
    COLLECTIONS,
    GEMX_EPOCH,
    SWIR_COLLECTION_KEYS,
    _granule_polygon,
    _in_gemx_era,
    aloh_band_depth,
    continuum_removed_band_depth,
    coverage_fraction,
    doublet_score,
    nearest_wavelength_index,
    summarize_by_site,
)


def test_collections_have_unique_keys_and_concept_ids():
    keys = [c.key for c in COLLECTIONS]
    ids = [c.concept_id for c in COLLECTIONS]
    assert len(keys) == len(set(keys))
    assert len(ids) == len(set(ids))


def test_only_aviris_is_swir_capable():
    """MASTER spans SWIR but with ~50 bands across 0.46-12.9 um, too coarse
    to resolve the 2.2 um Al-OH doublet, so it must not count as full SWIR."""
    assert set(SWIR_COLLECTION_KEYS) == {
        "aviris5_l2a_rfl",
        "avirisc_l2_rfl",
    }
    assert not any(c.sensor == "AVIRIS-3" for c in COLLECTIONS)
    master = [c for c in COLLECTIONS if c.sensor == "MASTER"]
    assert master and all(c.swir == "partial" for c in master)


def test_granule_polygon_from_gpolygon():
    umm = {
        "SpatialExtent": {
            "HorizontalSpatialDomain": {
                "Geometry": {
                    "GPolygons": [
                        {
                            "Boundary": {
                                "Points": [
                                    {"Longitude": -110.0, "Latitude": 31.0},
                                    {"Longitude": -109.0, "Latitude": 31.0},
                                    {"Longitude": -109.0, "Latitude": 32.0},
                                    {"Longitude": -110.0, "Latitude": 32.0},
                                    {"Longitude": -110.0, "Latitude": 31.0},
                                ]
                            }
                        }
                    ]
                }
            }
        }
    }
    poly = _granule_polygon(umm)
    assert isinstance(poly, Polygon)
    assert poly.bounds == (-110.0, 31.0, -109.0, 32.0)


def test_granule_polygon_from_bounding_rectangle():
    umm = {
        "SpatialExtent": {
            "HorizontalSpatialDomain": {
                "Geometry": {
                    "BoundingRectangles": [
                        {
                            "WestBoundingCoordinate": -110.0,
                            "SouthBoundingCoordinate": 31.0,
                            "EastBoundingCoordinate": -109.0,
                            "NorthBoundingCoordinate": 32.0,
                        }
                    ]
                }
            }
        }
    }
    assert _granule_polygon(umm).bounds == (-110.0, 31.0, -109.0, 32.0)


def test_granule_polygon_missing_geometry():
    assert _granule_polygon({}) is None
    assert _granule_polygon({"SpatialExtent": {}}) is None


def test_coverage_fraction_full_and_empty():
    bbox = (-110.15, 31.25, -109.65, 31.7)
    assert coverage_fraction([box(*bbox)], bbox) == pytest.approx(1.0, abs=1e-6)
    assert coverage_fraction([], bbox) == 0.0
    assert coverage_fraction([None], bbox) == 0.0


def test_coverage_fraction_partial_overlap():
    """A flight line clipping half the bbox must not read as full coverage."""
    bbox = (-110.0, 31.0, -109.0, 32.0)
    half = box(-110.0, 31.0, -109.5, 32.0)
    assert coverage_fraction([half], bbox) == pytest.approx(0.5, abs=0.01)


def test_coverage_fraction_is_capped_by_bbox():
    bbox = (-110.0, 31.0, -109.0, 32.0)
    huge = box(-120.0, 25.0, -100.0, 40.0)
    assert coverage_fraction([huge], bbox) == pytest.approx(1.0, abs=1e-6)


def test_coverage_fraction_unions_overlapping_lines():
    """Overlapping granules must not double-count toward the fraction."""
    bbox = (-110.0, 31.0, -109.0, 32.0)
    a = box(-110.0, 31.0, -109.5, 32.0)
    b = box(-109.75, 31.0, -109.25, 32.0)
    frac = coverage_fraction([a, b], bbox)
    assert frac == pytest.approx(0.75, abs=0.01)


def test_in_gemx_era_filters_legacy_flights():
    """AVIRIS collections are facility-wide archives going back to 2006, so
    pre-campaign flights must be excluded from GEMx-era counts."""
    assert _in_gemx_era("2025-06-01T00:00:00.000Z")
    assert _in_gemx_era(f"{GEMX_EPOCH}T00:00:00.000Z")
    assert not _in_gemx_era("2014-06-01T00:00:00.000Z")
    assert not _in_gemx_era(None)


def test_summarize_by_site_picks_best_sensor_and_tiers():
    coverage = pd.DataFrame([
        {"site_id": "a", "collection": "aviris5_l2a_rfl", "sensor": "AVIRIS-5",
         "coverage_frac": 0.9, "n_granules_gemx_era": 10},
        {"site_id": "a", "collection": "avirisc_l2_rfl", "sensor": "AVIRIS-Classic",
         "coverage_frac": 0.2, "n_granules_gemx_era": 3},
        {"site_id": "a", "collection": "master_gemx_2025", "sensor": "MASTER",
         "coverage_frac": 1.0, "n_granules_gemx_era": 40},
        {"site_id": "b", "collection": "aviris5_l2a_rfl", "sensor": "AVIRIS-5",
         "coverage_frac": 0.0, "n_granules_gemx_era": 0},
        {"site_id": "b", "collection": "master_gemx_2025", "sensor": "MASTER",
         "coverage_frac": 0.0, "n_granules_gemx_era": 0},
    ])
    out = summarize_by_site(coverage).set_index("site_id")

    assert out.loc["a", "swir_coverage_frac"] == pytest.approx(0.9)
    assert out.loc["a", "swir_best_sensor"] == "AVIRIS-5"
    assert out.loc["a", "swir_granules"] == 13
    # MASTER is tracked separately and must not inflate the SWIR figure.
    assert out.loc["a", "master_coverage_frac"] == pytest.approx(1.0)
    assert bool(out.loc["a", "has_swir"])
    assert out.loc["a", "swir_tier"] == "good"

    assert not bool(out.loc["b", "has_swir"])
    assert out.loc["b", "swir_tier"] == "none"


def test_nearest_wavelength_index():
    w = [380.0, 385.0, 2200.0, 2205.0, 2500.0]
    assert nearest_wavelength_index(w, 2200.0) == 2
    assert nearest_wavelength_index(w, 2202.0) == 2


def test_band_depth_is_zero_on_a_flat_spectrum():
    depth = continuum_removed_band_depth(0.4, 0.4, 0.4, 2120.0, 2200.0, 2270.0)
    assert depth == pytest.approx(0.0)


def test_band_depth_positive_for_an_absorption():
    """A 50% dip at the feature against a flat 0.4 continuum is depth 0.5."""
    depth = continuum_removed_band_depth(0.4, 0.2, 0.4, 2120.0, 2200.0, 2270.0)
    assert depth == pytest.approx(0.5)


def test_band_depth_nan_on_nodata_and_zero():
    assert np.isnan(continuum_removed_band_depth(0.0, 0.2, 0.4, 2120, 2200, 2270))
    assert np.isnan(continuum_removed_band_depth(-9999, 0.2, 0.4, 2120, 2200, 2270))


def test_band_depth_rejects_feature_outside_shoulders():
    with pytest.raises(ValueError):
        continuum_removed_band_depth(0.4, 0.2, 0.4, 2200.0, 2120.0, 2270.0)


def test_aloh_band_depth_finds_the_2200_nm_channel():
    w = np.arange(380.0, 2505.0, 5.0)
    cube = np.full((w.size, 3, 3), 0.4)
    i = nearest_wavelength_index(w, 2200.0)
    cube[i] = 0.2
    depth = aloh_band_depth(cube, w)
    assert depth.shape == (3, 3)
    assert depth[0, 0] == pytest.approx(0.5, abs=0.02)


def test_doublet_score_requires_both_features():
    """A deep dip at only one wavelength is not a kaolinite/alunite doublet."""
    both = doublet_score(0.20, 0.10)
    assert both == pytest.approx(0.10)
    only_one = doublet_score(0.20, 0.00)
    assert np.isnan(only_one)
    missing = doublet_score(0.20, np.nan)
    assert np.isnan(missing)
