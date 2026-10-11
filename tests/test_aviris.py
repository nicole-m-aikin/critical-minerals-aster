import sys
from pathlib import Path

import numpy as np
import pytest
from rasterio.transform import from_origin

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.aviris import (  # noqa: E402
    crop_to_valid,
    mosaic_to_bbox,
    select_granules_painting_bbox,
    valid_window,
)


def test_valid_window_trims_nan_margins():
    arr = np.full((20, 20), np.nan)
    arr[5:12, 3:10] = 0.2
    r0, r1, c0, c1 = valid_window(arr, pad=1)
    assert (r0, r1, c0, c1) == (4, 13, 2, 11)


def test_crop_to_valid_updates_affine():
    arr = np.full((10, 10), np.nan)
    arr[2:6, 4:8] = 1.0
    transform = from_origin(500_000.0, 4_200_000.0, 5.0, 5.0)
    cropped, new_t = crop_to_valid(arr, transform, pad=0)
    assert cropped.shape == (4, 4)
    # Window starts at col 4, row 2 → origin shifts 20 m east, 10 m south.
    assert new_t.c == 500_000.0 + 4 * 5.0
    assert new_t.f == 4_200_000.0 - 2 * 5.0


def test_mosaic_to_bbox_keeps_first_valid_and_matches_window():
    """A 2×2 patch at known UTM coords must land inside the WGS84 bbox grid."""
    # Pixel at UTM 11N ~ Goldfield town (479200 E, 4173500 N).
    arr = np.array([[0.4]], dtype=np.float32)
    transform = from_origin(479_195.0, 4_173_505.0, 10.0, 10.0)
    bbox = (-117.24, 37.705, -117.23, 37.712)
    mosaic, new_t, crs = mosaic_to_bbox(
        [arr], [transform], "EPSG:32611", bbox, resolution=10.0,
    )
    assert "32611" in crs or "UTM zone 11" in crs
    assert mosaic.shape[0] > 1 and mosaic.shape[1] > 1
    assert np.isfinite(mosaic).any()
    assert float(np.nanmax(mosaic)) == pytest.approx(0.4, abs=1e-6)


def test_select_granules_painting_bbox_greedy_and_seeds(tmp_path):
    """A redundant overlapping strip must lose to two complementary strips."""
    import json

    cache = tmp_path / "cache.json"
    cache.write_text(json.dumps([
        {"granule_ur": "west",
         "geometry_wkt": "POLYGON ((-2 0, -1 0, -1 1, -2 1, -2 0))"},
        {"granule_ur": "east",
         "geometry_wkt": "POLYGON ((1 0, 2 0, 2 1, 1 1, 1 0))"},
        {"granule_ur": "west_dup",
         "geometry_wkt": "POLYGON ((-1.9 0.1, -1.1 0.1, -1.1 0.9, -1.9 0.9, -1.9 0.1))"},
        {"granule_ur": "seed_only",
         "geometry_wkt": "POLYGON ((-0.2 0.4, 0.2 0.4, 0.2 0.6, -0.2 0.6, -0.2 0.4))"},
    ]))
    bbox = (-2.0, 0.0, 2.0, 1.0)
    picked = select_granules_painting_bbox(
        cache, bbox, must_contain=((0.0, 0.5),), min_gain=0.10,
    )
    assert picked[0] == "seed_only"
    assert "west" in picked
    assert "east" in picked
    assert "west_dup" not in picked
