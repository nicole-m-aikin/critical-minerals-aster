"""HyperspectralScene / ThermalScene contract smoke tests."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_bounds

from critical_minerals_aster.prospectivity.adapters import (
    load_emit_hyperspectral_scene,
    percentile_rank,
)

ROOT = Path(__file__).resolve().parents[1]


def test_percentile_rank_is_within_unit_interval():
    arr = np.array([[1.0, 2.0], [3.0, np.nan]])
    mask = np.isfinite(arr)
    ranks = percentile_rank(arr, mask)
    assert ranks.shape == arr.shape
    assert np.isnan(ranks[1, 1])
    finite = ranks[mask]
    assert finite.min() >= 0.0 and finite.max() <= 1.0
    assert ranks[0, 0] < ranks[0, 1] < ranks[1, 0]


def test_emit_synthetic_adapter_contract(tmp_path):
    """Synthetic EMIT stack loads through the shared HyperspectralScene contract."""
    site_id = "demo"
    out_dir = tmp_path / "data" / "emit" / site_id
    out_dir.mkdir(parents=True)
    (tmp_path / "results").mkdir()
    west, south, east, north = -118.0, 41.0, -117.9, 41.1
    width, height = 16, 16
    transform = from_bounds(west, south, east, north, width, height)
    aloh = np.linspace(0, 0.5, width * height, dtype=np.float32).reshape(height, width)
    profile = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": transform,
        "nodata": np.nan,
    }
    for key, arr in (("aloh", aloh), ("kaolinite", aloh * 0.7), ("alunite", aloh * 0.4)):
        with rasterio.open(out_dir / f"{site_id}_{key}.tif", "w", **profile) as dst:
            dst.write(arr.astype(np.float32), 1)
    (tmp_path / "results" / f"emit_{site_id}_swir_provenance.json").write_text(
        json.dumps({
            "sensor": "EMIT_L2A",
            "features_used_nm": [2165.0, 2200.0, 2210.0],
            "mosaic_resolution_m": 60,
        })
    )

    scene = load_emit_hyperspectral_scene(site_id, repo_root=tmp_path)
    assert scene is not None
    assert scene.sensor_id.startswith("EMIT")
    assert scene.gsd_m == 60.0
    assert scene.aloh.shape == scene.kaolinite.shape == scene.alunite.shape
    assert scene.coverage_frac > 0
    assert len(scene.features_used_nm) == 3
