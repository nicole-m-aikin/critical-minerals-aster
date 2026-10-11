"""Tests for ag100.py — tile keys, HDF5 reading, mosaicking and mask reprojection (offline, synthetic)."""

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from rasterio.transform import from_origin
from rasterio.warp import transform as warp_transform

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

netCDF4 = pytest.importorskip("netCDF4")

from critical_minerals_aster.ag100 import (  # noqa: E402
    FILL_VALUE,
    NDVI_VEGETATION_THRESHOLD,
    ag100_masks,
    fetch_ag100_tiles,
    geolocation_transform,
    read_ag100_mosaic,
    read_ag100_tile,
    tile_key_from_name,
    tiles_for_bbox,
)

N = 11  # synthetic tile: 11 × 11 pixel centres over 1° → 0.1° spacing


def write_tile(
    repo_root: Path,
    lat_nw: int,
    lon_nw: int,
    *,
    ndvi_raw: int = 10,
    water_cells: tuple[tuple[int, int], ...] = (),
    fill_cells: tuple[tuple[int, int], ...] = (),
) -> Path:
    """AG100-layout HDF5 tile with pixel-centre geolocation NW→SE inclusive."""
    cache = repo_root / "data" / "ag100"
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"AG100.v003.{lat_nw}.{lon_nw}.0001.h5"
    lats = np.linspace(lat_nw, lat_nw - 1, N, dtype=np.float32)
    lons = np.linspace(lon_nw, lon_nw + 1, N, dtype=np.float32)
    lon2d, lat2d = np.meshgrid(lons, lats)

    emis = np.empty((5, N, N), dtype=np.int16)
    for i, v in enumerate((950, 960, 970, 980, 990)):
        emis[i] = v
    ndvi = np.full((N, N), ndvi_raw, dtype=np.int16)
    lw = np.zeros((N, N), dtype=np.int16)
    obs = np.full((N, N), 7, dtype=np.int16)
    for r, c in water_cells:
        lw[r, c] = 1
    for r, c in fill_cells:
        emis[:, r, c] = FILL_VALUE
        ndvi[r, c] = FILL_VALUE
        obs[r, c] = FILL_VALUE

    with netCDF4.Dataset(path, "w", format="NETCDF4") as ds:
        def var(group, name, data, dtype):
            g = ds.groups.get(group) or ds.createGroup(group)
            dims = []
            for k, n in enumerate(data.shape):
                d = f"{name}_d{k}"
                g.createDimension(d, n)
                dims.append(d)
            g.createVariable(name, dtype, tuple(dims))[:] = data

        var("Emissivity", "Mean", emis, "i2")
        var("NDVI", "Mean", ndvi, "i2")
        var("Land Water Map", "LWmap", lw, "i2")
        var("Observations", "NumObs", obs, "i2")
        var("Geolocation", "Latitude", lat2d, "f4")
        var("Geolocation", "Longitude", lon2d, "f4")
    return path


def test_tile_key_parsing():
    assert tile_key_from_name("AG100.v003.43.-118.0001.h5") == (43, -118)
    assert tile_key_from_name("AG100.v003.42.-113.0001") == (42, -113)
    assert tile_key_from_name("AG100.v003.40.-077.0001.h5") == (40, -77)
    assert tile_key_from_name("AG1km.v003.43.-118.0010.h5") is None


def test_tiles_for_bbox_uses_nw_corner():
    # McDermitt bbox spans 41.8–42.4 N, 118.1–117.3 W → four 1° tiles
    assert sorted(tiles_for_bbox((-118.1, 41.8, -117.3, 42.4))) == [
        (42, -119), (42, -118), (43, -119), (43, -118)
    ]
    # Integer edges do not pull in the neighbouring tile
    assert tiles_for_bbox((-118.0, 42.0, -117.0, 43.0)) == [(43, -118)]


def test_geolocation_transform_pixel_centres():
    lats = np.linspace(43, 42, N)
    lons = np.linspace(-118, -117, N)
    lon2d, lat2d = np.meshgrid(lons, lats)
    t = geolocation_transform(lat2d, lon2d)
    assert t.a == pytest.approx(0.1) and t.e == pytest.approx(-0.1)
    assert t.c == pytest.approx(-118.05) and t.f == pytest.approx(43.05)
    with pytest.raises(ValueError):
        geolocation_transform(lat2d + np.random.default_rng(0).normal(0, 0.05, lat2d.shape), lon2d)


def test_read_tile_scales_and_fill(tmp_path):
    p = write_tile(tmp_path, 43, -118, ndvi_raw=45, water_cells=((2, 3),), fill_cells=((5, 5),))
    arrs, t = read_ag100_tile(p, ("emissivity", "ndvi", "water", "obs_count"))
    assert arrs["emissivity"].shape == (5, N, N)
    assert arrs["emissivity"][:, 0, 0] == pytest.approx([0.95, 0.96, 0.97, 0.98, 0.99])
    assert arrs["ndvi"][0, 0] == pytest.approx(0.45)
    assert arrs["water"][2, 3] == 1.0 and arrs["water"][0, 0] == 0.0
    assert arrs["obs_count"][0, 0] == 7
    assert np.isnan(arrs["emissivity"][:, 5, 5]).all()
    assert np.isnan(arrs["ndvi"][5, 5]) and np.isnan(arrs["obs_count"][5, 5])
    assert t.c == pytest.approx(-118.05)


def test_masks_on_native_grid_exact(tmp_path):
    write_tile(tmp_path, 43, -118, ndvi_raw=45, water_cells=((2, 3),), fill_cells=((5, 5),))
    dst_t = from_origin(-118.05, 43.05, 0.1, 0.1)
    masks = ag100_masks((-118, 42, -117, 43), dst_t, (N, N), "EPSG:4326", repo_root=tmp_path, login=False)
    assert set(masks) == {"water", "ndvi", "obs_count"}
    assert masks["water"].dtype == bool
    assert masks["water"][2, 3] and masks["water"].sum() == 1
    assert masks["ndvi"][0, 0] == pytest.approx(0.45)
    assert masks["obs_count"][0, 0] == 7
    assert np.isnan(masks["obs_count"][5, 5])
    assert (masks["ndvi"][np.isfinite(masks["ndvi"])] > NDVI_VEGETATION_THRESHOLD).all()


def test_mosaic_two_tiles_shares_seam(tmp_path):
    write_tile(tmp_path, 43, -118, ndvi_raw=10)
    write_tile(tmp_path, 43, -117, ndvi_raw=50)
    m = read_ag100_mosaic((-118, 42, -116, 43), tmp_path, ("ndvi",), login=False, pad_pixels=0)
    ndvi = m.layers["ndvi"]
    assert len(m.tiles) == 2
    # 11 + 11 columns sharing the -117° centre column → 21 columns
    assert ndvi.shape == (N, 2 * N - 1)
    assert not np.isnan(ndvi).any()
    assert ndvi[0, 0] == pytest.approx(0.10) and ndvi[0, -1] == pytest.approx(0.50)
    assert m.transform.c == pytest.approx(-118.05)


def test_masks_reproject_to_utm(tmp_path):
    write_tile(tmp_path, 43, -118, water_cells=((5, 5),))
    # 100 m UTM 11N grid centred on the water pixel (42.5 N, -117.5 W)
    (x,), (y,) = warp_transform("EPSG:4326", "EPSG:32611", [-117.5], [42.5])
    dst_t = from_origin(x - 2000, y + 2000, 100, 100)
    masks = ag100_masks((-118, 42, -117, 43), dst_t, (40, 40), "EPSG:32611", repo_root=tmp_path, login=False)
    assert masks["water"][20, 20]  # centre of the grid lies in the 0.1° water pixel
    assert masks["water"].mean() > 0.5
    assert masks["obs_count"][20, 20] == 7
    assert np.isfinite(masks["ndvi"]).all()


def test_missing_tiles_return_none_with_warning(tmp_path):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        assert fetch_ag100_tiles((-118, 42, -117, 43), tmp_path, login=False) == []
        out = ag100_masks(
            (-118, 42, -117, 43), from_origin(-118, 43, 0.1, 0.1), (10, 10), "EPSG:4326",
            repo_root=tmp_path, login=False,
        )
    assert out is None
    assert any("AG100" in str(x.message) for x in w)


def test_absent_tile_registry_skips_search(tmp_path):
    cache = tmp_path / "data" / "ag100"
    cache.mkdir(parents=True)
    (cache / "absent_tiles.json").write_text("[[43, -118]]")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        # Absent tiles are neither searched nor warned about as missing
        assert fetch_ag100_tiles((-118, 42, -117, 43), tmp_path, login=False) == []


def test_vegetation_threshold_constant():
    assert NDVI_VEGETATION_THRESHOLD == 0.3
