"""ASTER Global Emissivity Dataset AG100 v003: water / vegetation / observation masks.

AG100 (Hulley et al. 2015, doi:10.1002/2015GL065564; LP DAAC short name
``AG100``, version ``003``) is the 2000–2008 mean of every clear-sky ASTER
TES retrieval on 1°×1° geographic tiles at ~100 m (1000 × 1000 pixels).
Each HDF5 tile is named ``AG100.v003.{NW lat}.{NW lon}.0001.h5`` and holds
(verified on AG100.v003.43.-118.0001.h5 and AG100.v003.42.-113.0001.h5):

=============================  =====  ===========  =======================
dataset                        dtype  scale        notes
=============================  =====  ===========  =======================
``Emissivity/Mean`` (5, y, x)  int16  0.001        bands 10–14
``Emissivity/SDev``            int16  0.0001
``NDVI/Mean``                  int16  0.01         from TOA reflectance
``Land Water Map/LWmap``       int16  n/a          1 = water, 0 = land
``Observations/NumObs``        int16  n/a          scenes in the mean
``Temperature/Mean``           int32  0.01 K
``Geolocation/Latitude``       f32    n/a          pixel centres
``Geolocation/Longitude``      f32    n/a          pixel centres
=============================  =====  ===========  =======================

No scale / fill attributes are stored in the file; scales and the −9999
fill are from the ASTER GED v3 user guide (LP DAAC document 120). The guide
does not state the LWmap coding: 1 = water was checked against the Great
Salt Lake (south arm = 1, flat ~0.985 emissivity, NDVI 0.03; shore / Ogden
= 0). Geolocation runs from the NW integer corner to the SE corner
inclusive (1000 centres over 1°, spacing 1/999°), so neighbouring tiles
share their edge row / column.

HDF5 is read through ``netCDF4`` (the GDAL bundled with the rasterio wheels
has no HDF5 driver); plain HDF5 files open with phony dimensions.
"""

from __future__ import annotations

import json
import math
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from rasterio.crs import CRS
from rasterio.transform import Affine
from rasterio.warp import Resampling, reproject

AG100_SHORT_NAME = "AG100"
AG100_VERSION = "003"
AG100_CRS = CRS.from_epsg(4326)

FILL_VALUE = -9999
EMISSIVITY_SCALE = 0.001
NDVI_SCALE = 0.01
LWMAP_WATER = 1
EMISSIVITY_BANDS: tuple[int, ...] = (10, 11, 12, 13, 14)

NDVI_VEGETATION_THRESHOLD = 0.3
"""Mean-NDVI cutoff above which a pixel is treated as vegetated and masked.

TIR emissivity indices (quartz, carbonate, mafic) are mineral-spectrum
contrasts; green vegetation is a near-blackbody (ε ≈ 0.97–0.99, spectrally
flat in 8–12 µm) that dilutes or erases them, so index values over
vegetated pixels are unreliable. The NDVI-threshold emissivity method
(Sobrino et al. 2008) treats NDVI < 0.2 as bare soil and > 0.5 as full
cover; 0.3 sits at the low end of the mixed range and is the common
cutoff for masking vegetation before geological band-ratio mapping.
AG100 NDVI is from TOA reflectance, which typically reads slightly lower
than surface NDVI, so the effective cutoff is a little permissive.
"""

MIN_OBS_COUNT = 3
"""Pixels averaged from fewer than this many ASTER scenes are treated as unreliable."""

_DATASETS: dict[str, str] = {
    "emissivity": "Emissivity/Mean",
    "ndvi": "NDVI/Mean",
    "water": "Land Water Map/LWmap",
    "obs_count": "Observations/NumObs",
}
_BILINEAR = frozenset({"emissivity", "ndvi"})
_TILE_RE = re.compile(r"^AG100\.v003\.(-?\d+)\.(-?\d+)\.0001(?:\.h5)?$")
_ABSENT_FILE = "absent_tiles.json"

TileKey = tuple[int, int]


def ag100_cache_dir(repo_root: Path) -> Path:
    return Path(repo_root) / "data" / "ag100"


def tile_key_from_name(name: str) -> TileKey | None:
    """(NW lat, NW lon) from an AG100 granule id or file name; None if not AG100 v003."""
    m = _TILE_RE.match(Path(name).name)
    return (int(m.group(1)), int(m.group(2))) if m else None


def tiles_for_bbox(bbox_wgs84: Sequence[float]) -> list[TileKey]:
    """NW-corner keys of every 1° tile intersecting ``(west, south, east, north)``."""
    west, south, east, north = bbox_wgs84
    lats = range(math.floor(south) + 1, math.ceil(north) + 1)
    lons = range(math.floor(west), math.ceil(east))
    return [(lat, lon) for lat in lats for lon in lons]


def _cached_tiles(cache_dir: Path) -> dict[TileKey, Path]:
    out: dict[TileKey, Path] = {}
    if cache_dir.is_dir():
        for p in cache_dir.glob("AG100.v003.*.h5"):
            key = tile_key_from_name(p.name)
            if key is not None:
                out[key] = p
    return out


def _load_absent(cache_dir: Path) -> set[TileKey]:
    path = cache_dir / _ABSENT_FILE
    if not path.is_file():
        return set()
    return {tuple(k) for k in json.loads(path.read_text())}


def _save_absent(cache_dir: Path, absent: set[TileKey]) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / _ABSENT_FILE).write_text(json.dumps(sorted(list(k) for k in absent)))


def fetch_ag100_tiles(
    bbox_wgs84: Sequence[float],
    repo_root: Path,
    login: bool = True,
) -> list[Path]:
    """Cache-first AG100 tiles covering ``bbox_wgs84`` (downloads missing ones).

    Tiles live in ``data/ag100/``. Tiles CMR has no granule for (open ocean)
    are recorded in ``data/ag100/absent_tiles.json`` so they are not searched
    again. With ``login=False`` only cached tiles are returned. The list may
    be shorter than the tile set if some tiles are absent or unobtainable.
    """
    cache_dir = ag100_cache_dir(repo_root)
    needed = tiles_for_bbox(bbox_wgs84)
    have = _cached_tiles(cache_dir)
    absent = _load_absent(cache_dir)
    missing = [k for k in needed if k not in have and k not in absent]

    if missing and login:
        try:
            import earthaccess

            earthaccess.login(strategy="netrc")
            results = earthaccess.search_data(
                short_name=AG100_SHORT_NAME, version=AG100_VERSION, bounding_box=tuple(bbox_wgs84)
            )
            by_key = {tile_key_from_name(g["umm"]["GranuleUR"]): g for g in results}
            wanted = [by_key[k] for k in missing if k in by_key]
            if wanted:
                cache_dir.mkdir(parents=True, exist_ok=True)
                earthaccess.download(wanted, str(cache_dir), threads=2, show_progress=False)
            newly_absent = {k for k in missing if k not in by_key}
            if newly_absent:
                _save_absent(cache_dir, absent | newly_absent)
        except Exception as exc:  # network / auth / CMR failure — fall back to cache
            warnings.warn(f"AG100 download failed for bbox {tuple(bbox_wgs84)}: {exc!r}", RuntimeWarning, stacklevel=2)
        have = _cached_tiles(cache_dir)
    elif missing:
        warnings.warn(
            f"AG100 tiles {missing} not cached and login=False; using cached tiles only",
            RuntimeWarning,
            stacklevel=2,
        )
    return [have[k] for k in needed if k in have]


def geolocation_transform(lat: np.ndarray, lon: np.ndarray) -> Affine:
    """Affine (pixel-edge origin) from 2-D pixel-centre latitude / longitude grids.

    Raises ``ValueError`` if the grid is not a regular lat/lon lattice.
    """
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    ny, nx = lat.shape
    lat_col, lon_row = lat[:, 0], lon[0, :]
    dy = (lat_col[-1] - lat_col[0]) / (ny - 1)
    dx = (lon_row[-1] - lon_row[0]) / (nx - 1)
    tol_y, tol_x = 0.25 * abs(dy), 0.25 * abs(dx)
    if (
        np.ptp(lat, axis=1).max() > tol_y
        or np.ptp(lon, axis=0).max() > tol_x
        or np.abs(lat_col - (lat_col[0] + dy * np.arange(ny))).max() > tol_y
        or np.abs(lon_row - (lon_row[0] + dx * np.arange(nx))).max() > tol_x
    ):
        raise ValueError("AG100 geolocation is not a regular lat/lon grid")
    return Affine(dx, 0.0, lon_row[0] - dx / 2, 0.0, dy, lat_col[0] - dy / 2)


def _scale(name: str, raw: np.ndarray) -> np.ndarray:
    fill = raw == FILL_VALUE
    if name == "emissivity":
        arr = raw.astype(np.float64) * EMISSIVITY_SCALE
        arr[fill | (raw <= 0)] = np.nan
    elif name == "ndvi":
        arr = raw.astype(np.float64) * NDVI_SCALE
        arr[fill] = np.nan
    elif name == "water":
        arr = (raw == LWMAP_WATER).astype(np.float64)
        arr[fill] = np.nan
    elif name == "obs_count":
        arr = raw.astype(np.float64)
        arr[fill | (raw < 0)] = np.nan
    else:
        raise KeyError(name)
    return arr


def read_ag100_tile(path: Path, layers: Iterable[str]) -> tuple[dict[str, np.ndarray], Affine]:
    """Scaled float layers (NaN = fill) and the tile's EPSG:4326 transform.

    ``emissivity`` is (5, y, x) for bands 10–14; ``water`` is 1.0 water /
    0.0 land; ``ndvi`` is unitless; ``obs_count`` is a scene count.
    """
    import netCDF4

    with netCDF4.Dataset(str(path)) as ds:
        ds.set_auto_maskandscale(False)
        transform = geolocation_transform(ds["Geolocation/Latitude"][:], ds["Geolocation/Longitude"][:])
        out = {name: _scale(name, np.asarray(ds[_DATASETS[name]][:])) for name in layers}
    return out, transform


@dataclass
class AG100Mosaic:
    """Native-resolution AG100 layers mosaicked over a bbox (EPSG:4326)."""

    layers: dict[str, np.ndarray]
    transform: Affine
    crs: CRS = field(default_factory=lambda: AG100_CRS)
    tiles: list[Path] = field(default_factory=list)

    @property
    def shape(self) -> tuple[int, int]:
        return next(iter(self.layers.values())).shape[-2:]


def _paste(dst: np.ndarray, src: np.ndarray, row_off: int, col_off: int) -> None:
    """Copy ``src`` into ``dst`` at offset, filling only NaN destination pixels."""
    h, w = dst.shape[-2:]
    sh, sw = src.shape[-2:]
    r0, c0 = max(row_off, 0), max(col_off, 0)
    r1, c1 = min(row_off + sh, h), min(col_off + sw, w)
    if r0 >= r1 or c0 >= c1:
        return
    d = dst[..., r0:r1, c0:c1]
    s = src[..., r0 - row_off:r1 - row_off, c0 - col_off:c1 - col_off]
    take = np.isnan(d) & ~np.isnan(s)
    d[take] = s[take]


def read_ag100_mosaic(
    bbox_wgs84: Sequence[float],
    repo_root: Path,
    layers: Sequence[str] = ("water", "ndvi", "obs_count"),
    *,
    login: bool = True,
    pad_pixels: int = 2,
) -> AG100Mosaic | None:
    """Mosaic AG100 tiles on their native lat/lon lattice, cropped to the bbox.

    ``pad_pixels`` extra pixels on every side keep bilinear resampling at
    the bbox edge supported. Returns None (with a warning) if no tile could
    be read.
    """
    paths = fetch_ag100_tiles(bbox_wgs84, repo_root, login=login)
    tiles: list[tuple[dict[str, np.ndarray], Affine, Path]] = []
    for p in paths:
        try:
            arrs, t = read_ag100_tile(p, layers)
        except Exception as exc:  # truncated download / unexpected layout
            warnings.warn(f"AG100 tile {p.name} unreadable: {exc!r}", RuntimeWarning, stacklevel=2)
            continue
        tiles.append((arrs, t, p))
    if not tiles:
        warnings.warn(f"No AG100 tiles available for bbox {tuple(bbox_wgs84)}", RuntimeWarning, stacklevel=2)
        return None

    t0 = tiles[0][1]
    dx, dy = t0.a, t0.e
    west, south, east, north = bbox_wgs84
    col0 = math.floor((west - t0.c) / dx) - pad_pixels
    col1 = math.ceil((east - t0.c) / dx) + pad_pixels
    row0 = math.floor((north - t0.f) / dy) - pad_pixels
    row1 = math.ceil((south - t0.f) / dy) + pad_pixels
    transform = Affine(dx, 0.0, t0.c + col0 * dx, 0.0, dy, t0.f + row0 * dy)
    h, w = row1 - row0, col1 - col0

    mosaic: dict[str, np.ndarray] = {}
    for name in layers:
        nb = tiles[0][0][name].shape[:-2]
        mosaic[name] = np.full((*nb, h, w), np.nan)
    for arrs, t, p in tiles:
        roff = (t.f - transform.f) / dy
        coff = (t.c - transform.c) / dx
        if abs(roff - round(roff)) > 0.01 or abs(coff - round(coff)) > 0.01 or not (
            math.isclose(t.a, dx, rel_tol=1e-4) and math.isclose(t.e, dy, rel_tol=1e-4)
        ):
            raise ValueError(f"AG100 tile {p.name} is not on the mosaic lattice")
        for name in layers:
            _paste(mosaic[name], arrs[name], int(round(roff)), int(round(coff)))
    return AG100Mosaic(layers=mosaic, transform=transform, tiles=[p for _, _, p in tiles])


def mosaic_to_grid(
    mosaic: AG100Mosaic,
    dst_transform: Affine,
    dst_shape: tuple[int, int],
    dst_crs,
) -> dict[str, np.ndarray]:
    """Reproject mosaic layers to a destination grid (NaN where uncovered).

    Nearest for categorical / count layers (``water``, ``obs_count``),
    bilinear for continuous ones (``ndvi``, ``emissivity``).
    """
    out: dict[str, np.ndarray] = {}
    for name, src in mosaic.layers.items():
        method = Resampling.bilinear if name in _BILINEAR else Resampling.nearest
        dst = np.full((*src.shape[:-2], *dst_shape), np.nan)
        reproject(
            source=src,
            destination=dst,
            src_transform=mosaic.transform,
            src_crs=mosaic.crs,
            src_nodata=np.nan,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            dst_nodata=np.nan,
            resampling=method,
        )
        out[name] = dst
    return out


def ag100_masks(
    bbox_wgs84: Sequence[float],
    dst_transform: Affine,
    dst_shape: tuple[int, int],
    dst_crs,
    *,
    repo_root: Path,
    login: bool = True,
) -> dict[str, np.ndarray] | None:
    """AG100 water / NDVI / observation-count layers on a destination grid.

    Returns ``{"water": bool, "ndvi": float, "obs_count": float}`` arrays of
    ``dst_shape`` (NaN where AG100 has no data; ``water`` False there), or
    None with a warning if no tile can be obtained. Only tiles intersecting
    ``bbox_wgs84`` are used, so the bbox should cover the destination grid.
    """
    mosaic = read_ag100_mosaic(bbox_wgs84, repo_root, ("water", "ndvi", "obs_count"), login=login)
    if mosaic is None:
        return None
    return masks_from_layers(mosaic_to_grid(mosaic, dst_transform, dst_shape, dst_crs))


def masks_from_layers(grids: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """``ag100_masks`` output from already-resampled ``water`` / ``ndvi`` / ``obs_count`` grids."""
    return {
        "water": np.nan_to_num(grids["water"], nan=0.0) == 1.0,
        "ndvi": grids["ndvi"],
        "obs_count": grids["obs_count"],
    }
