"""Satellite imagery basemap fetch and cache for figure 03 context panel."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    import rasterio
    from critical_minerals_aster.config import SiteConfig
    from critical_minerals_aster.paths import SitePaths

FIG03_BASEMAP_CACHE_VERSION = 2


def _hillshade_bounds_wgs84(
    hs_transform: rasterio.Affine,
    hs_shape: tuple[int, int],
    crs: rasterio.crs.CRS,
) -> tuple[float, float, float, float]:
    from rasterio.warp import transform_bounds

    rows, cols = hs_shape
    t = hs_transform
    left = t.c
    right = t.c + t.a * cols
    bottom = t.f + t.e * rows
    top = t.f
    return transform_bounds(crs, "EPSG:4326", left, bottom, right, top)


def _fetch_tiles(
    west_3857: float,
    south_3857: float,
    east_3857: float,
    north_3857: float,
    zoom: int,
) -> tuple[np.ndarray, str, tuple[float, float, float, float]] | None:
    """Download RGB tiles in EPSG:3857.

    Accepts and returns bounds in EPSG:3857.  contextily always tiles in Web
    Mercator; using EPSG:3857 bounds directly (ll=False) avoids the ambiguity
    of how the returned extent is expressed when ll=True.

    Returns ``(rgb_array, provider_name, extent_3857)`` where
    ``extent_3857 = (left, right, bottom, top)`` in EPSG:3857.
    """
    import contextily as cx

    providers: list[tuple[str, Any]] = [
        ("Esri.WorldImagery", cx.providers.Esri.WorldImagery),
        ("CartoDB.Positron", cx.providers.CartoDB.Positron),
    ]
    for name, source in providers:
        try:
            img, extent = cx.bounds2img(
                west_3857, south_3857, east_3857, north_3857,
                zoom=zoom,
                source=source,
                ll=False,  # input is already EPSG:3857
            )
            if img is None or img.size == 0:
                continue
            arr = np.asarray(img)
            if arr.ndim == 2:
                arr = np.stack([arr, arr, arr], axis=-1)
            elif arr.shape[-1] == 4:
                arr = arr[..., :3]
            # extent from contextily is (left, right, bottom, top) in EPSG:3857
            return arr.astype(np.uint8), name, extent
        except Exception as exc:
            print(f"  [basemap] {name} failed: {exc}", file=sys.stderr)
    return None


def _reproject_rgb_to_grid(
    rgb_src: np.ndarray,
    src_bounds_lrbt: tuple[float, float, float, float],
    src_crs: str,
    dst_transform: rasterio.Affine,
    dst_crs: rasterio.crs.CRS,
    dst_shape: tuple[int, int],
) -> np.ndarray:
    """Reproject an (H, W, 3) RGB array onto the destination raster grid.

    Parameters
    ----------
    src_bounds_lrbt:
        Source image bounds as ``(left, right, bottom, top)`` in ``src_crs``.
    src_crs:
        CRS of the source image (e.g. ``"EPSG:3857"`` for contextily tiles).
    """
    from rasterio.transform import from_bounds
    from rasterio.warp import Resampling, reproject

    rows_s, cols_s = rgb_src.shape[0], rgb_src.shape[1]
    left, right, bottom, top = src_bounds_lrbt
    src_transform = from_bounds(left, bottom, right, top, cols_s, rows_s)
    dst_rows, dst_cols = dst_shape
    out = np.zeros((3, dst_rows, dst_cols), dtype=np.uint8)

    for band in range(3):
        reproject(
            source=rgb_src[:, :, band],
            destination=out[band],
            src_transform=src_transform,
            src_crs=src_crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.bilinear,
        )
    return np.transpose(out, (1, 2, 0))


def fetch_satellite_basemap_for_site(
    site: SiteConfig,
    paths: SitePaths,
    hs_transform: rasterio.Affine,
    hs_shape: tuple[int, int],
    crs: rasterio.crs.CRS,
    *,
    zoom: int = 12,
) -> tuple[np.ndarray, str, bool] | None:
    """Fetch or load cached satellite RGB aligned to the hillshade grid.

    Tiles are fetched in EPSG:3857 (Web Mercator) — the native CRS of all web
    map tile providers — and reprojected to the hillshade UTM grid.  Earlier
    versions incorrectly treated the EPSG:3857 tile data as equirectangular
    WGS84, causing systematic north-south misalignment (proportional to the
    Mercator stretch at the site's latitude, ~12% at 32°N, ~29% at 39°N).

    Cache version 2 stores only correctly reprojected basemaps.  Any cached
    file from version 1 is discarded and regenerated on the next call.

    Returns ``(rgb, provider_name, from_cache)`` where *rgb* is ``(rows, cols, 3)``
    uint8, or ``None`` when tiles cannot be retrieved.
    """
    import rasterio as rio
    from rasterio.warp import transform_bounds

    basemap_dir = paths.repo_root / "data" / "basemap" / site.id
    rgb_path = basemap_dir / "rgb.tif"
    meta_path = basemap_dir / "rgb.meta.json"

    dst_rows, dst_cols = hs_shape

    # Load from cache only when version matches and shape is consistent.
    if rgb_path.is_file() and meta_path.is_file():
        try:
            meta = json.loads(meta_path.read_text())
            if meta.get("cache_version") == FIG03_BASEMAP_CACHE_VERSION:
                with rio.open(rgb_path) as ds:
                    if (ds.height, ds.width) == (dst_rows, dst_cols):
                        rgb = np.transpose(ds.read(), (1, 2, 0))
                        provider = meta.get("provider", "cached")
                        return rgb.astype(np.uint8), provider, True
        except Exception:
            pass
        # Version mismatch or shape mismatch — discard stale cache.
        rgb_path.unlink(missing_ok=True)
    elif rgb_path.is_file():
        # Meta missing — can't verify version; discard.
        rgb_path.unlink(missing_ok=True)

    # Compute WGS84 extent of the hillshade grid, add a small margin so tile
    # edges don't clip the figure corners, then convert to EPSG:3857.
    west, south, east, north = _hillshade_bounds_wgs84(hs_transform, hs_shape, crs)
    pad = max((east - west), (north - south)) * 0.02
    west, south = west - pad, south - pad
    east, north = east + pad, north + pad

    # Convert padded WGS84 bounds to EPSG:3857 for contextily.
    w_3857, s_3857, e_3857, n_3857 = transform_bounds(
        "EPSG:4326", "EPSG:3857", west, south, east, north
    )

    fetched = _fetch_tiles(w_3857, s_3857, e_3857, n_3857, zoom)
    if fetched is None:
        return None

    rgb_tiles, provider, extent_3857 = fetched
    # Reproject from EPSG:3857 (actual tile CRS) to the hillshade UTM grid.
    rgb = _reproject_rgb_to_grid(
        rgb_tiles, extent_3857, "EPSG:3857",
        hs_transform, crs, hs_shape,
    )

    basemap_dir.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff",
        "dtype": "uint8",
        "count": 3,
        "height": dst_rows,
        "width": dst_cols,
        "crs": crs,
        "transform": hs_transform,
    }
    with rio.open(rgb_path, "w", **profile) as dst:
        for b in range(3):
            dst.write(rgb[:, :, b], b + 1)

    meta = {
        "provider": provider,
        "zoom": zoom,
        "cache_version": FIG03_BASEMAP_CACHE_VERSION,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "bounds_wgs84": [west, south, east, north],
        "bounds_3857": list(extent_3857),
    }
    meta_path.write_text(json.dumps(meta, indent=2))

    return rgb, provider, False
