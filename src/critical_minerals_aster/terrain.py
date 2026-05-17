"""DEM download and hillshade computation for structural geology context."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import rasterio
    from critical_minerals_aster.config import SiteConfig
    from critical_minerals_aster.paths import SitePaths


def compute_hillshade_for_site(
    site: SiteConfig,
    paths: SitePaths,
    transform: rasterio.Affine,
    shape: tuple[int, int],
    crs: rasterio.crs.CRS,
) -> "tuple[np.ndarray, rasterio.Affine, tuple[int, int]] | None":
    """Download NASADEM and compute hillshade covering the full site.bbox_wgs84.

    Returns ``(hillshade, hs_transform, hs_shape)`` where:
    - *hillshade* is a float32 array (NaN = nodata/transparent)
    - *hs_transform* is the affine transform of the hillshade grid
    - *hs_shape* is ``(rows, cols)`` of the hillshade grid

    The hillshade grid covers ``site.bbox_wgs84`` (not just the ASTER footprint)
    at the same pixel resolution as *transform*, ensuring terrain context even
    where ASTER has swath-rotation coverage gaps.

    Returns ``None`` if the DEM cannot be downloaded or computed.

    Cache invalidation
    ------------------
    ``dem_merged.tif``      — raw NASADEM tiles merged in WGS84; re-downloaded
                              when the site bbox changes (coverage check).
    ``dem_reprojected.tif`` — DEM reprojected to the hillshade grid; invalidated
                              when the pixel resolution or site bbox changes.
    """
    import tempfile

    import rasterio as rio
    import scipy.ndimage
    from rasterio.merge import merge as rmerge
    from rasterio.transform import from_bounds as transform_from_bounds
    from rasterio.warp import Resampling, reproject
    from rasterio.warp import transform_bounds as warp_transform_bounds

    dem_dir: Path = paths.repo_root / "data" / "dem" / site.id
    dem_merged = dem_dir / "dem_merged.tif"
    dem_reprojected = dem_dir / "dem_reprojected.tif"

    # Target coverage: full site bbox (not just ASTER raster footprint) so the
    # hillshade fills the figure even where ASTER has swath-rotation gaps.
    # Add 0.1° buffer so NASADEM tile boundaries never clip the visible area.
    _buf = 0.10
    west, south, east, north = site.bbox_wgs84
    dl_bbox = (west - _buf, south - _buf, east + _buf, north + _buf)

    # Pixel size from the ASTER TIR raster (metres per pixel in CRS units).
    px_size = abs(transform.a)  # assumes square pixels

    # Compute the hillshade grid: same CRS as ASTER, covers site.bbox_wgs84
    # exactly (no buffer) at the same resolution.  This is the canonical output
    # grid — dem_reprojected.tif is always built to match it.
    try:
        hs_bounds_utm = warp_transform_bounds("EPSG:4326", crs, west, south, east, north)
    except Exception as exc:
        print(f"  [terrain] bbox projection failed for {site.id}: {exc}", file=sys.stderr)
        return None
    hs_left, hs_bottom, hs_right, hs_top = hs_bounds_utm
    hs_cols = max(1, int(round((hs_right - hs_left) / px_size)))
    hs_rows = max(1, int(round((hs_top - hs_bottom) / px_size)))
    hs_transform = transform_from_bounds(hs_left, hs_bottom, hs_right, hs_top, hs_cols, hs_rows)

    # ---- Cache invalidation -------------------------------------------------
    # dem_reprojected: invalidate when pixel size or grid shape changed.
    if dem_reprojected.is_file():
        try:
            with rio.open(dem_reprojected) as _ds:
                cached_px = abs(_ds.transform.a)
                shape_ok = (_ds.height, _ds.width) == (hs_rows, hs_cols)
                res_ok = abs(cached_px - px_size) / px_size < 0.01
                if not (shape_ok and res_ok):
                    print(
                        f"  [terrain] DEM cache stale for {site.id} "
                        f"(cached {_ds.height}×{_ds.width}@{cached_px:.0f}m "
                        f"vs expected {hs_rows}×{hs_cols}@{px_size:.0f}m); regenerating.",
                        file=sys.stderr,
                    )
                    dem_reprojected.unlink()
        except Exception:
            dem_reprojected.unlink(missing_ok=True)

    # dem_merged: keep across runs, but re-download if corrupt or too small.
    if dem_merged.is_file() and not dem_reprojected.is_file():
        try:
            with rio.open(dem_merged):
                pass
        except Exception:
            dem_merged.unlink(missing_ok=True)

    if dem_merged.is_file() and not dem_reprojected.is_file():
        try:
            from shapely.geometry import box as _box
            req_box = _box(west, south, east, north)
            with rio.open(dem_merged) as _ds:
                m_bounds = warp_transform_bounds(_ds.crs, "EPSG:4326", *_ds.bounds)
            merged_box = _box(*m_bounds)
            overlap = req_box.intersection(merged_box).area / req_box.area
            if overlap < 0.95:
                print(
                    f"  [terrain] dem_merged covers only {overlap:.0%} of site bbox "
                    f"for {site.id}; re-downloading.",
                    file=sys.stderr,
                )
                dem_merged.unlink(missing_ok=True)
        except Exception:
            pass

    # ---- Download -----------------------------------------------------------
    if not dem_reprojected.is_file():
        dem_dir.mkdir(parents=True, exist_ok=True)

        if not dem_merged.is_file():
            try:
                import earthaccess

                earthaccess.login(strategy="netrc")
                results = earthaccess.search_data(
                    short_name="NASADEM_HGT",
                    bounding_box=dl_bbox,  # buffered bbox catches all edge tiles
                    count=20,
                )
                if not results:
                    print(f"  [terrain] No NASADEM tiles found for {site.id}", file=sys.stderr)
                    return None

                with tempfile.TemporaryDirectory() as tmpdir:
                    earthaccess.download(results, tmpdir)
                    tmppath = Path(tmpdir)

                    import zipfile
                    for zf in tmppath.glob("**/*.zip"):
                        try:
                            with zipfile.ZipFile(zf) as z:
                                z.extractall(zf.parent)
                        except Exception:
                            pass

                    candidates = (
                        list(tmppath.glob("**/*.hgt"))
                        + list(tmppath.glob("**/*.tif"))
                        + list(tmppath.glob("**/*dem.tif"))
                    )
                    valid_files: list[Path] = []
                    for f in candidates:
                        try:
                            with rio.open(f) as ds:
                                if ds.count >= 1 and ds.width > 0 and ds.height > 0:
                                    valid_files.append(f)
                        except Exception:
                            pass

                    if not valid_files:
                        print(f"  [terrain] No readable DEM files for {site.id}", file=sys.stderr)
                        return None

                    # Merge tiles, clipping to the buffered bbox.
                    datasets = [rio.open(f) for f in valid_files]
                    try:
                        merged, merged_t = rmerge(datasets, bounds=dl_bbox)
                        profile = datasets[0].profile.copy()
                        profile.update(
                            height=merged.shape[1],
                            width=merged.shape[2],
                            transform=merged_t,
                            count=1,
                            driver="GTiff",
                            dtype="float32",
                        )
                        with rio.open(dem_merged, "w", **profile) as dst:
                            dst.write(merged[0].astype("float32"), 1)
                    finally:
                        for ds in datasets:
                            ds.close()

            except Exception as exc:
                print(f"  [terrain] DEM download failed for {site.id}: {exc}", file=sys.stderr)
                return None

        # ---- Reproject to hillshade grid ------------------------------------
        try:
            with rio.open(dem_merged) as src:
                profile = src.profile.copy()
                profile.update(
                    crs=crs,
                    transform=hs_transform,
                    width=hs_cols,
                    height=hs_rows,
                    dtype="float32",
                    count=1,
                    nodata=-9999.0,
                )
                with rio.open(dem_reprojected, "w", **profile) as dst:
                    reproject(
                        source=rio.band(src, 1),
                        destination=rio.band(dst, 1),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=hs_transform,
                        dst_crs=crs,
                        resampling=Resampling.bilinear,
                    )
        except Exception as exc:
            print(f"  [terrain] DEM reprojection failed for {site.id}: {exc}", file=sys.stderr)
            return None

    # ---- Load and compute hillshade -----------------------------------------
    try:
        with rio.open(dem_reprojected) as src:
            dem = src.read(1).astype(float)
            # Read actual transform from cached file (may differ from hs_transform
            # if the cache predates this run's pixel-exact grid).
            hs_transform = src.transform
            hs_rows, hs_cols = src.height, src.width
    except Exception as exc:
        print(f"  [terrain] Could not read reprojected DEM for {site.id}: {exc}", file=sys.stderr)
        return None

    dem[dem < -9000] = np.nan  # NASADEM nodata sentinel

    dzdx = scipy.ndimage.uniform_filter1d(np.gradient(dem, axis=1), 3)
    dzdy = scipy.ndimage.uniform_filter1d(np.gradient(dem, axis=0), 3)

    az = np.radians(315)   # sun azimuth: NW (standard cartographic)
    el = np.radians(45)    # sun elevation: 45°
    slope = np.arctan(np.sqrt(dzdx**2 + dzdy**2))
    aspect = np.arctan2(-dzdx, dzdy)
    hs = np.cos(el) * np.cos(slope) + np.sin(el) * np.sin(slope) * np.cos(az - aspect)
    hs = np.clip(hs, 0, 1)
    hs[~np.isfinite(hs)] = np.nan

    return hs.astype(np.float32), hs_transform, (hs_rows, hs_cols)
