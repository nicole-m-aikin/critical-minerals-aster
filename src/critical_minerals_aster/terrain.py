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
) -> np.ndarray | None:
    """Download NASADEM and compute hillshade matching the TIR raster extent.

    Returns a (rows, cols) uint8 hillshade array, or None if DEM unavailable.
    The result is cached under data/dem/{site_id}/ so subsequent runs skip the
    download and reprojection steps.

    Cache invalidation
    ------------------
    ``dem_merged.tif``    — raw NASADEM tiles merged to a single WGS84 GeoTIFF.
                            Kept across runs; only re-downloaded if missing.
    ``dem_reprojected.tif`` — DEM reprojected to the TIR raster pixel grid.
                            Invalidated whenever the TIR shape changes (e.g.
                            after switching from a single granule to a mosaic).
    """
    import tempfile

    import rasterio as rio
    import scipy.ndimage
    from rasterio.merge import merge as rmerge
    from rasterio.warp import Resampling, reproject

    from critical_minerals_aster.spectral import raster_bbox_wgs84

    bbox = raster_bbox_wgs84(transform, shape, crs)
    dem_dir: Path = paths.repo_root / "data" / "dem" / site.id
    dem_merged = dem_dir / "dem_merged.tif"
    dem_reprojected = dem_dir / "dem_reprojected.tif"

    # Invalidate the cached reprojection when the TIR raster dimensions have
    # changed — e.g. after switching from a single granule to a mosaic, or
    # after modifying the site bbox.  dem_merged.tif is kept so only the
    # cheap reprojection step needs to re-run, not the full DEM download.
    if dem_reprojected.is_file():
        try:
            with rio.open(dem_reprojected) as _ds:
                if (_ds.height, _ds.width) != shape:
                    print(
                        f"  [terrain] DEM cache shape mismatch for {site.id} "
                        f"(cached {_ds.height}×{_ds.width} vs expected "
                        f"{shape[0]}×{shape[1]}); regenerating reprojection.",
                        file=sys.stderr,
                    )
                    dem_reprojected.unlink()
        except Exception:
            dem_reprojected.unlink(missing_ok=True)

    # Purge a corrupt cached merge (e.g. an HGT binary that was wrongly
    # saved with a .tif extension on a previous run).
    if dem_merged.is_file() and not dem_reprojected.is_file():
        try:
            with rio.open(dem_merged):
                pass
        except Exception:
            dem_merged.unlink(missing_ok=True)

    # Validate that dem_merged covers the required bbox sufficiently.
    # A stale merge built from a different granule/mosaic may cover only a tiny
    # fraction of the current raster extent — reprojecting it would fill most of
    # the output with nodata and produce a flat-grey hillshade.
    # Threshold: the merged DEM must overlap at least 80 % of the required bbox.
    if dem_merged.is_file() and not dem_reprojected.is_file():
        try:
            from rasterio.warp import transform_bounds as _tb
            from shapely.geometry import box as _box
            lon0_req, lat0_req, lon1_req, lat1_req = bbox
            with rio.open(dem_merged) as _ds:
                lon0_m, lat0_m, lon1_m, lat1_m = _tb(
                    _ds.crs, "EPSG:4326", *_ds.bounds
                )
            req_box = _box(lon0_req, lat0_req, lon1_req, lat1_req)
            merged_box = _box(lon0_m, lat0_m, lon1_m, lat1_m)
            overlap_frac = req_box.intersection(merged_box).area / req_box.area
            if overlap_frac < 0.80:
                print(
                    f"  [terrain] dem_merged coverage too small for {site.id} "
                    f"({overlap_frac:.0%} of required bbox); re-downloading.",
                    file=sys.stderr,
                )
                dem_merged.unlink(missing_ok=True)
        except Exception:
            pass  # if check fails, keep the file and let reprojection proceed

    if not dem_reprojected.is_file():
        dem_dir.mkdir(parents=True, exist_ok=True)

        if not dem_merged.is_file():
            # No merged DEM on disk yet — download NASADEM tiles.
            try:
                import earthaccess

                earthaccess.login(strategy="netrc")
                results = earthaccess.search_data(
                    short_name="NASADEM_HGT",
                    bounding_box=bbox,
                    count=20,
                )
                if not results:
                    print(
                        f"  [terrain] No NASADEM tiles found for {site.id}",
                        file=sys.stderr,
                    )
                    return None

                with tempfile.TemporaryDirectory() as tmpdir:
                    earthaccess.download(results, tmpdir)
                    tmppath = Path(tmpdir)

                    # NASADEM_HGT tiles are delivered as .zip archives containing
                    # a single .hgt elevation file.  Extract all zips first.
                    import zipfile

                    for zf in tmppath.glob("**/*.zip"):
                        try:
                            with zipfile.ZipFile(zf) as z:
                                z.extractall(zf.parent)
                        except Exception:
                            pass

                    # Collect every candidate file and validate with rasterio
                    # (avoids misreading .num/.lsm auxiliary files or empty zips).
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
                        print(
                            f"  [terrain] No readable DEM files after download for {site.id}",
                            file=sys.stderr,
                        )
                        return None

                    # Merge (clipping to site bbox so the output stays small even
                    # when many 1°×1° NASADEM tiles are required).
                    lon0, lat0, lon1, lat1 = bbox
                    datasets = [rio.open(f) for f in valid_files]
                    try:
                        merged, merged_t = rmerge(
                            datasets, bounds=(lon0, lat0, lon1, lat1)
                        )
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
                print(
                    f"  [terrain] DEM download failed for {site.id}: {exc}",
                    file=sys.stderr,
                )
                return None

        # Reproject / resample dem_merged to match TIR raster pixel grid exactly.
        # This runs whether dem_merged was just downloaded or already existed.
        rows, cols = shape
        try:
            with rio.open(dem_merged) as src:
                profile = src.profile.copy()
                profile.update(
                    crs=crs,
                    transform=transform,
                    width=cols,
                    height=rows,
                    dtype="float32",
                    count=1,
                )
                with rio.open(dem_reprojected, "w", **profile) as dst:
                    reproject(
                        source=rio.band(src, 1),
                        destination=rio.band(dst, 1),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=crs,
                        resampling=Resampling.bilinear,
                    )
        except Exception as exc:
            print(
                f"  [terrain] DEM reprojection failed for {site.id}: {exc}",
                file=sys.stderr,
            )
            return None

    # Load reprojected DEM and compute hillshade.
    try:
        with rio.open(dem_reprojected) as src:
            dem = src.read(1).astype(float)
    except Exception as exc:
        print(
            f"  [terrain] Could not read reprojected DEM for {site.id}: {exc}",
            file=sys.stderr,
        )
        return None

    dem[dem < -9000] = np.nan  # nodata mask (NASADEM uses -9999)

    # Smooth gradients slightly to reduce tile-boundary noise.
    dzdx = scipy.ndimage.uniform_filter1d(np.gradient(dem, axis=1), 3)
    dzdy = scipy.ndimage.uniform_filter1d(np.gradient(dem, axis=0), 3)

    # Sun azimuth 315° (NW), elevation 45° — standard cartographic convention.
    az = np.radians(315)
    el = np.radians(45)
    slope = np.arctan(np.sqrt(dzdx**2 + dzdy**2))
    aspect = np.arctan2(-dzdx, dzdy)
    hs = np.cos(el) * np.cos(slope) + np.sin(el) * np.sin(slope) * np.cos(az - aspect)
    hs = np.clip(hs, 0, 1)
    hs[~np.isfinite(hs)] = 0.5  # fill nodata pixels with neutral grey

    return (hs * 255).astype(np.uint8)
