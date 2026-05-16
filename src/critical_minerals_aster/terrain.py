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

    if not dem_reprojected.is_file():
        dem_dir.mkdir(parents=True, exist_ok=True)
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

                hgt_files = list(tmppath.glob("**/*.hgt")) + list(
                    tmppath.glob("**/*dem.tif")
                )
                if not hgt_files:
                    print(
                        f"  [terrain] No .hgt or dem.tif files after download for {site.id}",
                        file=sys.stderr,
                    )
                    return None

                if len(hgt_files) == 1:
                    import shutil

                    shutil.copy(hgt_files[0], dem_merged)
                else:
                    datasets = [rio.open(f) for f in hgt_files]
                    try:
                        merged, merged_t = rmerge(datasets)
                        profile = datasets[0].profile.copy()
                        profile.update(
                            height=merged.shape[1],
                            width=merged.shape[2],
                            transform=merged_t,
                            count=1,
                        )
                        with rio.open(dem_merged, "w", **profile) as dst:
                            dst.write(merged[0], 1)
                    finally:
                        for ds in datasets:
                            ds.close()

        except Exception as exc:
            print(
                f"  [terrain] DEM download failed for {site.id}: {exc}",
                file=sys.stderr,
            )
            return None

        # Reproject / resample to match TIR raster pixel grid exactly.
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
