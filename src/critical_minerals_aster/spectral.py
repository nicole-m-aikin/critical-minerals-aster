"""ASTER TIR band I/O, granule selection, and alteration ratios."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import rasterio
from rasterio.profiles import Profile
from shapely.geometry import Polygon, box

from critical_minerals_aster.config import BBox

REQUIRED_TIR_BANDS = (10, 11, 12, 13, 14)
_GRANULE_ID_RE = re.compile(r"(AST_L1T_\d+_\d+)")


def extract_granule_id(granule: Any) -> str:
    """Parse AST_L1T granule id prefix from earthaccess data links."""
    for link in granule.data_links():
        match = _GRANULE_ID_RE.search(str(link))
        if match:
            return match.group(1)
    raise ValueError("Could not determine granule id from data links")


def _granule_footprint(granule: Any) -> Polygon:
    """Footprint polygon from CMR UMM metadata (lon/lat order)."""
    points = granule["umm"]["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"][
        "GPolygons"
    ][0]["Boundary"]["Points"]
    lons = [float(p["Longitude"]) for p in points]
    lats = [float(p["Latitude"]) for p in points]
    poly = Polygon(zip(lons, lats))
    if not poly.is_valid:
        poly = poly.buffer(0)
    return poly


def _bbox_coverage_fraction(bbox_wgs84: BBox, footprint: Polygon) -> float:
    """Fraction of study bbox area intersecting the granule footprint (0–1)."""
    lon0, lat0, lon1, lat1 = bbox_wgs84
    study = box(lon0, lat0, lon1, lat1)
    if study.area == 0:
        return 0.0
    return float(study.intersection(footprint).area / study.area)


def _tir_band_count(granule: Any, bands: Sequence[int] = REQUIRED_TIR_BANDS) -> int:
    """Count TIR band files present in granule data links."""
    urls = [str(u).upper() for u in granule.data_links()]
    return sum(
        1 for band in bands if any(f"TIR_B{band}" in url for url in urls)
    )


def score_granule(granule: Any, bbox_wgs84: BBox) -> tuple[float, float, int]:
    """
    Score one granule for ranking.

    Returns (coverage_fraction, composite_score, tir_band_count).
    Composite favors bbox coverage, then band completeness.
    """
    footprint = _granule_footprint(granule)
    coverage = _bbox_coverage_fraction(bbox_wgs84, footprint)
    band_count = _tir_band_count(granule)
    composite = coverage * 10.0 + band_count
    return coverage, composite, band_count


def select_granule(
    results: Sequence[Any],
    bbox_wgs84: BBox,
    granule_id_override: str | None = None,
) -> Any:
    """
    Choose the best earthaccess ASTER L1T granule for a study bbox.

    If ``granule_id_override`` is set (non-empty), return the matching granule
    from ``results``. Otherwise rank by bbox coverage and TIR B10–B14 availability.
    """
    if not results:
        raise ValueError("No granules in search results")

    if granule_id_override:
        for granule in results:
            granule_id = extract_granule_id(granule)
            if granule_id == granule_id_override or granule_id_override in granule_id:
                return granule
        raise ValueError(
            f"granule_id override {granule_id_override!r} not found in search results"
        )

    best: Any | None = None
    best_score = -1.0
    for granule in results:
        try:
            band_count = _tir_band_count(granule)
        except (AttributeError, TypeError):
            continue
        if band_count == 0:
            continue
        try:
            _, composite, _ = score_granule(granule, bbox_wgs84)
        except (KeyError, IndexError, TypeError, ValueError):
            continue
        if composite > best_score:
            best_score = composite
            best = granule

    if best is None:
        raise ValueError("No granule with TIR bands found in search results")
    return best


def clip_bands_to_bbox(
    bands: list[np.ndarray],
    transform: rasterio.Affine,
    crs: rasterio.crs.CRS,
    bbox_wgs84: BBox,
) -> tuple[list[np.ndarray], rasterio.Affine]:
    """Clip a list of co-registered arrays to a WGS84 bounding box.

    Reprojects the bbox to the raster CRS, computes the pixel window, and
    returns sliced arrays with an updated affine transform.  Any band pixels
    outside the bbox are not included in downstream percentile statistics or
    zone vectorisation, making results site-specific rather than whole-scene.

    Returns (clipped_bands, new_transform).  If the bbox window is empty or
    entirely outside the raster the original arrays and transform are returned
    unchanged so the pipeline degrades gracefully.
    """
    from rasterio.crs import CRS as RasterioCRS
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds, Window

    if not bands:
        return bands, transform

    rows, cols = bands[0].shape
    try:
        dst_crs = RasterioCRS.from_epsg(4326)
        lon0, lat0, lon1, lat1 = bbox_wgs84
        # transform bbox from WGS84 into the raster CRS
        x0, y0, x1, y1 = transform_bounds(dst_crs, crs, lon0, lat0, lon1, lat1)
        window = from_bounds(x0, y0, x1, y1, transform)
        # clamp to actual raster extent
        row_off = max(0, int(window.row_off))
        col_off = max(0, int(window.col_off))
        row_end = min(rows, int(window.row_off + window.height))
        col_end = min(cols, int(window.col_off + window.width))
        if row_end <= row_off or col_end <= col_off:
            return bands, transform
        clipped = [b[row_off:row_end, col_off:col_end] for b in bands]
        new_transform = transform * rasterio.Affine.translation(col_off, row_off)
        return clipped, new_transform
    except Exception:
        return bands, transform


def band_ratio(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(b != 0, a / b, np.nan)


def load_tir_band(
    aster_dir: os.PathLike[str] | str,
    granule_id: str,
    band_num: int,
) -> tuple[np.ndarray, Profile, rasterio.Affine, rasterio.crs.CRS]:
    """Load one ASTER L1T TIR band (float32), masking zero to NaN."""
    path = Path(aster_dir) / f"{granule_id}_TIR_B{band_num}.tif"
    with rasterio.open(path) as src:
        data = src.read(1).astype(float)
        profile = src.profile
        transform = src.transform
        crs = src.crs
    data[data == 0] = np.nan
    return data, profile, transform, crs


def load_tir_bands_10_14(
    aster_dir: os.PathLike[str] | str,
    granule_id: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    Profile,
    rasterio.Affine,
    rasterio.crs.CRS,
]:
    """Load TIR bands B10–B14; profile/transform/crs come from B10."""
    b10, profile, transform, crs = load_tir_band(aster_dir, granule_id, 10)
    b11, _, _, _ = load_tir_band(aster_dir, granule_id, 11)
    b12, _, _, _ = load_tir_band(aster_dir, granule_id, 12)
    b13, _, _, _ = load_tir_band(aster_dir, granule_id, 13)
    b14, _, _, _ = load_tir_band(aster_dir, granule_id, 14)
    return b10, b11, b12, b13, b14, profile, transform, crs


def alteration_ratios(
    b12: np.ndarray, b13: np.ndarray, b14: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Silica, carbonate, and mafic ratios (same definitions as notebooks)."""
    silica = band_ratio(b13, b14)
    carbonate = band_ratio(b13, b12)
    mafic = band_ratio(b12, b13)
    return silica, carbonate, mafic
