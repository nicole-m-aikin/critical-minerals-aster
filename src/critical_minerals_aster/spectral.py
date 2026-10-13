"""ASTER TIR band I/O, granule selection, and alteration ratios."""

from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import rasterio
from rasterio.profiles import Profile
from shapely.geometry import Polygon, box

from critical_minerals_aster.config import BBox

REQUIRED_TIR_BANDS = (10, 11, 12, 13, 14)
DEFAULT_MAX_CLOUD_PCT = 30.0
CLOUD_SCORE_WEIGHT = 5.0
# May–September daytime scenes, as in the Nevada ASTER TIR mosaic of
# Rockwell & Hofstra (2008), limit snow and seasonal thermal contrast.
PREFERRED_MONTHS: tuple[int, ...] = (5, 6, 7, 8, 9)
_GRANULE_ID_RE = re.compile(r"(AST_(?:L1T|05|08)_\d+_\d+)")
_ACQUISITION_RE = re.compile(r"AST_(?:L1T|05|08)_(\d{17})_")

# AST_05 v004: uint16 emissivity x 0.001, 0 = fill (LP DAAC catalogue).
EMISSIVITY_SCALE = 0.001
EMISSIVITY_VALID_RANGE = (0.5, 1.0)


def extract_granule_id(granule: Any) -> str:
    """Parse the AST_L1T / AST_05 / AST_08 granule id prefix from earthaccess data links."""
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
    """Count TIR band files (L1T ``TIR_B10`` or AST_05 ``SRE_TIR_B10``) in granule data links."""
    urls = [str(u).upper() for u in granule.data_links()]
    return sum(
        1 for band in bands if any(f"TIR_B{band}." in url or f"TIR_B{band}_" in url for url in urls)
    )


def granule_cloud_cover(granule: Any) -> float | None:
    """Scene cloud cover (0–100 %) from CMR UMM ``CloudCover``; None if absent or invalid."""
    try:
        value = granule["umm"].get("CloudCover")
    except (KeyError, TypeError, AttributeError):
        return None
    if value is None:
        return None
    try:
        pct = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(pct) or pct < 0 or pct > 100:
        return None
    return pct


def granule_month(granule: Any) -> int | None:
    """Acquisition month (1–12) from CMR UMM ``TemporalExtent``; None if absent."""
    try:
        temporal = granule["umm"].get("TemporalExtent", {})
        stamp = temporal.get("SingleDateTime") or temporal.get("RangeDateTime", {}).get(
            "BeginningDateTime"
        )
        return int(str(stamp)[5:7])
    except (KeyError, TypeError, AttributeError, ValueError):
        return None


def granule_day_night(granule: Any) -> str | None:
    """``"Day"`` or ``"Night"`` from CMR UMM ``DataGranule.DayNightFlag``; None otherwise."""
    try:
        flag = granule["umm"].get("DataGranule", {}).get("DayNightFlag")
    except (KeyError, TypeError, AttributeError):
        return None
    if not isinstance(flag, str):
        return None
    flag = flag.strip().capitalize()
    return flag if flag in ("Day", "Night") else None


def _selection_tier(
    cloud: float | None, day_night: str | None, max_cloud_pct: float | None
) -> int:
    """0 clear day, 1 clear night, 2 cloudy day, 3 cloudy night (unknown = day / clear).

    Radiance ratios are temperature-sensitive, so a clear night scene is preferred
    over a cloudy day scene, but day is preferred when both are clear.
    """
    cloudy = max_cloud_pct is not None and cloud is not None and cloud > max_cloud_pct
    night = day_night == "Night"
    return 2 * int(cloudy) + int(night)


_TIER_NAMES = ("day_clear", "night_clear", "day_cloudy", "night_cloudy")


def score_granule(
    granule: Any,
    bbox_wgs84: BBox,
    cloud_weight: float = CLOUD_SCORE_WEIGHT,
) -> tuple[float, float, int]:
    """
    Score one granule for ranking.

    Returns (coverage_fraction, composite_score, tir_band_count).
    Composite favors bbox coverage, then band completeness, minus a cloud
    penalty of ``cloud_weight * cloud_pct / 100`` (with the default weight,
    10 points of scene cloud cost as much as 5 points of bbox coverage).
    Unknown cloud cover is not penalized.
    """
    footprint = _granule_footprint(granule)
    coverage = _bbox_coverage_fraction(bbox_wgs84, footprint)
    band_count = _tir_band_count(granule)
    cloud = granule_cloud_cover(granule)
    cloud_penalty = cloud_weight * cloud / 100.0 if cloud is not None else 0.0
    composite = coverage * 10.0 + band_count - cloud_penalty
    return coverage, composite, band_count


def select_granule(
    results: Sequence[Any],
    bbox_wgs84: BBox,
    granule_id_override: str | None = None,
    max_cloud_pct: float | None = DEFAULT_MAX_CLOUD_PCT,
    prefer_day: bool = True,
) -> Any:
    """
    Choose the best earthaccess ASTER L1T granule for a study bbox.

    If ``granule_id_override`` is set (non-empty), return the matching granule
    from ``results`` regardless of cloud cover. Otherwise granules are ranked
    first by tier (clear day, clear night, cloudy day, cloudy night; "cloudy"
    means scene cloud above ``max_cloud_pct``), then by the cloud-penalized
    composite from ``score_granule``. Unknown cloud counts as clear and unknown
    day/night as day. Falling back past clear-day emits a warning.
    ``max_cloud_pct=None`` disables the cloud cutoff; ``prefer_day=False``
    ignores day/night.
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
    best_key: tuple[int, float] = (99, np.inf)
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
        day_night = granule_day_night(granule) if prefer_day else None
        tier = _selection_tier(granule_cloud_cover(granule), day_night, max_cloud_pct)
        key = (tier, -composite)
        if key < best_key:
            best_key = key
            best = granule

    if best is None:
        raise ValueError("No granule with TIR bands found in search results")
    if best_key[0] > 0:
        cloud = granule_cloud_cover(best)
        cloud_txt = f"{cloud:g}% cloud" if cloud is not None else "unknown cloud"
        warnings.warn(
            f"No clear-sky day granule (<= {max_cloud_pct:g}% cloud); falling back to "
            f"{extract_granule_id(best)} ({_TIER_NAMES[best_key[0]]}, {cloud_txt})",
            stacklevel=2,
        )
    return best


@dataclass
class MosaicSelection:
    """Granules chosen for a site mosaic, best (reference) first."""

    granules: list[Any]
    tier: str
    union_coverage: float
    records: list[dict[str, Any]] = field(default_factory=list)


def select_mosaic_granules(
    results: Sequence[Any],
    bbox_wgs84: BBox,
    *,
    max_cloud_pct: float | None = DEFAULT_MAX_CLOUD_PCT,
    min_coverage: float = 0.30,
    target_coverage: float = 0.98,
    min_acceptable_coverage: float = 0.90,
    fallback_coverage_tolerance: float = 0.05,
    min_gain: float = 0.02,
    max_granules: int = 6,
    prefer_day: bool = True,
    max_size_mb: float | None = None,
    preferred_months: tuple[int, ...] | None = PREFERRED_MONTHS,
) -> MosaicSelection:
    """
    Choose a small, single-tier set of granules that covers the study bbox.

    Tiers are tried in order (clear day, clear night, cloudy day, cloudy night).
    Day and night scenes are never mixed in one mosaic. Within each tier,
    ``preferred_months`` scenes are tried alone before all months (tier name
    gets an ``_offseason`` suffix if any off-season granule is used). Within a tier the
    reference granule maximises ``coverage - 0.5 * cloud/100``; further
    granules are added greedily by new bbox area gained (same cloud penalty)
    until ``target_coverage`` or ``max_granules`` is reached. The first tier
    whose union reaches ``min_acceptable_coverage`` wins; otherwise the earliest
    tier within ``fallback_coverage_tolerance`` of the largest union is used
    (some bboxes sit on a swath edge that no day pass fully covers) and a
    warning is emitted.
    """
    lon0, lat0, lon1, lat1 = bbox_wgs84
    study = box(lon0, lat0, lon1, lat1)
    if study.area == 0:
        raise ValueError("Degenerate bbox")

    by_tier: dict[int, list[dict[str, Any]]] = {i: [] for i in range(4)}
    for granule in results:
        try:
            if _tir_band_count(granule) < 3:
                continue
            if max_size_mb is not None and granule.size() >= max_size_mb:
                continue
            footprint = _granule_footprint(granule)
        except Exception:
            continue
        clipped = study.intersection(footprint)
        coverage = float(clipped.area / study.area)
        if coverage <= min_coverage:
            continue
        cloud = granule_cloud_cover(granule)
        day_night = granule_day_night(granule)
        tier = _selection_tier(cloud, day_night if prefer_day else None, max_cloud_pct)
        by_tier[tier].append({
            "granule": granule,
            "geom": clipped,
            "coverage": coverage,
            "cloud": cloud,
            "day_night": day_night,
            "month": granule_month(granule),
        })

    def _penalty(c: dict[str, Any]) -> float:
        return 0.5 * (c["cloud"] or 0.0) / 100.0

    def _greedy(cands: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], float]:
        if not cands:
            return [], 0.0
        chosen = [max(cands, key=lambda c: c["coverage"] - _penalty(c))]
        union = chosen[0]["geom"]
        remaining = [c for c in cands if c is not chosen[0]]
        while remaining and len(chosen) < max_granules and union.area / study.area < target_coverage:
            gains = [(c["geom"].difference(union).area / study.area, c) for c in remaining]
            gain, pick = max(gains, key=lambda gc: gc[0] - _penalty(gc[1]))
            if gain < min_gain:
                break
            chosen.append(pick)
            union = union.union(pick["geom"])
            remaining.remove(pick)
        return chosen, float(union.area / study.area)

    def _in_season(c: dict[str, Any]) -> bool:
        return preferred_months is None or c["month"] is None or c["month"] in preferred_months

    best: tuple[str, list[dict[str, Any]], float] | None = None
    fallbacks: list[tuple[str, list[dict[str, Any]], float]] = []
    for tier in range(4):
        pools = [[c for c in by_tier[tier] if _in_season(c)]]
        if preferred_months is not None:
            pools.append(by_tier[tier])
        for pool in pools:
            chosen, union_cov = _greedy(pool)
            if not chosen:
                continue
            name = _TIER_NAMES[tier]
            if not prefer_day:
                dn = {c["day_night"] for c in chosen}
                prefix = "night" if dn == {"Night"} else "day" if dn <= {"Day", None} else "mixed"
                name = f"{prefix}_{'cloudy' if tier >= 2 else 'clear'}"
            if not all(_in_season(c) for c in chosen):
                name += "_offseason"
            if union_cov >= min_acceptable_coverage:
                best = (name, chosen, union_cov)
                break
            fallbacks.append((name, chosen, union_cov))
        if best is not None:
            break
    if best is None and fallbacks:
        top = max(f[2] for f in fallbacks)
        best = next(f for f in fallbacks if f[2] >= top - fallback_coverage_tolerance)

    if best is None:
        return MosaicSelection(granules=[], tier="none", union_coverage=0.0)

    tier_name, chosen, union_cov = best
    if tier_name != "day_clear" or union_cov < min_acceptable_coverage:
        warnings.warn(
            f"Mosaic uses tier {tier_name} with {union_cov:.0%} bbox coverage",
            stacklevel=2,
        )
    records = [
        {
            "granule_id": extract_granule_id(c["granule"]),
            "cloud_pct": c["cloud"],
            "day_night": c["day_night"],
            "month": c["month"],
            "bbox_coverage": round(c["coverage"], 4),
        }
        for c in chosen
    ]
    return MosaicSelection(
        granules=[c["granule"] for c in chosen],
        tier=tier_name,
        union_coverage=union_cov,
        records=records,
    )


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


def raster_bbox_wgs84(
    transform: rasterio.Affine,
    shape: tuple[int, int],
    crs: rasterio.crs.CRS,
) -> BBox:
    """Return the WGS84 (lon0, lat0, lon1, lat1) extent of a raster.

    Derived from the affine transform and pixel dimensions so it correctly
    represents the area actually analysed by the pipeline, which may be a
    clipped subset of the original ASTER scene.
    """
    from rasterio.warp import transform_bounds

    rows, cols = shape
    x_min = transform.c
    y_max = transform.f
    x_max = transform.c + cols * transform.a
    y_min = transform.f + rows * transform.e  # transform.e is negative (north-up)
    lon0, lat0, lon1, lat1 = transform_bounds(crs, "EPSG:4326", x_min, y_min, x_max, y_max)
    return (lon0, lat0, lon1, lat1)


def band_ratio(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(b != 0, a / b, np.nan)


def load_classification_mosaic(
    aster_dir: os.PathLike[str] | str,
    granule_id: str,
) -> tuple[
    np.ndarray,  # silica_cls  (uint8: 0=bg, 1=moderate, 2=strong)
    np.ndarray,  # carbonate_cls
    np.ndarray,  # mafic_cls
    np.ndarray,  # valid_mask  (float: 1.0=valid, NaN=nodata)
    "rasterio.Affine",
    "rasterio.crs.CRS",
]:
    """Load per-granule classification mosaic files produced by _build_ratio_mosaic_from_paths.

    Classifications were applied per-granule before merging, so percentile
    thresholds reflect each granule's local scene statistics rather than the
    pooled mosaic distribution.  Returns classified uint8 arrays ready to pass
    directly to combined_score(), bypassing classify_percentiles().

    valid_mask has 1.0 where any granule had real data and NaN elsewhere;
    pass it to compute_valid_data_footprint() in place of B10.
    """
    aster_dir = Path(aster_dir)
    transform = crs = None
    classes: dict[str, np.ndarray] = {}
    for name in ("silica", "carbonate", "mafic"):
        path = aster_dir / f"{granule_id}_cls_{name}.tif"
        with rasterio.open(path) as src:
            arr = src.read(1).astype(np.uint8)
            if transform is None:
                transform = src.transform
                crs = src.crs
        classes[name] = arr

    valid_path = aster_dir / f"{granule_id}_cls_valid.tif"
    with rasterio.open(valid_path) as src:
        valid_mask = src.read(1).astype(float)
    valid_mask[valid_mask == 0] = np.nan

    return classes["silica"], classes["carbonate"], classes["mafic"], valid_mask, transform, crs


def load_ratio_mosaic(
    aster_dir: os.PathLike[str] | str,
    granule_id: str,
) -> tuple[
    np.ndarray,  # silica ratio
    np.ndarray,  # carbonate ratio
    np.ndarray,  # mafic ratio
    np.ndarray,  # B10 (for valid-pixel footprint only)
    "rasterio.Affine",
    "rasterio.crs.CRS",
]:
    """Load pre-computed ratio mosaic files produced by download_and_mosaic_aster.

    Returns ratio arrays (silica, carbonate, mafic) ready for classify_percentiles,
    plus B10 for footprint computation and the shared transform/CRS.
    Nodata (stored as 0 in the GeoTIFF) is restored to NaN.
    """
    aster_dir = Path(aster_dir)
    transform = crs = None
    ratios: dict[str, np.ndarray] = {}
    for name in ("silica", "carbonate", "mafic"):
        path = aster_dir / f"{granule_id}_ratio_{name}.tif"
        with rasterio.open(path) as src:
            arr = src.read(1).astype(float)
            if transform is None:
                transform = src.transform
                crs = src.crs
        arr[arr == 0] = np.nan
        ratios[name] = arr

    b10_path = aster_dir / f"{granule_id}_TIR_B10.tif"
    with rasterio.open(b10_path) as src:
        b10 = src.read(1).astype(float)
    b10[b10 == 0] = np.nan

    return ratios["silica"], ratios["carbonate"], ratios["mafic"], b10, transform, crs


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
    """Legacy radiance recipe, kept for GEE parity and archived comparisons.

    The labels are wrong by the TIR literature (B13/B14 is the carbonate index)
    and "carbonate" and "mafic" are exact reciprocals. Use ``alteration_indices``.
    """
    silica = band_ratio(b13, b14)
    carbonate = band_ratio(b13, b12)
    mafic = band_ratio(b12, b13)
    return silica, carbonate, mafic


def alteration_indices(
    b10: np.ndarray, b11: np.ndarray, b12: np.ndarray, b13: np.ndarray, b14: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Quartz, carbonate and mafic indices of Ninomiya (2003) / Ninomiya et al. (2005).

    QI = B11² / (B10·B12), CI = B13 / B14, MI = B12 / B13. Intended for
    emissivity (AST_05); on radiance they are temperature-dependent. Returned in
    the pipeline's (silica, carbonate, mafic) slot order.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = b10 * b12
        qi = np.where(denom != 0, b11 * b11 / denom, np.nan)
    ci = band_ratio(b13, b14)
    mi = band_ratio(b12, b13)
    return qi, ci, mi


# ASTER TIR has 10 detectors per band (whiskbroom). Mis-gain between detectors
# produces periodic along-track striping. Moment-matching equalizes each
# detector class to the scene mean/std (Iwasaki & Tonooka-style; see roadmap).
ASTER_TIR_N_DETECTORS = 10


def _detector_slices(length: int, detector: int, n_detectors: int) -> np.ndarray:
    return np.arange(detector, length, n_detectors)


def destripe_amplitude(
    arr: np.ndarray, n_detectors: int = ASTER_TIR_N_DETECTORS, axis: int = 0
) -> float:
    """Peak-to-peak difference of per-detector means (finite pixels only)."""
    if arr.ndim != 2:
        raise ValueError("destripe_amplitude expects a 2-D array")
    means: list[float] = []
    length = arr.shape[axis]
    for d in range(n_detectors):
        idx = _detector_slices(length, d, n_detectors)
        slab = arr[idx, :] if axis == 0 else arr[:, idx]
        if np.isfinite(slab).any():
            means.append(float(np.nanmean(slab)))
    if len(means) < 2:
        return 0.0
    return float(max(means) - min(means))


def destripe_scan_axis(
    arr: np.ndarray, n_detectors: int = ASTER_TIR_N_DETECTORS
) -> int:
    """Axis (0=rows, 1=cols) with stronger period-``n_detectors`` mean modulation."""
    a0 = destripe_amplitude(arr, n_detectors=n_detectors, axis=0)
    a1 = destripe_amplitude(arr, n_detectors=n_detectors, axis=1)
    return 0 if a0 >= a1 else 1


def destripe_moment_match(
    arr: np.ndarray,
    n_detectors: int = ASTER_TIR_N_DETECTORS,
    axis: int | None = None,
) -> np.ndarray:
    """Equalize mean and std of each detector class to the scene statistics.

    ``axis`` is the along-track / detector-cycle axis (0 = rows, 1 = columns).
    When ``axis`` is None, the stronger of the two is chosen. NaNs are preserved.
    Constant or all-NaN inputs are returned unchanged.
    """
    if arr.ndim != 2:
        raise ValueError("destripe_moment_match expects a 2-D array")
    out = np.array(arr, dtype=np.float64, copy=True)
    finite = np.isfinite(out)
    if not finite.any():
        return out
    global_vals = out[finite]
    g_mean = float(global_vals.mean())
    g_std = float(global_vals.std())
    if g_std <= 1e-12 * max(abs(g_mean), 1.0):
        return out
    if axis is None:
        axis = destripe_scan_axis(out, n_detectors=n_detectors)
    length = out.shape[axis]
    for d in range(n_detectors):
        idx = _detector_slices(length, d, n_detectors)
        if axis == 0:
            slab = out[idx, :]
            mask = finite[idx, :]
        else:
            slab = out[:, idx]
            mask = finite[:, idx]
        if not mask.any():
            continue
        vals = slab[mask]
        d_mean = float(vals.mean())
        d_std = float(vals.std())
        if d_std <= 1e-12 * max(abs(d_mean), 1.0):
            corrected = slab - d_mean + g_mean
        else:
            corrected = (slab - d_mean) / d_std * g_std + g_mean
        if axis == 0:
            out[idx, :] = np.where(mask, corrected, out[idx, :])
        else:
            out[:, idx] = np.where(mask, corrected, out[:, idx])
    return out


def destripe_emissivity_bands(
    bands: dict[int, np.ndarray],
    *,
    band_ids: Sequence[int] = (10, 11, 12),
    n_detectors: int = ASTER_TIR_N_DETECTORS,
    axis: int | None = None,
) -> tuple[dict[int, np.ndarray], dict[str, Any]]:
    """Destripe selected emissivity bands in place-safe copy; return stats.

    Axis is chosen from B11 (QI's squared numerator) when ``axis`` is None, then
    applied to every band so the indices stay internally consistent.
    """
    if axis is None and 11 in bands:
        axis = destripe_scan_axis(bands[11], n_detectors=n_detectors)
    elif axis is None:
        axis = 0
    out = {b: a.copy() for b, a in bands.items()}
    stats: dict[str, Any] = {"axis": axis, "n_detectors": n_detectors, "bands": {}}
    for b in band_ids:
        if b not in out:
            continue
        before = destripe_amplitude(out[b], n_detectors=n_detectors, axis=axis)
        out[b] = destripe_moment_match(out[b], n_detectors=n_detectors, axis=axis)
        after = destripe_amplitude(out[b], n_detectors=n_detectors, axis=axis)
        stats["bands"][str(b)] = {
            "amplitude_before": round(before, 6),
            "amplitude_after": round(after, 6),
        }
    return out, stats


def ast05_qa_good(qa: np.ndarray, max_adjacency: int = 1) -> np.ndarray:
    """Pixels usable from the AST_05 ``QA_DataPlane``.

    Bits 7–4 quality (0000 good; 0001–0111 suspect; 1xxx bad), bits 3–2 cloud
    (00 clear, 01 thin, 10 thick), bits 1–0 cloud adjacency (00 far … 11 very
    near) per the ASTER Level-2 product specification. Keeps good, clear pixels
    no closer to cloud than ``max_adjacency``.
    """
    qa = qa.astype(np.uint8)
    quality = qa >> 4
    cloud = (qa >> 2) & 0b11
    adjacency = qa & 0b11
    return (quality == 0) & (cloud == 0) & (adjacency <= max_adjacency)


def acquisition_token(granule_id: str) -> str | None:
    """Acquisition datetime (MMDDYYYYhhmmss) shared by AST_L1T / AST_05 / AST_08 of one scene."""
    m = _ACQUISITION_RE.search(granule_id or "")
    return m.group(1)[3:] if m else None


def load_ast05_granule(
    directory: os.PathLike[str] | str,
    granule_id: str,
) -> tuple[dict[int, np.ndarray], np.ndarray, rasterio.Affine, rasterio.crs.CRS, Profile]:
    """Load AST_05 emissivity B10–B14 (float, NaN outside QA-good valid pixels).

    Returns (bands, qa_good_mask, transform, crs, profile). Pixels failing the QA
    rule or with any band outside ``EMISSIVITY_VALID_RANGE`` are NaN in every band.
    """
    directory = Path(directory)
    bands: dict[int, np.ndarray] = {}
    transform = crs = profile = None
    for b in REQUIRED_TIR_BANDS:
        with rasterio.open(directory / f"{granule_id}_SRE_TIR_B{b}.tif") as src:
            raw = src.read(1)
            if transform is None:
                transform, crs, profile = src.transform, src.crs, src.profile.copy()
        arr = raw.astype(np.float64) * EMISSIVITY_SCALE
        arr[raw == 0] = np.nan
        bands[b] = arr
    with rasterio.open(directory / f"{granule_id}_SRE_TIR_QA_DataPlane.tif") as src:
        qa_good = ast05_qa_good(src.read(1))
    lo, hi = EMISSIVITY_VALID_RANGE
    valid = qa_good.copy()
    for arr in bands.values():
        valid &= np.isfinite(arr) & (arr >= lo) & (arr <= hi)
    for arr in bands.values():
        arr[~valid] = np.nan
    return bands, qa_good, transform, crs, profile


def compute_valid_data_footprint(
    band: np.ndarray,
    transform: rasterio.Affine,
    crs: rasterio.crs.CRS,
) -> "gpd.GeoDataFrame | None":
    """Return a GeoDataFrame containing the polygon of valid (non-NaN) pixels.

    The polygon is in the raster's native CRS, simplified to ~5 pixel widths to
    remove pixel-level jaggedness while preserving scene-boundary shape.  Returns
    None if no valid pixels exist or if the computation fails.
    """
    try:
        import geopandas as gpd
        from rasterio.features import shapes
        from shapely.geometry import shape
        from shapely.ops import unary_union

        valid_mask = np.isfinite(band).astype(np.uint8)
        if not valid_mask.any():
            return None

        geoms = [
            shape(geom)
            for geom, val in shapes(valid_mask, transform=transform)
            if val == 1
        ]
        if not geoms:
            return None

        footprint = unary_union(geoms)
        # Simplify to eliminate pixel-level staircasing; preserve scene boundary shape.
        _tol = 0.005 if crs.is_geographic else 450.0  # degrees or ~5 × 90 m pixels
        footprint = footprint.simplify(_tol)

        return gpd.GeoDataFrame(geometry=[footprint], crs=crs)
    except Exception:
        return None
