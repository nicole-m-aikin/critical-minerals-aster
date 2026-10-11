"""Box-independence sensitivity helpers: bbox variants, reclassification, critical hit test.

Supports ``scripts/bbox_independence_sensitivity.py``, which asks whether a
site's critical-mineral enrichment depends on the analyst-drawn
``bbox_wgs84``. The helpers mirror production code rather than reimplementing
it:

- :func:`classify_strong_zones` reproduces the ratio-mosaic branch of
  ``pipeline.run_classification`` (percentiles -> combined score -> strong
  zones, valid footprint from B10) on whatever extent it is given.
- :func:`critical_hit_test` reproduces ``metrics.compute_site_summary`` +
  ``significance.add_uncertainty_columns`` restricted to critical minerals,
  i.e. the quantity tabulated by ``scripts/site_specific_null_significance.py``.
"""

from __future__ import annotations

from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box

from critical_minerals_aster.classification import (
    classify_percentiles,
    combined_score,
    vectorize_strong_zones,
)
from critical_minerals_aster.config import BBox, ClassificationParams
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    is_critical_mineral,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.significance import coverage_fraction_footprint, run_binomial
from critical_minerals_aster.spectral import compute_valid_data_footprint, raster_bbox_wgs84

# Expansion-template box side (docs/expansion_plan_100_sites.md: "roughly 0.4-0.5
# degrees on a side centered on the district").
TEMPLATE_BOX_DEG = 0.46


def perturb_bbox(
    bbox: BBox,
    dx_frac: float = 0.0,
    dy_frac: float = 0.0,
    scale: float = 1.0,
) -> BBox:
    """Shift the box centre by fractions of its width/height and rescale its sides.

    ``scale`` multiplies each side length about the (shifted) centre, so
    ``scale=0.85`` trims 7.5% from every edge.
    """
    lon0, lat0, lon1, lat1 = bbox
    w, h = lon1 - lon0, lat1 - lat0
    cx = (lon0 + lon1) / 2 + dx_frac * w
    cy = (lat0 + lat1) / 2 + dy_frac * h
    hw, hh = w * scale / 2, h * scale / 2
    return (cx - hw, cy - hh, cx + hw, cy + hh)


def perturbation_variants(bbox: BBox, frac: float = 0.15) -> dict[str, BBox]:
    """Eight shifted / rescaled copies of ``bbox`` (six shifts, one shrink, one expand)."""
    pct = int(round(frac * 100))
    specs: dict[str, tuple[float, float, float]] = {
        f"shift_n_{pct}": (0.0, frac, 1.0),
        f"shift_s_{pct}": (0.0, -frac, 1.0),
        f"shift_e_{pct}": (frac, 0.0, 1.0),
        f"shift_w_{pct}": (-frac, 0.0, 1.0),
        f"shift_ne_{pct}": (frac, frac, 1.0),
        f"shift_sw_{pct}": (-frac, -frac, 1.0),
        f"shrink_{pct}": (0.0, 0.0, 1.0 - frac),
        f"expand_{pct}": (0.0, 0.0, 1.0 + frac),
    }
    return {name: perturb_bbox(bbox, dx, dy, s) for name, (dx, dy, s) in specs.items()}


def centred_box(lon: float, lat: float, size_deg: float = TEMPLATE_BOX_DEG) -> BBox:
    """Square (in degrees) box of side ``size_deg`` centred on (lon, lat)."""
    half = size_deg / 2
    return (lon - half, lat - half, lon + half, lat + half)


def bbox_overlap_fraction(a: BBox, b: BBox) -> float:
    """Fraction of box ``a``'s (degree) area that lies inside box ``b``."""
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    area_a = (ax1 - ax0) * (ay1 - ay0)
    if area_a <= 0:
        return 0.0
    ix = max(0.0, min(ax1, bx1) - max(ax0, bx0))
    iy = max(0.0, min(ay1, by1) - max(ay0, by0))
    return ix * iy / area_a


def raster_extent_polygon(transform: Any, shape: tuple[int, int]):
    """Rectangle (raster CRS) covered by a north-up grid with this transform and shape."""
    rows, cols = shape
    x0, y1 = transform.c, transform.f
    x1 = x0 + cols * transform.a
    y0 = y1 + rows * transform.e
    return box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))


def footprint_centroid_wgs84(footprint: gpd.GeoDataFrame) -> tuple[float, float]:
    """(lon, lat) of the footprint polygon centroid, computed in the footprint's own CRS."""
    c = footprint.geometry.union_all().centroid
    pt = gpd.GeoSeries([c], crs=footprint.crs).to_crs("EPSG:4326").iloc[0]
    return float(pt.x), float(pt.y)


def classify_strong_zones(
    silica: np.ndarray,
    carbonate: np.ndarray,
    mafic: np.ndarray,
    b10: np.ndarray,
    transform: Any,
    crs: Any,
    cp: ClassificationParams,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame | None, BBox]:
    """Percentile-classify the three ratios over the given extent and vectorize strong zones.

    Returns ``(zones, footprint, raster_bbox_wgs84)``. Percentile thresholds are
    computed over exactly the pixels passed in, so the extent *is* the
    experimental variable.
    """
    silica_cls, _, _ = classify_percentiles(silica, cp.low_pct, cp.high_pct)
    carbonate_cls, _, _ = classify_percentiles(carbonate, cp.low_pct, cp.high_pct)
    mafic_cls, _, _ = classify_percentiles(mafic, cp.low_pct, cp.high_pct)
    combined = combined_score(silica_cls, carbonate_cls, mafic_cls)
    zones = vectorize_strong_zones(combined, transform, crs, min_score=cp.strong_score_min)
    footprint = compute_valid_data_footprint(b10, transform, crs)
    return zones, footprint, raster_bbox_wgs84(transform, silica.shape, crs)


def restrict_to_region(
    zones: gpd.GeoDataFrame,
    footprint: gpd.GeoDataFrame,
    region: Any,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Clip zones and footprint to ``region`` (a shapely geometry in the zones' CRS).

    Used for the "outside the analyst box" variant (region = scene footprint
    minus the analyst raster window) and for no-reclassification box
    perturbation of cached zones.
    """
    fp_geom = footprint.geometry.union_all().intersection(region)
    fp_out = gpd.GeoDataFrame(geometry=[fp_geom], crs=footprint.crs)
    if len(zones):
        z = gpd.clip(zones, gpd.GeoDataFrame(geometry=[fp_geom], crs=zones.crs))
        z = z[~z.geometry.is_empty].copy()
    else:
        z = zones.copy()
    z["area_km2"] = z.geometry.area / 1e6
    return z, fp_out


def summarize_hits(n: int, hits: int, p0: float) -> dict[str, float | int]:
    """Hit rate, expected hits, enrichment and one-sided binomial p against null ``p0``."""
    p_value, expected = run_binomial(hits, n, p0)
    rate = hits / n if n else float("nan")
    enrichment = rate / p0 if (n and p0 > 0) else float("nan")
    return {
        "n_crit": int(n),
        "hits_crit": int(hits),
        "hit_rate_crit_pct": round(rate * 100, 2) if n else float("nan"),
        "p0_null": round(p0, 4),
        "expected_hits_null": round(expected, 1),
        "enrichment": round(enrichment, 2) if np.isfinite(enrichment) else float("nan"),
        "p_binomial": p_value,
    }


def critical_hit_test(
    mrds: pd.DataFrame,
    zones: gpd.GeoDataFrame,
    footprint: gpd.GeoDataFrame | None,
    raster_bbox: BBox,
) -> dict[str, float | int]:
    """Critical-mineral hits vs the site-specific null p0 = zone area / footprint area.

    Same deposit set and null as the published per-site table: MRDS rows in
    ``raster_bbox`` (WGS84 pre-filter), clipped to the valid footprint,
    ``within`` zone polygons, critical = ``is_critical_mineral(earth_mri_category)``.
    """
    local = filter_mrds_bbox(mrds, raster_bbox)
    deposits = mrds_to_points_gdf(local, zones.crs)
    if footprint is not None and len(footprint):
        fp = footprint if footprint.crs == zones.crs else footprint.to_crs(zones.crs)
        deposits = gpd.clip(deposits, fp)
    if len(zones) and len(deposits):
        joined, _, _ = spatial_join_deposits_zones(deposits, zones)
        hit_ids = set(joined.index[joined["index_right"].notna()])
    else:
        hit_ids = set()
    deposits = reclassify_mrds_earth_mri(deposits)
    crit = deposits[deposits["earth_mri_category"].map(is_critical_mineral).astype(bool)]
    hits = int(crit.index.isin(list(hit_ids)).sum())
    p0 = coverage_fraction_footprint(zones, footprint, raster_bbox) if len(zones) else 0.0

    out = summarize_hits(len(crit), hits, p0)
    fp_area = (
        float(footprint.geometry.area.sum()) / 1e6
        if footprint is not None and len(footprint) else float("nan")
    )
    out["footprint_km2"] = round(fp_area, 1)
    out["zone_km2"] = round(float(zones.geometry.area.sum()) / 1e6, 2) if len(zones) else 0.0
    out["n_zones"] = int(len(zones))
    return out
