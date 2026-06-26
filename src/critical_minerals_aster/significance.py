"""
Statistical significance tests for ASTER anomaly-zone hit rates.

Null hypothesis (both tests): MRDS deposits are spatially uniform within the
site bounding box.  Under H₀ each deposit has probability p = (zone area /
bbox area) of landing inside a strong-anomaly zone.

Functions
---------
coverage_fraction   — zone area / bbox area in the zones' projected CRS
run_binomial        — exact one-sided p-value via scipy.stats.binomtest
run_permutation     — Monte Carlo p-value via rasterised grid sampling
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.stats import binomtest
from shapely.geometry import box
from shapely.ops import unary_union


def _bbox_in_crs(bbox_wgs84: tuple, crs) -> tuple[float, float, float, float]:
    """Return (minx, miny, maxx, maxy) for bbox_wgs84 reprojected to crs."""
    minlon, minlat, maxlon, maxlat = bbox_wgs84
    gdf = gpd.GeoDataFrame(geometry=[box(minlon, minlat, maxlon, maxlat)], crs="EPSG:4326")
    return gdf.to_crs(crs).geometry[0].bounds


def coverage_fraction(zones_gdf: gpd.GeoDataFrame, bbox_wgs84: tuple) -> float:
    """Fraction of projected bbox area covered by strong-anomaly zones (0–1)."""
    minx, miny, maxx, maxy = _bbox_in_crs(bbox_wgs84, zones_gdf.crs)
    bbox_area = (maxx - minx) * (maxy - miny)
    if bbox_area == 0:
        return 0.0
    union = unary_union(zones_gdf.geometry)
    return min(union.area / bbox_area, 1.0)


def run_binomial(n_hits: int, n_deposits: int, p: float) -> tuple[float, float]:
    """
    One-sided exact binomial test: H₀ hit_rate ≤ p.

    Returns
    -------
    p_value : float
    expected_hits : float  — n_deposits × p
    """
    if n_deposits == 0 or p <= 0:
        return 1.0, 0.0
    result = binomtest(n_hits, n_deposits, p, alternative="greater")
    return float(result.pvalue), n_deposits * p


def run_permutation(
    zones_gdf: gpd.GeoDataFrame,
    bbox_wgs84: tuple,
    n_deposits: int,
    n_hits: int,
    n_iter: int = 10_000,
    seed: int = 42,
    grid_res: int = 1_000,
) -> float:
    """
    Monte Carlo spatial permutation p-value.

    Rasterises the zone union onto a ``grid_res × grid_res`` grid covering
    the bbox, then for each iteration samples ``n_deposits`` grid cells
    uniformly at random and counts how many are inside a zone.

    Returns P(random_hits ≥ n_hits) over ``n_iter`` iterations.

    Mathematically equivalent to placing n_deposits random points uniformly
    in the bbox and checking containment, but runs in milliseconds via numpy
    array indexing rather than repeated point-in-polygon tests.
    """
    if n_deposits == 0:
        return 1.0

    minx, miny, maxx, maxy = _bbox_in_crs(bbox_wgs84, zones_gdf.crs)

    import rasterio.features
    from rasterio.transform import from_bounds as rio_from_bounds

    union = unary_union(zones_gdf.geometry)
    transform = rio_from_bounds(minx, miny, maxx, maxy, grid_res, grid_res)
    mask = rasterio.features.rasterize(
        [(union, 1)],
        out_shape=(grid_res, grid_res),
        transform=transform,
        fill=0,
        dtype=np.uint8,
    ).ravel()

    n_cells = len(mask)
    rng = np.random.default_rng(seed)
    # Sample all iterations at once: shape (n_iter, n_deposits)
    indices = rng.integers(0, n_cells, size=(n_iter, n_deposits))
    null_hits = mask[indices].sum(axis=1)
    return float((null_hits >= n_hits).mean())


# ---------------------------------------------------------------------------
# Wilson score confidence interval
# ---------------------------------------------------------------------------


def wilson_ci(k: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Wilson score confidence interval for a proportion k/n.

    Returns (lower, upper) as fractions in [0, 1].

    Why Wilson over the naive Wald interval (p ± z·se)?
    The Wald interval can produce negative lower bounds or bounds > 1 when
    the observed rate is extreme or n is small — both common in this dataset.
    Wilson is algebraically constrained to [0, 1] and has substantially
    better coverage properties for proportions near 0 or 1.

    Reference: Wilson (1927), J. Amer. Statist. Assoc.
    """
    if n == 0:
        return 0.0, 1.0
    import math
    # z-score for two-sided interval (e.g. 1.96 for 95%)
    z = {0.90: 1.645, 0.95: 1.960, 0.99: 2.576}.get(confidence, 1.960)
    p_hat = k / n
    denom = 1.0 + z ** 2 / n
    center = (p_hat + z ** 2 / (2 * n)) / denom
    margin = z * math.sqrt(p_hat * (1 - p_hat) / n + z ** 2 / (4 * n ** 2)) / denom
    return max(0.0, center - margin), min(1.0, center + margin)


# ---------------------------------------------------------------------------
# Footprint-aware coverage fraction (exposed for pipeline use)
# ---------------------------------------------------------------------------


def coverage_fraction_footprint(
    zones: gpd.GeoDataFrame,
    footprint: gpd.GeoDataFrame | None,
    bbox: tuple,
) -> float:
    """Zone area / denominator area, where denominator is the TIR footprint polygon
    (if available) or the rectangular bbox.

    Preferred over ``coverage_fraction`` when the TIR scene has a diagonal
    boundary (standard ASTER swath geometry) — the footprint polygon is more
    accurate than the bounding rectangle.
    """
    union = unary_union(zones.geometry)
    if footprint is not None and len(footprint):
        fp = footprint if footprint.crs == zones.crs else footprint.to_crs(zones.crs)
        denom = float(fp.iloc[0].geometry.area)
    else:
        minx, miny, maxx, maxy = _bbox_in_crs(bbox, zones.crs)
        denom = (maxx - minx) * (maxy - miny)
    if denom == 0:
        return 0.0
    return min(float(union.area / denom), 1.0)


# ---------------------------------------------------------------------------
# Annotate a site-summary DataFrame with uncertainty columns
# ---------------------------------------------------------------------------


def add_uncertainty_columns(
    summary: pd.DataFrame,
    zones: gpd.GeoDataFrame,
    footprint: gpd.GeoDataFrame | None,
    bbox: tuple,
    confidence: float = 0.95,
) -> pd.DataFrame:
    """Add Wilson CIs, null hit rate, and binomial p-value to every summary row.

    Works on the DataFrame produced by ``metrics.compute_site_summary()``.
    All per-row calculations use the same null probability (zone coverage
    fraction of the site footprint), because the null model is: a deposit
    is placed uniformly at random within the TIR footprint, so its probability
    of landing in a zone equals the fraction of the footprint covered by zones.

    New columns added
    -----------------
    null_hit_rate_pct : float
        Expected hit rate under H₀ = p_cover × 100.
    hit_rate_ci_low : float
        Wilson lower bound on observed hit rate (%), at *confidence* level.
    hit_rate_ci_high : float
        Wilson upper bound on observed hit rate (%), at *confidence* level.
    p_binomial : float
        One-sided exact binomial p-value: P(hits ≥ observed | n, p_cover).
        1.0 when n_deposits_bbox == 0 or p_cover == 0.
    """
    out = summary.copy()

    p_cover = coverage_fraction_footprint(zones, footprint, bbox)
    null_pct = round(p_cover * 100, 2)

    ci_lows, ci_highs, p_vals = [], [], []

    for _, row in out.iterrows():
        n = int(row["n_deposits_bbox"])
        k = int(row["n_deposits_in_zones"])
        lo, hi = wilson_ci(k, n, confidence)
        ci_lows.append(round(lo * 100, 2))
        ci_highs.append(round(hi * 100, 2))
        p_val, _ = run_binomial(k, n, p_cover)
        p_vals.append(round(p_val, 4))

    out["null_hit_rate_pct"] = null_pct
    out["hit_rate_ci_low"] = ci_lows
    out["hit_rate_ci_high"] = ci_highs
    out["p_binomial"] = p_vals
    return out
