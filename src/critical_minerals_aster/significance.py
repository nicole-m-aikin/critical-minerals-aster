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
