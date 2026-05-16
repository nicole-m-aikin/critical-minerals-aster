"""Percentile classification and vectorization of anomaly zones."""

from __future__ import annotations

import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import shapes
from shapely.geometry import shape


def classify_percentiles(
    ratio: np.ndarray,
    low_pct: float = 70.0,
    high_pct: float = 90.0,
) -> tuple[np.ndarray, float, float]:
    """3-class map: 0 background, 1 moderate, 2 strong; NaN stays 0."""
    low = float(np.nanpercentile(ratio, low_pct))
    high = float(np.nanpercentile(ratio, high_pct))
    classes = np.zeros_like(ratio, dtype=np.uint8)
    classes[ratio >= low] = 1
    classes[ratio >= high] = 2
    classes[np.isnan(ratio)] = 0
    return classes, low, high


def combined_score(
    silica_cls: np.ndarray,
    carbonate_cls: np.ndarray,
    mafic_cls: np.ndarray,
) -> np.ndarray:
    """Sum of three 0–2 class maps → combined score 0–6."""
    return (
        silica_cls.astype(np.uint8)
        + carbonate_cls.astype(np.uint8)
        + mafic_cls.astype(np.uint8)
    )


def vectorize_strong_zones(
    combined: np.ndarray,
    transform: rasterio.Affine,
    crs: rasterio.crs.CRS,
    min_score: int = 3,
) -> gpd.GeoDataFrame:
    """Polygonize pixels where combined score >= min_score (default: strong anomalies)."""
    mask = (combined >= min_score).astype(np.uint8)
    geoms = []
    for geom, val in shapes(mask, transform=transform):
        if val == 1:
            geoms.append(shape(geom))
    zones = gpd.GeoDataFrame(geometry=geoms, crs=crs)
    zones["area_km2"] = zones.geometry.area / 1e6
    return zones
