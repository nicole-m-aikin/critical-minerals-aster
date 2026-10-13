"""Pixel labels: near critical MRDS deposits + stratum-matched hard negatives."""

from __future__ import annotations

from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import xy
from shapely.geometry import Point, box

from critical_minerals_aster.config import SiteConfig
from critical_minerals_aster.lithology import (
    STRATUM_CODES,
    assign_cover_stratum,
    cover_strata,
    rasterize_cover_strata,
)
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    is_critical_mineral,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
)
from critical_minerals_aster.paths import SitePaths

DEFAULT_POS_RADIUS_M = 180.0  # 2 × 90 m
ERA_CUTOFF = 1950


def _best_year(row: pd.Series) -> float:
    dy = row.get("disc_yr")
    if pd.notna(dy):
        try:
            return float(dy)
        except (TypeError, ValueError):
            pass
    yfp = row.get("yr_fst_prd")
    if pd.notna(yfp):
        try:
            return float(yfp)
        except (TypeError, ValueError):
            pass
    return float("nan")


def load_critical_deposits(
    site: SiteConfig,
    paths: SitePaths,
    *,
    target_crs,
    bbox_wgs84: tuple[float, float, float, float] | None = None,
) -> gpd.GeoDataFrame:
    """Critical-mineral MRDS points in site bbox, with ``disc_yr`` best-year."""
    mrds = read_mrds_national(paths)
    bbox = bbox_wgs84 or site.bbox_wgs84
    local = filter_mrds_bbox(mrds, bbox)
    deps = mrds_to_points_gdf(local, target_crs)
    if deps.empty:
        return deps
    deps = reclassify_mrds_earth_mri(deps)
    deps = deps[deps["earth_mri_category"].map(is_critical_mineral)].copy()
    if deps.empty:
        return deps
    deps["disc_yr"] = deps.apply(_best_year, axis=1)
    return deps


def raster_cell_centers(shape: tuple[int, int], transform) -> tuple[np.ndarray, np.ndarray]:
    rows, cols = np.indices(shape)
    xs, ys = xy(transform, rows, cols, offset="center")
    return np.asarray(xs), np.asarray(ys)


def label_grid(
    *,
    shape: tuple[int, int],
    transform,
    crs,
    deposits: gpd.GeoDataFrame,
    stratum_grid: np.ndarray,
    pos_radius_m: float = DEFAULT_POS_RADIUS_M,
    n_neg_per_pos: int = 5,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (label, disc_yr, sample_mask) grids.

    Positives: cells whose center is within ``pos_radius_m`` of a critical deposit.
    Hard negatives: same stratum as at least one positive deposit, outside the
    positive buffer, subsampled to ``n_neg_per_pos`` × n_pos (capped).
    """
    rng = rng or np.random.default_rng(0)
    h, w = shape
    label = np.zeros((h, w), dtype=np.int8)
    year = np.full((h, w), np.nan, dtype=np.float32)
    sample = np.zeros((h, w), dtype=bool)

    if deposits is None or deposits.empty:
        return label, year, sample

    xs, ys = raster_cell_centers(shape, transform)
    dep_x = deposits.geometry.x.to_numpy()
    dep_y = deposits.geometry.y.to_numpy()
    dep_yr = (
        deposits["disc_yr"].to_numpy(dtype=float)
        if "disc_yr" in deposits.columns
        else np.full(len(deposits), np.nan)
    )

    from scipy.spatial import cKDTree

    tree = cKDTree(np.column_stack([dep_x, dep_y]))
    flat_xy = np.column_stack([xs.ravel(), ys.ravel()])
    min_d, nn = tree.query(flat_xy, k=1, workers=-1)
    nearest_yr = dep_yr[nn].astype(np.float32)
    min_d = min_d.reshape(shape)
    nearest_yr = nearest_yr.reshape(shape)
    pos = min_d <= pos_radius_m
    label[pos] = 1
    year[pos] = nearest_yr[pos]
    sample[pos] = True

    # Hard negatives: same strata as positives, outside buffer
    pos_strata = set(int(v) for v in np.unique(stratum_grid[pos]) if v > 0)
    if not pos_strata:
        pos_strata = {1, 2, 3}
    far = (~pos) & np.isin(stratum_grid, list(pos_strata)) & (stratum_grid > 0)
    n_pos = int(pos.sum())
    n_neg_target = min(int(far.sum()), max(n_pos * n_neg_per_pos, n_pos))
    far_idx = np.flatnonzero(far.ravel())
    if far_idx.size and n_neg_target > 0:
        pick = rng.choice(far_idx, size=min(n_neg_target, far_idx.size), replace=False)
        sample.ravel()[pick] = True

    return label, year, sample


def build_stratum_grid(
    site_id: str,
    *,
    shape: tuple[int, int],
    transform,
    crs,
    bbox_wgs84: tuple[float, float, float, float],
    repo_root,
) -> np.ndarray:
    """Rasterize bedrock/cover/unmapped onto the study grid."""
    # Footprint = full raster bounds
    h, w = shape
    minx = transform.c
    maxy = transform.f
    maxx = minx + w * transform.a
    miny = maxy + h * transform.e  # e is typically negative
    footprint = box(min(minx, maxx), min(miny, maxy), max(minx, maxx), max(miny, maxy))
    strata = cover_strata(
        site_id,
        footprint,
        crs,
        bbox_wgs84=bbox_wgs84,
        repo_root=repo_root,
    )
    return rasterize_cover_strata(strata, transform, shape)
