"""MRDS point filtering and CRS-safe joins to anomaly polygons."""

from __future__ import annotations

from typing import Tuple

import geopandas as gpd
import pandas as pd

from critical_minerals_aster.config import BBox


def filter_mrds_bbox(mrds: pd.DataFrame, bbox: BBox) -> pd.DataFrame:
    """Subset MRDS rows whose longitude/latitude fall inside bbox (WGS84)."""
    lon0, lat0, lon1, lat1 = bbox
    return mrds.loc[
        (mrds["longitude"] >= lon0)
        & (mrds["longitude"] <= lon1)
        & (mrds["latitude"] >= lat0)
        & (mrds["latitude"] <= lat1)
    ].copy()


def mrds_to_points_gdf(mrds_local: pd.DataFrame, target_crs) -> gpd.GeoDataFrame:
    """Build point GeoDataFrame from MRDS lat/lon; reproject to target CRS (e.g. zones.crs)."""
    gdf = gpd.GeoDataFrame(
        mrds_local,
        geometry=gpd.points_from_xy(mrds_local.longitude, mrds_local.latitude),
        crs="EPSG:4326",
    )
    return gdf.to_crs(target_crs)


def spatial_join_deposits_zones(
    deposits: gpd.GeoDataFrame,
    zones: gpd.GeoDataFrame,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Join deposits to zones with predicate 'within'. Returns joined, hits, misses."""
    joined = gpd.sjoin(deposits, zones, how="left", predicate="within")
    hits = joined[joined["index_right"].notna()]
    misses = joined[joined["index_right"].isna()]
    return joined, hits, misses
