"""Structural geology layer helpers (faults, contacts, folds)."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd

from critical_minerals_aster.config import SiteConfig, StructureLayer
from critical_minerals_aster.paths import SitePaths


def load_structure_layers(
    site: SiteConfig,
    repo_root: Path,
    target_crs,
) -> gpd.GeoDataFrame:
    """Load and concatenate configured structure layers, reprojected to target_crs."""
    if not site.structure_layers:
        return gpd.GeoDataFrame(geometry=[], crs=target_crs)

    frames: list[gpd.GeoDataFrame] = []
    for layer in site.structure_layers:
        path = Path(layer.path)
        if not path.is_absolute():
            path = repo_root / path
        gdf = gpd.read_file(path)
        if gdf.crs is None:
            gdf = gdf.set_crs("EPSG:4326")
        gdf = gdf.to_crs(target_crs)
        gdf["structure_type"] = layer.type
        gdf["buffer_m"] = layer.buffer_m
        frames.append(gdf)

    if not frames:
        return gpd.GeoDataFrame(geometry=[], crs=target_crs)
    return pd.concat(frames, ignore_index=True)


def nearest_structure_distance_m(
    points: gpd.GeoDataFrame,
    structures: gpd.GeoDataFrame,
) -> pd.Series:
    """Minimum distance (m) from each point to any structure geometry."""
    if structures.empty or points.empty:
        return pd.Series([float("nan")] * len(points), index=points.index)

    distances = []
    for pt in points.geometry:
        dists = structures.geometry.distance(pt)
        distances.append(float(dists.min()))
    return pd.Series(distances, index=points.index)


def points_on_structure(
    points: gpd.GeoDataFrame,
    structures: gpd.GeoDataFrame,
    buffer_m: float | None = None,
) -> pd.Series:
    """True if point lies within buffer of any structure feature."""
    if structures.empty or points.empty:
        return pd.Series(False, index=points.index)

    buf = buffer_m if buffer_m is not None else 500.0
    buffered = structures.copy()
    buffered["geometry"] = buffered.geometry.buffer(buf)
    joined = gpd.sjoin(points, buffered, how="left", predicate="within")
    on_idx = joined[joined["index_right"].notna()].index.unique()
    return pd.Series(points.index.isin(on_idx), index=points.index)


def annotate_deposits_with_structure(
    deposits: gpd.GeoDataFrame,
    site: SiteConfig,
    paths: SitePaths,
) -> gpd.GeoDataFrame:
    """Add nearest_structure_m and on_structure columns when layers are configured."""
    if not site.structure_layers:
        return deposits

    structures = load_structure_layers(site, paths.repo_root, deposits.crs)
    out = deposits.copy()
    out["nearest_structure_m"] = nearest_structure_distance_m(out, structures)
    default_buffer = site.structure_layers[0].buffer_m
    out["on_structure"] = points_on_structure(out, structures, default_buffer)
    return out
