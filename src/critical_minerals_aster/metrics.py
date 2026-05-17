"""Per-site summary metrics and commodity grouping."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd

from critical_minerals_aster.config import BBox, SiteConfig
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    is_critical_mineral,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    reclassify_mrds_mineral_system,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.paths import SitePaths


def simplify_commodity(commod: Any) -> str:
    c = str(commod).lower()
    if "lithium" in c:
        return "Lithium"
    if "mercury" in c:
        return "Mercury"
    if "uranium" in c:
        return "Uranium"
    if "gold" in c or "silver" in c:
        return "Gold/Silver"
    if "antimony" in c:
        return "Antimony"
    if "sand" in c or "gravel" in c:
        return "Sand/Gravel"
    if "stone" in c:
        return "Stone"
    if "gemstone" in c or "semiprecious" in c:
        return "Gemstone"
    return "Other"


def read_mrds_national(paths: SitePaths) -> pd.DataFrame:
    csv_path = paths.mrds_csv
    if not csv_path.is_file():
        raise FileNotFoundError(
            f"MRDS CSV not found at {csv_path}. Download via notebook 04 or place under data/mrds/."
        )
    return pd.read_csv(csv_path, low_memory=False)


def compute_site_summary(
    site: SiteConfig,
    paths: SitePaths,
    zones: gpd.GeoDataFrame,
    granule_id: str,
    mrds_bbox: BBox | None = None,
    n_on_structure: int | None = None,
    mean_nearest_m: float | None = None,
    annotated_deposits: gpd.GeoDataFrame | None = None,
) -> pd.DataFrame:
    """Site-level summary row plus one row per commodity group.

    Parameters
    ----------
    mrds_bbox:
        WGS84 bounding box used to filter MRDS deposits.  Defaults to
        ``site.bbox_wgs84`` when not supplied, but callers should pass the
        actual raster extent (from ``run_classification``) so only deposits
        within the TIR coverage area are counted.
    n_on_structure:
        Count of MRDS deposits that lie within the configured buffer of any
        structure layer.  None when no structure layers are configured.
    mean_nearest_m:
        Mean nearest-structure distance (metres) across all deposits in the
        site bbox.  None when no structure layers are configured.
    annotated_deposits:
        GeoDataFrame already annotated with ``nearest_structure_m`` and
        ``on_structure`` columns (output of
        :func:`~critical_minerals_aster.structure.annotate_deposits_with_structure`).
        When supplied, per-group structure metrics are computed for
        commodity / earth_mri / mineral_system sub-rows.
    """
    mrds = read_mrds_national(paths)
    effective_bbox: BBox = mrds_bbox if mrds_bbox is not None else site.bbox_wgs84
    local = filter_mrds_bbox(mrds, effective_bbox)
    deposits = mrds_to_points_gdf(local, zones.crs)
    joined, _, _ = spatial_join_deposits_zones(deposits, zones)

    inside = joined[joined["index_right"].notna()]
    n_dep = len(deposits)
    n_hit = inside.index.nunique()
    hit_rate = (n_hit / n_dep * 100.0) if n_dep else 0.0

    total_area = float(zones["area_km2"].sum()) if len(zones) else 0.0
    median_zone = float(zones["area_km2"].median()) if len(zones) else 0.0

    base = {
        "site_id": site.id,
        "site_name": site.name,
        "granule_id": granule_id,
        "n_zones": len(zones),
        "total_anomaly_km2": round(total_area, 2),
        "median_zone_km2": round(median_zone, 2),
        "n_deposits_bbox": n_dep,
        "n_deposits_in_zones": int(n_hit),
        "hit_rate_pct": round(hit_rate, 1),
        "layout": site.layout,
    }

    rows: list[dict] = [
        {
            **base,
            "row_type": "site",
            "commodity_group": "",
            "earth_mri_category": "",
            "is_critical_mineral": False,
            "mineral_system": "",
            "n_deposits_on_structure": n_on_structure,
            "mean_nearest_structure_m": mean_nearest_m,
        }
    ]

    if n_dep and "commod1" in deposits.columns:
        hit_ids = set(inside.index.unique())
        deposits = deposits.copy()
        deposits["inside_zone"] = deposits.index.isin(hit_ids)
        deposits["commodity_group"] = deposits["commod1"].apply(simplify_commodity)
        deposits = reclassify_mrds_earth_mri(deposits)
        deposits = reclassify_mrds_mineral_system(deposits)

        # Merge structure annotation columns if provided so per-group metrics
        # are computed from the same deposit set used for site-level metrics.
        _has_struct = False
        if annotated_deposits is not None and not annotated_deposits.empty:
            if "nearest_structure_m" in annotated_deposits.columns:
                deposits = deposits.join(
                    annotated_deposits[["nearest_structure_m", "on_structure"]],
                    how="left",
                )
                _has_struct = True

        def _struct_metrics(
            grp_df: pd.DataFrame,
        ) -> tuple[int | None, float | None]:
            """Return (n_on_structure, mean_nearest_m) for a deposit sub-group."""
            if not _has_struct or "on_structure" not in grp_df.columns:
                return None, None
            n_on = int(grp_df["on_structure"].sum())
            valid = grp_df["nearest_structure_m"].dropna()
            mean_m = float(valid.mean()) if not valid.empty else None
            return n_on, mean_m

        for grp, grp_df in deposits.groupby("commodity_group"):
            inside_n = int(grp_df["inside_zone"].sum())
            total = len(grp_df)
            n_on, mean_m = _struct_metrics(grp_df)
            rows.append(
                {
                    **base,
                    "row_type": "commodity",
                    "commodity_group": grp,
                    "earth_mri_category": "",
                    "is_critical_mineral": False,
                    "mineral_system": "",
                    "n_deposits_bbox": total,
                    "n_deposits_in_zones": inside_n,
                    "hit_rate_pct": round(inside_n / total * 100, 1) if total else 0.0,
                    "n_deposits_on_structure": n_on,
                    "mean_nearest_structure_m": mean_m,
                }
            )

        for grp, grp_df in deposits.groupby("earth_mri_category"):
            inside_n = int(grp_df["inside_zone"].sum())
            total = len(grp_df)
            n_on, mean_m = _struct_metrics(grp_df)
            rows.append(
                {
                    **base,
                    "row_type": "earth_mri",
                    "commodity_group": "",
                    "earth_mri_category": grp,
                    "is_critical_mineral": is_critical_mineral(grp),
                    "mineral_system": "",
                    "n_deposits_bbox": total,
                    "n_deposits_in_zones": inside_n,
                    "hit_rate_pct": round(inside_n / total * 100, 1) if total else 0.0,
                    "n_deposits_on_structure": n_on,
                    "mean_nearest_structure_m": mean_m,
                }
            )

        for grp, grp_df in deposits.groupby("mineral_system"):
            inside_n = int(grp_df["inside_zone"].sum())
            total = len(grp_df)
            n_on, mean_m = _struct_metrics(grp_df)
            rows.append(
                {
                    **base,
                    "row_type": "mineral_system",
                    "commodity_group": "",
                    "earth_mri_category": "",
                    "is_critical_mineral": False,
                    "mineral_system": grp,
                    "n_deposits_bbox": total,
                    "n_deposits_in_zones": inside_n,
                    "hit_rate_pct": round(inside_n / total * 100, 1) if total else 0.0,
                    "n_deposits_on_structure": n_on,
                    "mean_nearest_structure_m": mean_m,
                }
            )

    return pd.DataFrame(rows)


def write_site_summary(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
