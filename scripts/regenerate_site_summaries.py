#!/usr/bin/env python
"""
Phase 2 — regenerate all 45 site summary CSVs with the site-specific null.

This is a *cache-only* re-run: it reuses already-computed zone polygons
(data/sites/{id}/vectors/strong_anomaly_zones.geojson), TIR footprints and
granule ids (results/{id}_provenance.json), MRDS records (data/mrds.csv),
and cached structure/fault layers (data/structures/*.geojson). It never
touches raw ASTER rasters and makes no network calls, so it can be re-run
offline any time the classification or MRDS logic changes.

It mirrors the tail of pipeline.run_site() (structure annotation ->
compute_site_summary -> add_uncertainty_columns -> write_site_summary)
exactly, so the output is byte-for-byte what a full pipeline re-run would
produce for the summary CSV, without needing to reprocess rasters.

Why this script exists (see docs/results.md 1.8): the
checked-in results/*_summary.csv files predate add_uncertainty_columns
(the site-specific null machinery) and results/significance_critical_only.csv
is an orphaned artifact of an even older pooled-null computation. Both need
a fresh, code-consistent regeneration before any paper number can be trusted.

Usage:
    conda run -n aster-minerals python scripts/regenerate_site_summaries.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import pandas as pd
import yaml

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import SitePaths, site_paths_for
from critical_minerals_aster.metrics import (
    compute_site_summary,
    write_site_summary,
    read_mrds_national,
)
from critical_minerals_aster.mrds import filter_mrds_bbox, mrds_to_points_gdf
from critical_minerals_aster.structure import (
    load_structure_layers,
    annotate_deposits_with_structure,
    nearest_structure_distance_m,
    points_on_structure,
)
from critical_minerals_aster.significance import add_uncertainty_columns

REPO_ROOT = Path(__file__).parent.parent


def load_site_ids() -> list[str]:
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        return yaml.safe_load(f)["sites"]


def load_cached_structs(site, target_crs) -> "gpd.GeoDataFrame | None":
    """Cache-only structure loader: configured layers, else a cached
    auto-fetch file if one already exists on disk. Never touches the network."""
    if site.structure_layers:
        return load_structure_layers(site, REPO_ROOT, target_crs)
    auto_path = REPO_ROOT / "data" / "structures" / f"{site.id}_faults_auto.geojson"
    if auto_path.is_file() and auto_path.stat().st_size > 100:
        gdf = gpd.read_file(auto_path)
        if not gdf.empty:
            return gdf.to_crs(target_crs)
    return None


def main() -> None:
    site_ids = load_site_ids()
    n_ok, n_skip = 0, 0

    for site_id in site_ids:
        print(f"{site_id:22s}", end=" ", flush=True)

        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        paths = site_paths_for(site, REPO_ROOT)

        zones_path = paths.strong_zones_geojson
        prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
        if not zones_path.exists() or not prov_path.exists():
            print("SKIP — missing zones or provenance")
            n_skip += 1
            continue

        zones = gpd.read_file(zones_path)
        if zones.empty:
            print("SKIP — empty zones")
            n_skip += 1
            continue

        with open(prov_path) as f:
            prov = json.load(f)
        granule_id = prov.get("granule_id", "unknown")
        raster_bbox = tuple(prov["raster_bbox_wgs84"]) if prov.get("raster_bbox_wgs84") else site.bbox_wgs84

        tir_footprint = None
        if prov.get("tir_footprint_wgs84") is not None:
            try:
                tir_footprint = gpd.GeoDataFrame.from_features(
                    prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
                ).to_crs(zones.crs)
            except Exception:
                tir_footprint = None

        try:
            mrds = read_mrds_national(paths)
        except FileNotFoundError as exc:
            print(f"SKIP — {exc}")
            n_skip += 1
            continue

        local = filter_mrds_bbox(mrds, raster_bbox)
        deposits = mrds_to_points_gdf(local, zones.crs)
        if tir_footprint is not None and len(tir_footprint):
            try:
                deposits = gpd.clip(deposits, tir_footprint)
            except Exception:
                pass

        structs = load_cached_structs(site, zones.crs)
        n_on_structure = None
        mean_nearest_m = None
        annotated = None
        if structs is not None and not structs.empty and len(deposits):
            if site.structure_layers:
                annotated = annotate_deposits_with_structure(deposits, site, paths, structs=structs)
            else:
                # No configured structure_layers -> this came from the cached
                # auto-fetch file; annotate manually with the same 500 m
                # default buffer pipeline.py uses for auto-fetched structures.
                annotated = deposits.copy()
                annotated["nearest_structure_m"] = nearest_structure_distance_m(annotated, structs)
                annotated["on_structure"] = points_on_structure(annotated, structs, 500.0)
            n_on_structure = int(annotated["on_structure"].sum())
            mean_nearest_m = float(annotated["nearest_structure_m"].mean())

        summary = compute_site_summary(
            site, paths, zones, granule_id,
            mrds_bbox=raster_bbox,
            tir_footprint=tir_footprint,
            n_on_structure=n_on_structure,
            mean_nearest_m=mean_nearest_m,
            annotated_deposits=annotated,
        )
        summary = add_uncertainty_columns(summary, zones, tir_footprint, raster_bbox)
        write_site_summary(summary, paths.site_summary_csv)

        site_row = summary[summary["row_type"] == "site"].iloc[0]
        print(
            f"n={int(site_row['n_deposits_bbox']):4d}  hits={int(site_row['n_deposits_in_zones']):3d}"
            f"  hit%={site_row['hit_rate_pct']:5.1f}  null%={site_row['null_hit_rate_pct']:5.2f}"
            f"  p={site_row['p_binomial']:.4g}"
        )
        n_ok += 1

    print(f"\n{n_ok} sites regenerated, {n_skip} skipped.")


if __name__ == "__main__":
    main()
