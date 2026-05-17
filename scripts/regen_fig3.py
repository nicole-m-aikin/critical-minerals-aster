"""Regenerate figure 03 (deposit overlay with hillshade) for all sites.

Uses only cached data — zones GeoJSON, MRDS CSV, and reprojected DEM — so the
full ASTER classification pipeline does not need to re-run.
"""

from __future__ import annotations

import sys
from pathlib import Path

import rasterio

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import pandas as pd

from critical_minerals_aster.config import load_site_config
from critical_minerals_aster.metrics import filter_mrds_bbox, read_mrds_national, simplify_commodity
from critical_minerals_aster.mrds import mrds_to_points_gdf, spatial_join_deposits_zones
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.pipeline import save_deposit_overlay_figure
from critical_minerals_aster.structure import annotate_deposits_with_structure, load_structure_layers
from critical_minerals_aster.terrain import compute_hillshade_for_site

sites_dir = REPO_ROOT / "sites"
site_yamls = sorted(p for p in sites_dir.glob("*.yaml") if p.stem != "index")

for yaml_file in site_yamls:
    site = load_site_config(yaml_file)
    paths = site_paths_for(site, REPO_ROOT)

    print(f"\n=== {site.id} ===")

    if not paths.strong_zones_geojson.exists():
        print("  zones GeoJSON missing — skipping")
        continue

    import geopandas as gpd
    zones = gpd.read_file(paths.strong_zones_geojson)
    print(f"  zones: {len(zones)}")

    # Read raster metadata from reprojected DEM (same grid as TIR raster).
    dem_path = REPO_ROOT / "data" / "dem" / site.id / "dem_reprojected.tif"
    if not dem_path.exists():
        print("  reprojected DEM missing — figure will have no hillshade")
        raster_transform = None
        raster_shape = None
        raster_crs = None
    else:
        with rasterio.open(dem_path) as ds:
            raster_transform = ds.transform
            raster_shape = (ds.height, ds.width)
            raster_crs = ds.crs

    # Compute hillshade.
    hillshade = None
    if raster_transform is not None:
        try:
            hillshade = compute_hillshade_for_site(
                site, paths, raster_transform, raster_shape, raster_crs
            )
            print(f"  hillshade: {'ok' if hillshade is not None else 'failed'}")
        except Exception as exc:
            print(f"  hillshade error: {exc}")

    # Load MRDS deposits and spatial-join to zones.
    deposits = gpd.GeoDataFrame()
    if paths.mrds_csv.exists():
        try:
            mrds = read_mrds_national(paths)
            from critical_minerals_aster.spectral import raster_bbox_wgs84
            if raster_transform is not None:
                bbox = raster_bbox_wgs84(raster_transform, raster_shape, raster_crs)
            else:
                bbox = site.bbox_wgs84
            local = filter_mrds_bbox(mrds, bbox)
            deposits = mrds_to_points_gdf(local, zones.crs)
            joined, hits, _ = spatial_join_deposits_zones(deposits, zones)
            hit_ids = joined[joined["index_right"].notna()].index.unique()
            deposits["inside_zone"] = deposits.index.isin(hit_ids)
            deposits["commodity_group"] = deposits["commod1"].apply(simplify_commodity)
            print(f"  deposits: {len(deposits)} ({len(hit_ids)} inside zones)")
        except Exception as exc:
            print(f"  deposits error: {exc}")

    # Load structure layers and compute on-structure deposit count.
    structs = None
    n_on_structure: int | None = None
    if site.structure_layers and len(deposits) > 0:
        try:
            target_crs = zones.crs if len(zones) else deposits.crs
            structs = load_structure_layers(site, REPO_ROOT, target_crs)
            annotated = annotate_deposits_with_structure(deposits, site, paths, structs=structs)
            n_on_structure = int(annotated["on_structure"].sum())
            print(f"  on_structure: {n_on_structure}/{len(deposits)}")
        except Exception as exc:
            print(f"  structure error: {exc}")

    save_deposit_overlay_figure(
        site, paths, zones, deposits, REPO_ROOT,
        hillshade=hillshade,
        raster_transform=raster_transform,
        raster_shape=raster_shape,
        structs=structs,
        n_on_structure=n_on_structure,
        n_total_deposits=len(deposits) if len(deposits) > 0 else None,
    )
    print(f"  saved → {paths.figures_dir / '03_deposit_overlay.png'}")

print("\nDone.")
