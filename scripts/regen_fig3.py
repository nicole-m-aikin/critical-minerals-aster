"""Regenerate figure 03 (dual-panel deposit overlay) for one or all sites.

Uses cached zones GeoJSON, MRDS, DEM/hillshade, and optional basemap tiles —
no full ASTER classification rerun required.
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def _regen_fig3_worker(site_id: str, repo_root_str: str) -> tuple[str, bool, str]:
    """Picklable worker: regenerate fig 03 for one site."""
    repo_root = Path(repo_root_str)
    sys.path.insert(0, str(repo_root / "src"))

    import geopandas as gpd
    import rasterio

    from critical_minerals_aster.basemap import fetch_satellite_basemap_for_site
    from critical_minerals_aster.config import load_site_by_id
    from critical_minerals_aster.metrics import filter_mrds_bbox, read_mrds_national, simplify_commodity
    from critical_minerals_aster.mrds import mrds_to_points_gdf, spatial_join_deposits_zones
    from critical_minerals_aster.paths import site_paths_for
    from critical_minerals_aster.pipeline import (
        auto_fetch_structure,
        load_tir_footprint_from_provenance,
        save_deposit_overlay_figure,
    )
    from critical_minerals_aster.spectral import raster_bbox_wgs84
    from critical_minerals_aster.structure import (
        annotate_deposits_with_structure,
        load_structure_layers,
        nearest_structure_distance_m,
        points_on_structure,
    )
    from critical_minerals_aster.terrain import compute_hillshade_for_site

    try:
        site = load_site_by_id(site_id, repo_root / "sites")
        paths = site_paths_for(site, repo_root)

        if not paths.strong_zones_geojson.exists():
            return site_id, False, "zones GeoJSON missing"

        zones = gpd.read_file(paths.strong_zones_geojson)

        dem_path = repo_root / "data" / "dem" / site.id / "dem_reprojected.tif"
        hillshade = None
        hs_transform = None
        hs_shape = None
        raster_crs = zones.crs

        if dem_path.exists():
            with rasterio.open(dem_path) as ds:
                raster_transform = ds.transform
                raster_shape = (ds.height, ds.width)
                raster_crs = ds.crs
            hs_result = compute_hillshade_for_site(
                site, paths, raster_transform, raster_shape, raster_crs,
            )
            if hs_result is not None:
                hillshade, hs_transform, hs_shape = hs_result

        deposits = gpd.GeoDataFrame()
        if paths.mrds_csv.exists():
            mrds = read_mrds_national(paths)
            if hs_transform is not None and hs_shape is not None:
                bbox = raster_bbox_wgs84(hs_transform, hs_shape, raster_crs)
            else:
                bbox = site.bbox_wgs84
            local = filter_mrds_bbox(mrds, bbox)
            deposits = mrds_to_points_gdf(local, zones.crs)
            tir_fp = load_tir_footprint_from_provenance(site.id, repo_root, zones.crs)
            if tir_fp is not None and len(tir_fp):
                try:
                    deposits = gpd.clip(deposits, tir_fp)
                except Exception:
                    pass
            joined, _, _ = spatial_join_deposits_zones(deposits, zones)
            hit_ids = joined[joined["index_right"].notna()].index.unique()
            deposits["inside_zone"] = deposits.index.isin(hit_ids)
            deposits["commodity_group"] = deposits["commod1"].apply(simplify_commodity)

        structs = None
        n_on_structure: int | None = None
        target_crs = zones.crs if len(zones) else (deposits.crs if len(deposits) else None)
        if len(deposits) > 0 and target_crs is not None:
            if site.structure_layers:
                structs = load_structure_layers(site, repo_root, target_crs)
                annotated = annotate_deposits_with_structure(
                    deposits, site, paths, structs=structs,
                )
                n_on_structure = int(annotated["on_structure"].sum())
            else:
                structs = auto_fetch_structure(site, repo_root, target_crs)
                if structs is not None and not structs.empty:
                    buffer_m = 500.0
                    annotated = deposits.copy()
                    annotated["nearest_structure_m"] = nearest_structure_distance_m(
                        annotated, structs,
                    )
                    annotated["on_structure"] = points_on_structure(
                        annotated, structs, buffer_m,
                    )
                    n_on_structure = int(annotated["on_structure"].sum())

        tir_footprint = load_tir_footprint_from_provenance(site.id, repo_root, zones.crs)

        basemap_rgb = None
        basemap_source = None
        basemap_cached = None
        if hs_transform is not None and hs_shape is not None:
            bm = fetch_satellite_basemap_for_site(
                site, paths, hs_transform, hs_shape, raster_crs,
            )
            if bm is not None:
                basemap_rgb, basemap_source, basemap_cached = bm

        fig03_meta = save_deposit_overlay_figure(
            site, paths, zones, deposits, repo_root,
            hillshade=hillshade,
            tir_footprint=tir_footprint,
            structs=structs,
            n_on_structure=n_on_structure,
            n_total_deposits=len(deposits) if len(deposits) > 0 else None,
            hs_transform=hs_transform,
            hs_shape=hs_shape,
            basemap_rgb=basemap_rgb,
            basemap_source=basemap_source,
            basemap_cached=basemap_cached,
        )

        prov_path = paths.site_provenance_json
        if prov_path.exists():
            prov = json.loads(prov_path.read_text())
        else:
            prov = {"site_id": site.id}
        prov.update(fig03_meta)
        prov_path.parent.mkdir(parents=True, exist_ok=True)
        prov_path.write_text(json.dumps(prov, indent=2))

        return site_id, True, str(paths.figures_dir / "03_deposit_overlay.png")

    except Exception as exc:
        return site_id, False, str(exc)


def _resolve_site_ids(args: argparse.Namespace, repo_root: Path) -> list[str]:
    if args.all_sites:
        from critical_minerals_aster.config import list_site_ids

        return list_site_ids(repo_root / "sites")
    if args.sites:
        return args.sites
    return sorted(
        p.stem for p in (repo_root / "sites").glob("*.yaml") if p.stem != "index"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Regenerate figure 03 for sites.")
    parser.add_argument("--all-sites", action="store_true", help="All sites in index.yaml")
    parser.add_argument("--sites", nargs="*", help="Specific site IDs")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers (default 4)")
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args()

    site_ids = _resolve_site_ids(args, args.repo_root)
    if not site_ids:
        print("No sites to process.", file=sys.stderr)
        return 1

    workers = max(1, args.workers)
    repo_str = str(args.repo_root.resolve())
    failed: list[str] = []

    if workers == 1 or len(site_ids) == 1:
        for sid in site_ids:
            site_id, ok, msg = _regen_fig3_worker(sid, repo_str)
            status = "OK" if ok else "FAILED"
            print(f"{site_id}: {status} — {msg}")
            if not ok:
                failed.append(site_id)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_regen_fig3_worker, sid, repo_str): sid for sid in site_ids
            }
            for fut in as_completed(futures):
                site_id, ok, msg = fut.result()
                status = "OK" if ok else "FAILED"
                print(f"{site_id}: {status} — {msg}")
                if not ok:
                    failed.append(site_id)

    print(f"\nDone: {len(site_ids) - len(failed)}/{len(site_ids)} succeeded.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
