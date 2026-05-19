"""
Re-run binomial + permutation significance tests excluding Non-Critical deposits.

Compares:
  - all-deposits baseline (existing headline hit rate)
  - critical-minerals-only (earth_mri_category != 'Non-Critical')

Outputs results/significance_critical_only.csv and prints a summary table.
"""

from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.significance import (
    coverage_fraction,
    run_binomial,
    run_permutation,
)


def load_site_bbox(site_id: str) -> tuple:
    yaml_path = ROOT / "sites" / f"{site_id}.yaml"
    with open(yaml_path) as f:
        cfg = yaml.safe_load(f)
    bb = cfg["bbox_wgs84"]
    return tuple(bb)  # [min_lon, min_lat, max_lon, max_lat]


def load_raster_bbox(site_id: str) -> tuple:
    """Return the TIR raster bbox from provenance JSON (rectangular fallback)."""
    import json
    prov_path = ROOT / "results" / f"{site_id}_provenance.json"
    if prov_path.exists():
        with open(prov_path) as f:
            prov = json.load(f)
        rb = prov.get("raster_bbox_wgs84")
        if rb is not None:
            return tuple(rb)
    return load_site_bbox(site_id)


def load_tir_footprint(site_id: str, zones_crs) -> gpd.GeoDataFrame | None:
    """Load the valid-pixel footprint polygon from provenance JSON.

    Returns a GeoDataFrame in *zones_crs*, or None if the footprint is not
    available (old provenance files or computation failure).  The footprint
    polygon reflects the actual diagonal ASTER scene boundary, not just the
    bounding rectangle, so it is a more accurate denominator for coverage_fraction.
    """
    import json
    prov_path = ROOT / "results" / f"{site_id}_provenance.json"
    if not prov_path.exists():
        return None
    with open(prov_path) as f:
        prov = json.load(f)
    geojson = prov.get("tir_footprint_wgs84")
    if geojson is None:
        return None
    try:
        gdf = gpd.GeoDataFrame.from_features(geojson["features"], crs="EPSG:4326")
        return gdf.to_crs(zones_crs)
    except Exception:
        return None


def coverage_fraction_from_footprint(
    zones: gpd.GeoDataFrame,
    footprint: gpd.GeoDataFrame | None,
    bbox: tuple,
) -> float:
    """Compute zone coverage fraction using footprint polygon as denominator.

    Falls back to the rectangular bbox if the footprint is unavailable.
    """
    from shapely.ops import unary_union
    union = unary_union(zones.geometry)
    if footprint is not None and len(footprint):
        _fp = footprint if footprint.crs == zones.crs else footprint.to_crs(zones.crs)
        denom_area = float(_fp.iloc[0].geometry.area)
    else:
        from critical_minerals_aster.significance import _bbox_in_crs
        minx, miny, maxx, maxy = _bbox_in_crs(bbox, zones.crs)
        denom_area = (maxx - minx) * (maxy - miny)
    if denom_area == 0:
        return 0.0
    return min(float(union.area / denom_area), 1.0)


def load_zones(site_id: str) -> gpd.GeoDataFrame | None:
    path = ROOT / "data" / "sites" / site_id / "vectors" / "strong_anomaly_zones.geojson"
    if not path.exists():
        return None
    gdf = gpd.read_file(path)
    if gdf.crs is None or gdf.crs.is_geographic:
        gdf = gdf.to_crs("EPSG:32611")
    return gdf


def main() -> None:
    results_dir = ROOT / "results"
    csvs = sorted(p for p in results_dir.glob("*_summary.csv") if "national" not in p.name)

    # Aggregate per site: all-deposits and critical-only counts
    site_rows: dict[str, dict] = {}
    for path in csvs:
        df = pd.read_csv(path)
        site_id = df["site_id"].iloc[0]

        site_row = df[df["row_type"] == "site"].iloc[0]
        n_all = int(site_row["n_deposits_bbox"])
        hits_all = int(site_row["n_deposits_in_zones"])

        em = df[df["row_type"] == "earth_mri"]
        crit = em[em["earth_mri_category"] != "Non-Critical"]
        n_crit = int(crit["n_deposits_bbox"].sum())
        hits_crit = int(crit["n_deposits_in_zones"].sum())

        site_rows[site_id] = {
            "site_id": site_id,
            "site_name": site_row["site_name"],
            "n_all": n_all,
            "hits_all": hits_all,
            "hr_all": round(hits_all / n_all * 100, 1) if n_all else 0.0,
            "n_crit": n_crit,
            "hits_crit": hits_crit,
            "hr_crit": round(hits_crit / n_crit * 100, 1) if n_crit else 0.0,
        }

    records = []
    for site_id, d in sorted(site_rows.items()):
        print(f"  {site_id} ...", flush=True)
        zones = load_zones(site_id)
        if zones is None or zones.empty:
            print(f"    no zones, skipping")
            continue

        try:
            bbox = load_raster_bbox(site_id)
        except FileNotFoundError:
            print(f"    no YAML, skipping")
            continue

        footprint = load_tir_footprint(site_id, zones.crs)
        p_cover = coverage_fraction_from_footprint(zones, footprint, bbox)

        # All-deposits binomial baseline
        binom_all, _ = run_binomial(d["hits_all"], d["n_all"], p_cover)

        # Critical-only tests
        binom_crit, exp_crit = run_binomial(d["hits_crit"], d["n_crit"], p_cover)
        if d["n_crit"] > 0:
            perm_crit = run_permutation(zones, bbox, d["n_crit"], d["hits_crit"])
        else:
            perm_crit = 1.0

        records.append(
            {
                "site_id": site_id,
                "site_name": d["site_name"],
                "coverage_pct": round(p_cover * 100, 1),
                "n_all": d["n_all"],
                "hits_all": d["hits_all"],
                "hr_all_pct": d["hr_all"],
                "binom_p_all": round(binom_all, 4),
                "n_crit": d["n_crit"],
                "hits_crit": d["hits_crit"],
                "hr_crit_pct": d["hr_crit"],
                "expected_crit": round(exp_crit, 1),
                "binom_p_crit": round(binom_crit, 4),
                "perm_p_crit": round(perm_crit, 4),
                "sig_all": binom_all < 0.05,
                "sig_crit": binom_crit < 0.05,
            }
        )

    out_df = pd.DataFrame(records).sort_values("hr_crit_pct", ascending=False)
    out_path = results_dir / "significance_critical_only.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}\n")

    # Print summary table
    cols = [
        "site_id", "coverage_pct",
        "hr_all_pct", "binom_p_all", "sig_all",
        "hr_crit_pct", "binom_p_crit", "perm_p_crit", "sig_crit",
        "n_crit", "hits_crit",
    ]
    print(out_df[cols].to_string(index=False))

    print("\n=== SIGNIFICANCE CHANGE SUMMARY ===")
    gained = out_df[~out_df["sig_all"] & out_df["sig_crit"]]
    lost = out_df[out_df["sig_all"] & ~out_df["sig_crit"]]
    kept = out_df[out_df["sig_all"] & out_df["sig_crit"]]
    print(f"Significant in BOTH (all & crit): {list(kept['site_id'])}")
    print(f"Gained significance (crit only):  {list(gained['site_id'])}")
    print(f"Lost significance (all only):     {list(lost['site_id'])}")
    sig_crit_all = out_df[out_df["sig_crit"]]
    print(f"\nTotal significant (critical-only test): {len(sig_crit_all)}/{len(out_df)} sites")
    print(sig_crit_all[["site_id", "hr_crit_pct", "binom_p_crit", "perm_p_crit"]].to_string(index=False))


if __name__ == "__main__":
    main()
