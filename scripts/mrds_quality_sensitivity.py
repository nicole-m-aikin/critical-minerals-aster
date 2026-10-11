#!/usr/bin/env python
"""
MRDS record-quality sensitivity — does the headline signal depend on
low-quality legacy records?

CIM (2018) §2.1.1.1 asks that legacy data be verified for "its degree of
quality before being included in the project database". MRDS carries two
quality indicators that ``validate_mrds()`` exposes as flags:

    qc_score_high     USGS record-quality grade A–C (ScienceBase
                      5d712a5ae4b0c4f70cfd4e54; D = summary only, E = generally
                      no bibliographic reference)
    qc_producer       dev_stat Producer / Past Producer (located on a working,
                      so usually better positioned than a reported occurrence)
    qc_shared_coords  shares exact coordinates with a differently-named record

For every site the critical-mineral binomial test (same definition as
scripts/site_specific_null_significance.py: validated MRDS, clipped to the
valid TIR footprint, critical minerals only, null p0 = zone area / footprint
area) is re-run on each quality subset, then BH-FDR (alpha 0.05) is applied
within each subset across all sites.

Outputs
-------
    results/mrds_quality_sensitivity.csv   one row per site x subset

Usage
-----
    conda run -n aster-minerals python scripts/mrds_quality_sensitivity.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    is_critical_mineral,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.significance import run_binomial

ALPHA = 0.05

SUBSETS: dict[str, str | None] = {
    "all_validated": None,
    "score_A_to_C": "qc_score_high",
    "producers": "qc_producer",
    "score_A_to_C_and_producer": "qc_score_high & qc_producer",
    "unique_coords_only": "~qc_shared_coords",
}


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg adjusted p-values (NaN-safe)."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = len(pv)
    if m == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def _subset_mask(df: pd.DataFrame, expr: str | None) -> pd.Series:
    if expr is None:
        return pd.Series(True, index=df.index)
    mask = pd.Series(True, index=df.index)
    for term in expr.split("&"):
        term = term.strip()
        if term.startswith("~"):
            mask &= ~df[term[1:]].astype(bool)
        else:
            mask &= df[term].astype(bool)
    return mask


def site_deposits(site_id: str) -> tuple[gpd.GeoDataFrame, float] | None:
    site = load_site_by_id(site_id, REPO_ROOT / "sites")
    paths = site_paths_for(site, REPO_ROOT)
    prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
    if not paths.strong_zones_geojson.exists() or not prov_path.exists():
        return None
    zones = gpd.read_file(paths.strong_zones_geojson)
    if zones.empty:
        return None
    prov = json.loads(prov_path.read_text())
    raster_bbox = tuple(prov.get("raster_bbox_wgs84") or site.bbox_wgs84)
    footprint = None
    if prov.get("tir_footprint_wgs84") is not None:
        footprint = gpd.GeoDataFrame.from_features(
            prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
        ).to_crs(zones.crs)

    local = filter_mrds_bbox(read_mrds_national(paths), raster_bbox)
    deposits = mrds_to_points_gdf(local, zones.crs)
    if footprint is not None and len(footprint):
        deposits = gpd.clip(deposits, footprint)
    deposits = reclassify_mrds_earth_mri(deposits)
    deposits = deposits[deposits["earth_mri_category"].map(is_critical_mineral)].copy()
    joined, _, _ = spatial_join_deposits_zones(deposits, zones)
    hit_ids = set(joined.index[joined["index_right"].notna()])
    deposits["inside_zone"] = deposits.index.isin(hit_ids)

    union_area = zones.geometry.union_all().area if hasattr(zones.geometry, "union_all") else zones.unary_union.area
    fp_area = footprint.geometry.union_all().area if footprint is not None else float("nan")
    p0 = union_area / fp_area if fp_area and fp_area > 0 else float("nan")
    return deposits, p0


def main() -> None:
    site_ids = yaml.safe_load((REPO_ROOT / "sites" / "index.yaml").read_text())["sites"]
    ref = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv").set_index("site_id")

    rows: list[dict] = []
    for site_id in site_ids:
        res = site_deposits(site_id)
        if res is None:
            print(f"{site_id:24s} SKIP")
            continue
        deposits, p0_geom = res
        if deposits.empty:
            print(f"{site_id:24s} SKIP (no critical-mineral deposits)")
            continue
        p0 = float(ref.loc[site_id, "p0_null"]) if site_id in ref.index else p0_geom
        for name, expr in SUBSETS.items():
            sub = deposits[_subset_mask(deposits, expr)]
            n, hits = len(sub), int(sub["inside_zone"].sum())
            p, expected = run_binomial(hits, n, p0)
            rows.append({
                "site_id": site_id,
                "subset": name,
                "n_crit": n,
                "hits_crit": hits,
                "hit_rate_pct": round(100 * hits / n, 2) if n else np.nan,
                "p0_null": round(p0, 4),
                "p0_geometry_check": round(p0_geom, 4),
                "enrichment": round((hits / n) / p0, 2) if n and p0 > 0 else np.nan,
                "p_binomial": p if n else np.nan,
            })
        print(f"{site_id:24s} n_all={len(deposits):5d}", flush=True)

    out = pd.DataFrame(rows)
    out["q_binomial"] = np.nan
    for name in SUBSETS:
        m = out["subset"] == name
        out.loc[m, "q_binomial"] = bh_qvalues(out.loc[m, "p_binomial"].to_numpy())
    out["sig_fdr"] = out["q_binomial"] < ALPHA

    out_path = REPO_ROOT / "results" / "mrds_quality_sensitivity.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved {out_path.relative_to(REPO_ROOT)}")

    wide = out.pivot(index="site_id", columns="subset", values="sig_fdr")
    print("\nFDR-significant sites per subset:")
    for name in SUBSETS:
        print(f"  {name:28s} {int(wide[name].sum()):3d}")
    base = set(wide.index[wide["all_validated"]])
    for name in SUBSETS:
        if name == "all_validated":
            continue
        lost = sorted(base - set(wide.index[wide[name]]))
        print(f"  lost under {name}: {lost}")


if __name__ == "__main__":
    main()
