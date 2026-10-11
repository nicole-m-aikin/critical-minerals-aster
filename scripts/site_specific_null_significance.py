#!/usr/bin/env python
"""
Phase 2 — primary significance table: BEDROCK-vs-COVER STRATIFIED null.

Why stratify
------------
Deposits sit in bedrock. ASTER TIR strong-anomaly zones move between outcrop
and alluvial valley fill depending on acquisition temperature (verified when
the archive mosaics were rebuilt from clear-day scenes), so a purely
geometric null — every footprint pixel equally likely to host a deposit —
rewards any detector that preferentially flags outcrop, and penalises one
that flags valley fill, independently of alteration. Both the original
enrichment and its later collapse tracked that bedrock/cover shift.

Primary null (per site)
-----------------------
The TIR footprint is partitioned into three strata
(``critical_minerals_aster.lithology.cover_strata``, USGS SGMC polygons):

    cover    = SGMC GENERALIZED_LITH "Unconsolidated…" (alluvium, valley fill)
    unmapped = footprint with no SGMC polygon (e.g. Mexico)
    bedrock  = everything else

Each critical-mineral deposit i gets the zone coverage fraction of its OWN
stratum, p_i = area(zone ∩ stratum) / area(stratum). Expected hits
E = Σ p_i; p_binomial is the exact Poisson-binomial tail P(X ≥ hits).
p0_null = E / n (the hit rate the null predicts), so enrichment = hits / E
keeps its meaning and every downstream consumer of p0_null / enrichment /
p_binomial reads the stratified values without code changes.

The previous geometric null (p0 = zone area / footprint area, exact binomial)
is kept as the ``*_geometric`` columns for comparison.

Inputs mirror scripts/regenerate_site_summaries.py (validated MRDS → raster
bbox → TIR footprint clip → critical minerals); n / hits are cross-checked
against the regenerated results/{site}_summary.csv.

Outputs
-------
    results/site_specific_null_significance.csv

Columns
-------
    site_id, site_name, null_model
    footprint_area_km2, zone_area_km2
    p0_null, expected_hits_null, enrichment, p_binomial      -- PRIMARY (stratified)
    n_crit, hits_crit, hit_rate_crit_pct, ci_low_pct, ci_high_pct
    p0_geometric, expected_hits_geometric, enrichment_geometric, p_binomial_geometric
    {bedrock,cover,unmapped}_frac_footprint                  -- stratum share of footprint
    p0_{bedrock,cover,unmapped}                              -- zone coverage within stratum
    zone_frac_on_cover, zone_frac_on_unmapped                -- share of zone area
    n_crit_{bedrock,cover,unmapped}, hits_{bedrock,cover,unmapped}
    deposit_frac_on_cover, deposit_frac_on_unmapped
    n_all, hits_all, hit_rate_all_pct                        -- all deposits (incl. non-critical)
    pooled_null_pct, expected_hits_pooled, p_binomial_pooled -- HISTORICAL comparison only
    sig_site_specific_05, sig_geometric_05, sig_pooled_05

Usage:
    conda run -n aster-minerals python scripts/site_specific_null_significance.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import numpy as np
import pandas as pd
import yaml
from shapely.geometry import box
from shapely.ops import unary_union
from shapely.validation import make_valid

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.lithology import (
    COVER,
    COVER_STRATA,
    UNMAPPED,
    assign_cover_stratum,
    cover_strata,
    stratified_null,
    stratum_zone_fractions,
)
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    is_critical_mineral,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.significance import (
    coverage_fraction_footprint,
    run_binomial,
    wilson_ci,
)

REPO_ROOT = Path(__file__).parent.parent
NULL_MODEL = "bedrock_cover_stratified"


def load_site_ids() -> list[str]:
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        return yaml.safe_load(f)["sites"]


def load_site_geometry(site_id: str):
    """Zones, footprint GeoDataFrame (or None), footprint geometry, raster bbox,
    critical deposits with ``hit`` — same deposit set as regenerate_site_summaries."""
    site = load_site_by_id(site_id, REPO_ROOT / "sites")
    paths = site_paths_for(site, REPO_ROOT)
    prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
    if not paths.strong_zones_geojson.exists() or not prov_path.exists():
        raise FileNotFoundError("missing zones or provenance")
    zones = gpd.read_file(paths.strong_zones_geojson)
    if zones.empty:
        raise ValueError("empty zones")
    prov = json.loads(prov_path.read_text())
    raster_bbox = tuple(prov["raster_bbox_wgs84"]) if prov.get("raster_bbox_wgs84") else site.bbox_wgs84

    footprint = None
    if prov.get("tir_footprint_wgs84") is not None:
        footprint = gpd.GeoDataFrame.from_features(
            prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
        ).to_crs(zones.crs)
    if footprint is not None and len(footprint):
        fp_geom = make_valid(unary_union(footprint.geometry.values))
    else:
        fp_geom = gpd.GeoSeries([box(*raster_bbox)], crs="EPSG:4326").to_crs(zones.crs).iloc[0]

    deposits = mrds_to_points_gdf(filter_mrds_bbox(read_mrds_national(paths), raster_bbox), zones.crs)
    if footprint is not None and len(footprint):
        deposits = gpd.clip(deposits, footprint)
    deposits = reclassify_mrds_earth_mri(deposits)
    deposits = deposits[deposits["earth_mri_category"].apply(is_critical_mineral)].copy()
    _, hits, _ = spatial_join_deposits_zones(deposits, zones)
    deposits["hit"] = deposits.index.isin(hits.index.unique())
    return site, zones, footprint, fp_geom, raster_bbox, deposits


def stratified_site(site_id: str, zones, fp_geom, raster_bbox, deposits):
    """(strata, per-stratum fractions table, per-deposit stratum labels, per-deposit p_i)."""
    strata = cover_strata(site_id, fp_geom, zones.crs, bbox_wgs84=raster_bbox)
    zone_union = make_valid(unary_union(zones.geometry.values))
    frac = stratum_zone_fractions(strata, zone_union)
    labels = assign_cover_stratum(deposits, strata) if len(deposits) else pd.Series(dtype=object)
    p_i = labels.map(frac["p"]).to_numpy(dtype=float)
    return strata, frac, labels, p_i


def main() -> None:
    site_ids = load_site_ids()
    rows: list[dict] = []
    skipped: list[tuple[str, str]] = []
    mismatches: list[tuple] = []

    for site_id in site_ids:
        path = REPO_ROOT / "results" / f"{site_id}_summary.csv"
        if not path.exists():
            skipped.append((site_id, "no regenerated summary"))
            continue
        df = pd.read_csv(path)
        if "null_hit_rate_pct" not in df.columns:
            skipped.append((site_id, "summary predates add_uncertainty_columns"))
            continue
        try:
            site, zones, footprint, fp_geom, raster_bbox, deposits = load_site_geometry(site_id)
        except Exception as exc:
            skipped.append((site_id, repr(exc)))
            continue

        site_row = df[df["row_type"] == "site"].iloc[0]
        em = df[df["row_type"] == "earth_mri"]
        crit = em[em["earth_mri_category"] != "Non-Critical"]
        n_ref, hits_ref = int(crit["n_deposits_bbox"].sum()), int(crit["n_deposits_in_zones"].sum())

        n_crit = len(deposits)
        hits_crit = int(deposits["hit"].sum())
        if (n_crit, hits_crit) != (n_ref, hits_ref):
            mismatches.append((site_id, n_crit, n_ref, hits_crit, hits_ref))
        hit_rate_crit = round(hits_crit / n_crit * 100, 2) if n_crit else 0.0
        ci_lo, ci_hi = wilson_ci(hits_crit, n_crit)

        # --- Geometric null (previous primary) ---
        p0_geo = coverage_fraction_footprint(zones, footprint, raster_bbox)
        p_geo, exp_geo = run_binomial(hits_crit, n_crit, p0_geo)

        # --- Bedrock / cover / unmapped stratified null (primary) ---
        _, frac, labels, p_i = stratified_site(site_id, zones, fp_geom, raster_bbox, deposits)
        strat = stratified_null(p_i, hits_crit)

        fp_area = float(frac["area"].sum())
        zone_fp_area = float(frac["zone_area"].sum())
        hit_arr = deposits["hit"].to_numpy(dtype=bool)
        lab_arr = labels.to_numpy(dtype=object)

        row = {
            "site_id": site_id,
            "site_name": site_row["site_name"],
            "null_model": NULL_MODEL,
            "footprint_area_km2": round(float(fp_geom.area) / 1e6, 1),
            "zone_area_km2": round(float(site_row["total_anomaly_km2"]), 2),
            "p0_null": round(strat["p0_eff"], 4),
            "n_crit": n_crit,
            "hits_crit": hits_crit,
            "hit_rate_crit_pct": hit_rate_crit,
            "expected_hits_null": round(strat["expected"], 1),
            "enrichment": round(strat["enrichment"], 2) if strat["expected"] > 0 else float("nan"),
            "p_binomial": strat["p"],
            "ci_low_pct": round(ci_lo * 100, 2),
            "ci_high_pct": round(ci_hi * 100, 2),
            "p0_geometric": round(p0_geo, 4),
            "expected_hits_geometric": round(exp_geo, 1),
            "enrichment_geometric": round((hits_crit / n_crit) / p0_geo, 2) if n_crit and p0_geo > 0 else float("nan"),
            "p_binomial_geometric": p_geo,
        }
        for s in COVER_STRATA:
            row[f"{s}_frac_footprint"] = round(frac.loc[s, "area"] / fp_area, 4) if fp_area > 0 else np.nan
        for s in COVER_STRATA:
            row[f"p0_{s}"] = round(frac.loc[s, "p"], 4)
        for s in (COVER, UNMAPPED):
            row[f"zone_frac_on_{s}"] = round(frac.loc[s, "zone_area"] / zone_fp_area, 4) if zone_fp_area > 0 else np.nan
        for s in COVER_STRATA:
            in_s = lab_arr == s
            row[f"n_crit_{s}"] = int(in_s.sum())
            row[f"hits_{s}"] = int((in_s & hit_arr).sum())
        for s in (COVER, UNMAPPED):
            row[f"deposit_frac_on_{s}"] = round(row[f"n_crit_{s}"] / n_crit, 4) if n_crit else np.nan
        row.update({
            "n_all": int(site_row["n_deposits_bbox"]),
            "hits_all": int(site_row["n_deposits_in_zones"]),
            "hit_rate_all_pct": float(site_row["hit_rate_pct"]),
        })
        rows.append(row)
        print(
            f"{site_id:24s} n={n_crit:4d} hits={hits_crit:3d}  "
            f"geo: p0={p0_geo:.3f} E={row['enrichment_geometric']} p={p_geo:.3g}  "
            f"strat: p0={strat['p0_eff']:.3f} E={row['enrichment']} p={strat['p']:.3g}  "
            f"cover fp/zone/dep={row['cover_frac_footprint']:.2f}/"
            f"{row['zone_frac_on_cover']:.2f}/{row['deposit_frac_on_cover']}"
        )

    out = pd.DataFrame(rows)

    # Historical pooled-null comparison: same pooled rate the paper previously
    # cited (sum of hits / sum of n across all sites), applied uniformly.
    pooled_p = out["hits_crit"].sum() / out["n_crit"].sum()
    out["pooled_null_pct"] = round(pooled_p * 100, 2)
    out["expected_hits_pooled"] = (out["n_crit"] * pooled_p).round(1)
    out["p_binomial_pooled"] = [
        round(run_binomial(int(h), int(n), pooled_p)[0], 6)
        for h, n in zip(out["hits_crit"], out["n_crit"])
    ]

    out["sig_site_specific_05"] = out["p_binomial"] < 0.05
    out["sig_geometric_05"] = out["p_binomial_geometric"] < 0.05
    out["sig_pooled_05"] = out["p_binomial_pooled"] < 0.05

    out = out.sort_values("p_binomial")
    out_path = REPO_ROOT / "results" / "site_specific_null_significance.csv"
    out.to_csv(out_path, index=False)

    print(f"\nSaved {out_path}\n")
    if skipped:
        print(f"Skipped ({len(skipped)}): {skipped}\n")
    if mismatches:
        print(f"n/hits differ from results/*_summary.csv at {len(mismatches)} sites "
              f"(site, n, n_ref, hits, hits_ref): {mismatches}\n")
    else:
        print(f"n/hits match results/*_summary.csv at all {len(out)} sites.\n")

    print(f"Historical pooled null (critical-mineral, this run): {pooled_p*100:.2f}%")
    print(f"Stratified p0 range: {out['p0_null'].min()*100:.2f}% - {out['p0_null'].max()*100:.2f}%  "
          f"(geometric: {out['p0_geometric'].min()*100:.2f}% - {out['p0_geometric'].max()*100:.2f}%)")
    print(f"Median cover share — footprint {out['cover_frac_footprint'].median():.2f}, "
          f"zone area {out['zone_frac_on_cover'].median():.2f}, "
          f"critical deposits {out['deposit_frac_on_cover'].median():.2f}\n")

    cols = ["site_id", "n_crit", "hits_crit", "p0_geometric", "p0_null",
            "enrichment_geometric", "enrichment", "p_binomial_geometric", "p_binomial",
            "cover_frac_footprint", "zone_frac_on_cover", "deposit_frac_on_cover"]
    print(out[cols].to_string(index=False))

    print(f"\nSignificant (stratified null, p<0.05): {int(out['sig_site_specific_05'].sum())}/{len(out)}")
    print(f"Significant (geometric null, p<0.05):  {int(out['sig_geometric_05'].sum())}/{len(out)}")
    print(f"Significant (pooled null, p<0.05):     {int(out['sig_pooled_05'].sum())}/{len(out)}")


if __name__ == "__main__":
    main()
