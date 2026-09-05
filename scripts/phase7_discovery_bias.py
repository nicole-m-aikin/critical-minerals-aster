#!/usr/bin/env python
"""
Phase 7 -- discovery-bias (ascertainment-bias) re-analysis, corrected.

Two problems with the original scripts/discovery_bias_analysis.py, per the
Phase 1 discovery-bias audit (docs/results.md):

1. It recomputed a per-site coverage fraction inline (`_coverage_fraction`)
   rather than sourcing it from the authoritative Phase 2 table
   (results/site_specific_null_significance.csv), and had already drifted
   from what generates the current significance numbers.
2. It restricted the era-split to only the 12 sites flagged significant
   under the OLD pooled null -- a selection on the outcome variable. The
   audit found the 33 excluded sites contain more post-1950 dated deposits
   (183, bbox-based estimate) than the entire pooled sample it was using
   (35) -- discarding most of the available statistical power for exactly
   the underpowered cohort that most needed it.

Fix: run the pre-/post-1950 split UNCONDITIONALLY across all 45 sites with
at least one dated critical-mineral deposit, using each site's own p0_i
from the Phase 2 table. Report the old-style "previously-flagged-sites-only"
version alongside as a sensitivity comparison -- using the CURRENT
FDR-significant site set (Phase 4's union of binomial + continuous-score
survivors), not the superseded pooled-null 12-site list, since re-using an
already-corrected list from Phase 4 is the only internally consistent
choice at this point in the pipeline.

This remains explicitly a SECONDARY / supplementary analysis, not a
headline claim -- see the caveat text this script prints and
docs/results.md for exactly why: ~14% of critical-mineral
deposits have any usable date at all, dated records are demonstrably skewed
toward historically-producing mines (not a random subsample), and the 1950
cutoff does not cleanly separate "surface-cue-driven" from
"instrument-driven" discovery (post-1950 exploration is still guided by
geological maps built substantially from earlier outcrop mapping).

Outputs
-------
    results/phase7_discovery_bias.csv          -- per-site x cohort table
    results/phase7_discovery_bias_pooled.csv   -- pooled summary (both framings)

Usage:
    conda run -n aster-minerals python scripts/phase7_discovery_bias.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.stats import norm as sp_norm

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.significance import run_binomial

REPO_ROOT = Path(__file__).parent.parent
CUTOFF = 1950


def _best_year(row: pd.Series) -> float:
    dy = row.get("disc_yr")
    if pd.notna(dy):
        return float(dy)
    yfp = row.get("yr_fst_prd")
    if pd.notna(yfp):
        return float(yfp)
    return float("nan")


def _cohort_stats(cohort: pd.DataFrame, label: str, p0: float) -> dict:
    n = len(cohort)
    hits = int(cohort["in_zone"].sum())
    hr = hits / n * 100 if n else 0.0
    p_val, expected = run_binomial(hits, n, p0)
    return {
        f"n_{label}": n,
        f"hits_{label}": hits,
        f"hr_{label}_pct": round(hr, 1),
        f"expected_{label}": round(expected, 2),
        f"binom_p_{label}": round(p_val, 4),
    }


def main() -> None:
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    fdr = pd.read_csv(REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv")
    fdr_sites = set(fdr.loc[fdr["sig_fdr_binomial"] | fdr["sig_fdr_mannwhitney"], "site_id"])
    print(f"Current FDR-significant sites (sensitivity subset, n={len(fdr_sites)}): "
          f"{sorted(fdr_sites)}\n")

    records = []
    n_total_crit, n_total_dated = 0, 0

    for _, ss_row in ss.iterrows():
        site_id = ss_row["site_id"]
        p0 = float(ss_row["p0_null"])

        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        paths = site_paths_for(site, REPO_ROOT)
        zones_path = paths.strong_zones_geojson
        prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
        if not zones_path.exists() or not prov_path.exists():
            continue
        zones = gpd.read_file(zones_path)
        with open(prov_path) as f:
            prov = json.load(f)
        raster_bbox = tuple(prov["raster_bbox_wgs84"])
        footprint = None
        if prov.get("tir_footprint_wgs84") is not None:
            try:
                footprint = gpd.GeoDataFrame.from_features(
                    prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
                ).to_crs(zones.crs)
            except Exception:
                footprint = None

        mrds = read_mrds_national(paths)
        local = filter_mrds_bbox(mrds, raster_bbox)
        deposits = mrds_to_points_gdf(local, zones.crs)
        if footprint is not None and len(footprint):
            try:
                deposits = gpd.clip(deposits, footprint)
            except Exception:
                pass
        deposits = reclassify_mrds_earth_mri(deposits)
        crit = deposits[deposits["earth_mri_category"] != "Non-Critical"].copy()
        if crit.empty:
            continue

        joined, _, _ = spatial_join_deposits_zones(crit, zones)
        hit_ids = set(joined[joined["index_right"].notna()].index.unique())
        crit["in_zone"] = crit.index.isin(hit_ids)
        crit["best_year"] = crit.apply(_best_year, axis=1)

        n_total_crit += len(crit)
        n_total_dated += int(crit["best_year"].notna().sum())

        dated = crit[crit["best_year"].notna()].copy()
        if dated.empty:
            continue
        pre = dated[dated["best_year"] < CUTOFF]
        post = dated[dated["best_year"] >= CUTOFF]

        row = {
            "site_id": site_id,
            "p0_null": p0,
            "n_crit_total": len(crit),
            "n_dated": len(dated),
            "pct_dated": round(len(dated) / len(crit) * 100, 1),
            "fdr_significant_current": site_id in fdr_sites,
            **_cohort_stats(pre, "pre50", p0),
            **_cohort_stats(post, "post50", p0),
        }
        records.append(row)

    out_df = pd.DataFrame(records)
    out_path = REPO_ROOT / "results" / "phase7_discovery_bias.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Saved {out_path}\n")
    print(f"Missingness across all 45 sites: {n_total_dated}/{n_total_crit} "
          f"({n_total_dated/n_total_crit*100:.1f}%) critical-mineral deposits have a usable "
          f"discovery/first-production year.\n")

    def _pooled_test(df: pd.DataFrame, n_col: str, hits_col: str) -> dict:
        obs_hits = df[hits_col].sum()
        exp_list, var_list = [], []
        for _, r in df.iterrows():
            n_i, p_i = r[n_col], r["p0_null"]
            exp_list.append(n_i * p_i)
            var_list.append(n_i * p_i * (1 - p_i))
        exp_hits, var_total = sum(exp_list), sum(var_list)
        if var_total == 0 or df.empty:
            z, p_val = 0.0, 1.0
        else:
            z = (obs_hits - exp_hits) / var_total ** 0.5
            p_val = float(sp_norm.sf(z))
        obs_n = int(df[n_col].sum())
        hr = obs_hits / obs_n * 100 if obs_n else 0.0
        null_hr = exp_hits / obs_n * 100 if obs_n else 0.0
        enrichment = hr / null_hr if null_hr else float("nan")
        return {
            "n": obs_n, "hits": int(obs_hits), "hr_pct": round(hr, 1),
            "null_hr_pct": round(null_hr, 1), "enrichment": round(enrichment, 2),
            "expected_hits": round(exp_hits, 1), "z": round(z, 2), "pooled_p": round(p_val, 4),
        }

    pooled_rows = []
    for label, subset in [
        ("all_45_sites (PRIMARY, this phase)", out_df),
        ("fdr_significant_sites_only (sensitivity)", out_df[out_df["fdr_significant_current"]]),
    ]:
        pre_stats = _pooled_test(subset, "n_pre50", "hits_pre50")
        post_stats = _pooled_test(subset, "n_post50", "hits_post50")
        pooled_rows.append({"framing": label, "cohort": "pre-1950", "n_sites": len(subset), **pre_stats})
        pooled_rows.append({"framing": label, "cohort": "post-1950", "n_sites": len(subset), **post_stats})

    pooled_df = pd.DataFrame(pooled_rows)
    pooled_path = REPO_ROOT / "results" / "phase7_discovery_bias_pooled.csv"
    pooled_df.to_csv(pooled_path, index=False)
    print(f"Saved {pooled_path}\n")
    print(pooled_df.to_string(index=False))

    print(
        "\nCAVEAT (repeat in every citation of this result): this is a SECONDARY / supplementary\n"
        "analysis. ~14% data coverage, dated records skew toward historically-producing mines\n"
        "(not a random subsample -- see docs/results.md), and the 1950 cutoff\n"
        "does not cleanly separate surface-cue-driven from instrument-driven discovery, since\n"
        "modern exploration still builds on earlier outcrop/alteration mapping. Do NOT claim\n"
        "discovery bias has been ruled out, eliminated, or confirmed absent."
    )


if __name__ == "__main__":
    main()
