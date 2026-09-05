#!/usr/bin/env python
"""
Phase 2 — primary significance table using the SITE-SPECIFIC null.

Null hypothesis, per site i:
    p0_i = anomaly_zone_area_i / valid_ASTER_footprint_area_i
i.e. the probability a uniformly-randomly-placed deposit lands inside a
strong-anomaly zone, given *that site's own* zone geometry — not a rate
pooled across all 45 sites. See docs/results.md for the
full derivation and why the previously-reported pooled null (8.5% critical /
9.5% all-deposit) is not a scientifically defensible denominator once
per-site coverage is shown to range 9.0-14.7%.

Reads the site summary CSVs written by scripts/regenerate_site_summaries.py
(which already ran add_uncertainty_columns() and therefore carries
null_hit_rate_pct = p0_i * 100 on every row of every site). This script's
only job is to aggregate the earth_mri sub-rows into a critical-vs-non-critical
split per site and compute the requested Phase 2 columns.

Outputs
-------
    results/site_specific_null_significance.csv

Columns
-------
    site_id, site_name
    footprint_area_km2, zone_area_km2, p0_null            -- site geometry
    n_crit, hits_crit, hit_rate_crit_pct                   -- critical-mineral observed
    expected_hits_null, enrichment                         -- vs site-specific null
    p_binomial, ci_low_pct, ci_high_pct                     -- Wilson CI on observed rate
    pooled_null_pct, expected_hits_pooled, p_binomial_pooled  -- HISTORICAL comparison only
    n_all, hits_all, hit_rate_all_pct                       -- all-deposit (incl. non-critical)

Usage:
    conda run -n aster-minerals python scripts/site_specific_null_significance.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pandas as pd
import yaml

from critical_minerals_aster.significance import run_binomial, wilson_ci

REPO_ROOT = Path(__file__).parent.parent


def load_site_ids() -> list[str]:
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        return yaml.safe_load(f)["sites"]


def main() -> None:
    site_ids = load_site_ids()
    rows: list[dict] = []
    skipped: list[str] = []

    for site_id in site_ids:
        path = REPO_ROOT / "results" / f"{site_id}_summary.csv"
        if not path.exists():
            skipped.append(site_id)
            continue
        df = pd.read_csv(path)
        if "null_hit_rate_pct" not in df.columns:
            skipped.append(site_id)
            continue

        site_row = df[df["row_type"] == "site"].iloc[0]
        p0 = float(site_row["null_hit_rate_pct"]) / 100.0
        zone_area_km2 = float(site_row["total_anomaly_km2"])
        footprint_area_km2 = zone_area_km2 / p0 if p0 > 0 else float("nan")

        em = df[df["row_type"] == "earth_mri"]
        crit = em[em["earth_mri_category"] != "Non-Critical"]
        n_crit = int(crit["n_deposits_bbox"].sum())
        hits_crit = int(crit["n_deposits_in_zones"].sum())
        hit_rate_crit = round(hits_crit / n_crit * 100, 2) if n_crit else 0.0

        p_binom, expected = run_binomial(hits_crit, n_crit, p0)
        ci_lo, ci_hi = wilson_ci(hits_crit, n_crit)
        enrichment = round((hit_rate_crit / 100) / p0, 2) if p0 > 0 else float("nan")

        rows.append({
            "site_id": site_id,
            "site_name": site_row["site_name"],
            "footprint_area_km2": round(footprint_area_km2, 1),
            "zone_area_km2": round(zone_area_km2, 2),
            "p0_null": round(p0, 4),
            "n_crit": n_crit,
            "hits_crit": hits_crit,
            "hit_rate_crit_pct": hit_rate_crit,
            "expected_hits_null": round(expected, 1),
            "enrichment": enrichment,
            "p_binomial": round(p_binom, 6),
            "ci_low_pct": round(ci_lo * 100, 2),
            "ci_high_pct": round(ci_hi * 100, 2),
            "n_all": int(site_row["n_deposits_bbox"]),
            "hits_all": int(site_row["n_deposits_in_zones"]),
            "hit_rate_all_pct": float(site_row["hit_rate_pct"]),
        })

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
    out["sig_pooled_05"] = out["p_binomial_pooled"] < 0.05

    out = out.sort_values("p_binomial")
    out_path = REPO_ROOT / "results" / "site_specific_null_significance.csv"
    out.to_csv(out_path, index=False)

    print(f"Saved {out_path}\n")
    if skipped:
        print(f"Skipped (no regenerated summary found): {skipped}\n")

    print(f"Historical pooled null (critical-mineral, this run): {pooled_p*100:.2f}%")
    print(f"Site-specific null range: {out['p0_null'].min()*100:.2f}% - {out['p0_null'].max()*100:.2f}%\n")

    cols = ["site_id", "p0_null", "n_crit", "hits_crit", "hit_rate_crit_pct",
            "enrichment", "p_binomial", "sig_site_specific_05",
            "p_binomial_pooled", "sig_pooled_05"]
    print(out[cols].to_string(index=False))

    n_sig_ss = int(out["sig_site_specific_05"].sum())
    n_sig_pooled = int(out["sig_pooled_05"].sum())
    flipped = out[out["sig_site_specific_05"] != out["sig_pooled_05"]]
    print(f"\nSignificant (site-specific null, p<0.05): {n_sig_ss}/{len(out)}")
    print(f"Significant (pooled null, p<0.05):        {n_sig_pooled}/{len(out)}")
    if not flipped.empty:
        print(f"\nSites where the two nulls disagree on significance:")
        print(flipped[["site_id", "p_binomial", "p_binomial_pooled",
                        "sig_site_specific_05", "sig_pooled_05"]].to_string(index=False))


if __name__ == "__main__":
    main()
