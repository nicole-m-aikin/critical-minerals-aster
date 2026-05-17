#!/usr/bin/env python
"""
Per-(site, Earth MRI category) significance tests.

For each site × Earth MRI category cell with n ≥ MIN_N deposits, runs a
one-sided binomial test asking whether that category's deposits hit the
anomaly zones more than expected by zone area alone.

Also reports a national roll-up per category: pooled hits vs. pooled expected
across all sites (each deposit weighted by its site's coverage fraction).

Usage:
    conda run -n aster-minerals python scripts/compute_significance_by_category.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import yaml
import geopandas as gpd
import pandas as pd
import numpy as np
from scipy.stats import binomtest

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import SitePaths
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    mrds_to_points_gdf,
    spatial_join_deposits_zones,
    reclassify_mrds_earth_mri,
)
from critical_minerals_aster.significance import coverage_fraction, run_binomial

REPO_ROOT = Path(__file__).parent.parent
MIN_N = 10   # minimum deposits per cell to include in per-site table

# Expected TIR detectability per Earth MRI category
TIR_EXPECTED = {
    "Base Metals":              "YES  — porphyry/skarn/VMS alteration halos",
    "REE":                      "YES  — carbonate signature (B13/B12)",
    "PGM":                      "YES  — mafic ratio (B12/B13)",
    "Mafic Magmatic":           "YES  — mafic ratio",
    "Battery Metals – Co/Ni":   "MAYBE — some magmatic/skarn settings",
    "Specialty/High-Tech":      "MAYBE — deposit-type dependent",
    "Gold/Silver":              "NO   — placer/Carlin/epithermal, buried/eroded",
    "Energy":                   "NO   — roll-front U, coal, no alteration",
    "Battery Metals – Li/Brine":"NO   — evaporite/brine hosted",
    "Industrial":               "NO   — sand/gravel/stone",
    "Non-Critical":             "NO   — mixed non-critical",
}


def load_site_ids() -> list[str]:
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        return yaml.safe_load(f)["sites"]


def collect_all_deposits() -> list[dict]:
    """Load deposits + zone join for all sites, tagged with site coverage fraction."""
    site_ids = load_site_ids()
    records = []

    for site_id in site_ids:
        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        paths = SitePaths(repo_root=REPO_ROOT, site=site)

        zones_path = paths.strong_zones_geojson
        if not zones_path.exists():
            continue
        zones = gpd.read_file(zones_path)
        if zones.empty:
            continue

        try:
            mrds_df = read_mrds_national(paths)
        except FileNotFoundError:
            continue

        local = filter_mrds_bbox(mrds_df, site.bbox_wgs84)
        if local.empty:
            continue

        local = reclassify_mrds_earth_mri(local)
        deposits = mrds_to_points_gdf(local, zones.crs)
        joined, _, _ = spatial_join_deposits_zones(deposits, zones)
        inside_ids = set(joined[joined["index_right"].notna()].index.unique())

        p_cov = coverage_fraction(zones, site.bbox_wgs84)

        for idx, row in deposits.iterrows():
            records.append({
                "site_id": site_id,
                "dep_id": idx,
                "earth_mri_category": row.get("earth_mri_category", "Unknown"),
                "in_zone": idx in inside_ids,
                "site_coverage": p_cov,
            })

        print(f"  {site_id}: {len(deposits)} deposits, coverage={p_cov:.1%}")

    return records


def per_site_category_table(df: pd.DataFrame, min_n: int = MIN_N) -> pd.DataFrame:
    rows = []
    for (site_id, cat), grp in df.groupby(["site_id", "earth_mri_category"]):
        n = len(grp)
        if n < min_n:
            continue
        n_hit = int(grp["in_zone"].sum())
        p_cov = grp["site_coverage"].iloc[0]
        p_binom, expected = run_binomial(n_hit, n, p_cov)
        rows.append({
            "site_id": site_id,
            "earth_mri_category": cat,
            "n": n,
            "n_hits": n_hit,
            "hit_rate_pct": round(n_hit / n * 100, 1),
            "zone_cov_pct": round(p_cov * 100, 1),
            "expected": round(expected, 1),
            "p_binom": round(p_binom, 4),
        })
    return pd.DataFrame(rows).sort_values(["earth_mri_category", "p_binom"])


def national_category_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pooled test per category across all sites.

    Each deposit has expected probability = its site's coverage fraction.
    Total expected hits = sum of site_coverage over all deposits in category.
    Observed hits = actual count in zone.
    Ratio = observed / expected (>1 means better than random).

    For the pooled p-value, approximate the sum of independent Bernoullis
    (different p_i per deposit) using the normal approximation to the
    Poisson binomial distribution.
    """
    rows = []
    for cat, grp in df.groupby("earth_mri_category"):
        n = len(grp)
        n_hit = int(grp["in_zone"].sum())
        p_i = grp["site_coverage"].values
        expected = p_i.sum()
        variance = (p_i * (1 - p_i)).sum()
        # Normal approximation (continuity correction)
        z = (n_hit - 0.5 - expected) / np.sqrt(variance) if variance > 0 else 0.0
        from scipy.stats import norm
        p_pooled = float(norm.sf(z))
        ratio = n_hit / expected if expected > 0 else float("nan")
        rows.append({
            "earth_mri_category": cat,
            "tir_expected": TIR_EXPECTED.get(cat, "?"),
            "n_deposits": n,
            "n_hits": n_hit,
            "expected_hits": round(expected, 1),
            "hit_vs_expected": round(ratio, 2),
            "p_pooled": round(p_pooled, 4),
        })
    return (
        pd.DataFrame(rows)
        .sort_values("p_pooled")
        .reset_index(drop=True)
    )


def sig(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "** "
    if p < 0.05:  return "*  "
    return "   "


def main() -> None:
    print("Loading deposits for all sites...")
    records = collect_all_deposits()
    df = pd.DataFrame(records)

    # ── National roll-up ──────────────────────────────────────────────────
    print("\n" + "═" * 90)
    print("NATIONAL: pooled hits vs expected, by Earth MRI category")
    print("═" * 90)
    nat = national_category_table(df)
    for _, row in nat.iterrows():
        print(
            f"  {row['earth_mri_category']:28s}"
            f"  n={row['n_deposits']:5d}"
            f"  hits={row['n_hits']:4d}"
            f"  exp={row['expected_hits']:6.1f}"
            f"  ratio={row['hit_vs_expected']:5.2f}"
            f"  p={row['p_pooled']:.4f}{sig(row['p_pooled'])}"
            f"  TIR: {row['tir_expected']}"
        )

    # ── Per-site × category ───────────────────────────────────────────────
    print("\n" + "═" * 90)
    print(f"PER-SITE × CATEGORY  (min n={MIN_N} deposits per cell)")
    print("═" * 90)
    site_cat = per_site_category_table(df)

    for cat, grp in site_cat.groupby("earth_mri_category"):
        tir = TIR_EXPECTED.get(cat, "?")
        print(f"\n  {cat}  [{tir}]")
        for _, row in grp.iterrows():
            ratio = row["n_hits"] / row["expected"] if row["expected"] > 0 else float("nan")
            direction = "▲" if ratio > 1.05 else ("▼" if ratio < 0.95 else "─")
            print(
                f"    {row['site_id']:22s}"
                f"  n={row['n']:4d}"
                f"  hits={row['n_hits']:3d}"
                f"  exp={row['expected']:5.1f}"
                f"  {direction} ratio={ratio:.2f}"
                f"  p={row['p_binom']:.4f}{sig(row['p_binom'])}"
            )

    # Save
    out_nat = REPO_ROOT / "results" / "significance_by_category_national.csv"
    out_site = REPO_ROOT / "results" / "significance_by_category_per_site.csv"
    nat.to_csv(out_nat, index=False)
    site_cat.to_csv(out_site, index=False)
    print(f"\nWrote {out_nat}")
    print(f"Wrote {out_site}")
    print("\n* p<0.05  ** p<0.01  *** p<0.001")


if __name__ == "__main__":
    main()
