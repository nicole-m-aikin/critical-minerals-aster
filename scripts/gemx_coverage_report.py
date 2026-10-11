"""Report GEMx / AVIRIS hyperspectral coverage for every study site.

Answers the prerequisite question for any hyperspectral extension of this
project: which of the 109 sites actually sit under an imaging-spectrometer
flight line, and with which sensor?

State-level coverage claims ("California, Nevada, Arizona, New Mexico") badly
overstate usable overlap, because GEMx flies discrete numbered boxes rather
than whole states.  This script resolves coverage per site by intersecting
real CMR granule footprints with each site bbox.

The report is oriented toward the SWIR detection gap named in
docs/results.md limitation 2: ASTER SWIR is unavailable here, so argillic /
phyllic alteration is invisible to the TIR method, which may explain the
below-chance epithermal and sediment-hosted medians.  AVIRIS resolves the
2.2 um Al-OH features that identify those minerals, so the sites that matter
most are the epithermal and porphyry nulls.

Outputs
-------
results/gemx_coverage.csv          — long form, one row per (site, collection)
results/gemx_coverage_by_site.csv  — one decision row per site

Usage
-----
    conda run -n aster-minerals python scripts/gemx_coverage_report.py

Requires no EarthData credentials: CMR metadata search is public.  Raw CMR
responses are cached under data/gemx/cmr_cache/ so re-runs are offline.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.config import list_site_ids, load_site_by_id
from critical_minerals_aster.hyperspectral import (
    COLLECTIONS,
    build_coverage_table,
    summarize_by_site,
)

SITES_DIR = ROOT / "sites"
RESULTS_DIR = ROOT / "results"
CACHE_DIR = ROOT / "data" / "gemx" / "cmr_cache"


def load_site_context() -> pd.DataFrame:
    """Per-site name, deposit-system category and significance tier.

    The category and tier columns are what make the coverage table
    actionable: they say whether a covered site is one of the epithermal /
    porphyry nulls that the SWIR-gap test needs.
    """
    rows = []
    for site_id in list_site_ids(SITES_DIR):
        cfg = load_site_by_id(site_id, SITES_DIR)
        rows.append({"site_id": site_id, "site_name": cfg.name})
    ctx = pd.DataFrame(rows)
    # Site names carry the state as a trailing suffix in three forms: single
    # ("Bisbee, AZ"), slash for a site straddling a border ("McDermitt
    # Caldera, NV/OR") and hyphen for tri-state districts ("..., MO-KS-OK").
    ctx["state"] = ctx["site_name"].str.extract(
        r",\s*([A-Z]{2}(?:[-/][A-Z]{2})*)\s*$"
    )[0]

    cat_path = RESULTS_DIR / "deposit_system_categories.csv"
    if cat_path.is_file():
        cats = pd.read_csv(cat_path).set_index("site_id")["category"]
        ctx["category"] = ctx["site_id"].map(cats)

    fdr_path = RESULTS_DIR / "phase4_fdr_corrected_significance.csv"
    clust_path = RESULTS_DIR / "phase5_clustering_sensitivity.csv"
    if fdr_path.is_file():
        fdr = pd.read_csv(fdr_path)
        sig = set(fdr.loc[fdr["sig_fdr_binomial"].fillna(False), "site_id"])
        robust: set[str] = set()
        if clust_path.is_file():
            clust = pd.read_csv(clust_path)
            all_radii = clust.groupby("site_id")["sig_05"].all()
            robust = set(all_radii[all_radii].index) & sig
        ctx["tier"] = ctx["site_id"].apply(
            lambda s: "robust" if s in robust else ("fdr_sig" if s in sig else "null")
        )
        ctx["enrichment"] = ctx["site_id"].map(
            fdr.set_index("site_id")["enrichment"]
        )
    return ctx


def main() -> None:
    site_ids = list_site_ids(SITES_DIR)
    bboxes = {
        sid: tuple(load_site_by_id(sid, SITES_DIR).bbox_wgs84) for sid in site_ids
    }

    print(
        f"Querying CMR for {len(bboxes)} sites x {len(COLLECTIONS)} collections "
        f"(cached in {CACHE_DIR.relative_to(ROOT)})\n"
    )
    coverage = build_coverage_table(bboxes, cache_dir=CACHE_DIR)

    ctx = load_site_context()
    coverage = coverage.merge(ctx, on="site_id", how="left")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    long_path = RESULTS_DIR / "gemx_coverage.csv"
    coverage.to_csv(long_path, index=False)

    by_site = summarize_by_site(coverage).merge(ctx, on="site_id", how="left")
    by_site = by_site.sort_values("swir_coverage_frac", ascending=False)
    site_path = RESULTS_DIR / "gemx_coverage_by_site.csv"
    by_site.to_csv(site_path, index=False)

    _print_summary(by_site)
    print(f"\nWrote {long_path.relative_to(ROOT)}")
    print(f"Wrote {site_path.relative_to(ROOT)}")


def _print_summary(by_site: pd.DataFrame) -> None:
    n = len(by_site)
    print("\n" + "=" * 78)
    print("  GEMx / AVIRIS coverage of the 109-site survey")
    print("=" * 78)

    print("\nSWIR-capable coverage (AVIRIS-C / -3 / -5, resolves 2.2 um Al-OH):")
    tiers = by_site["swir_tier"].value_counts().reindex(
        ["good", "partial", "marginal", "none"]
    ).fillna(0).astype(int)
    for tier, count in tiers.items():
        print(f"  {tier:<10} {count:>4} sites  ({count / n * 100:>5.1f}%)")

    print("\nBy significance tier — sites with any SWIR coverage:")
    if "tier" in by_site.columns:
        grp = by_site.groupby("tier")["has_swir"].agg(["sum", "count"])
        for tier, row in grp.iterrows():
            print(f"  {tier:<10} {int(row['sum']):>3} of {int(row['count']):>3}")

    if "category" in by_site.columns:
        print("\nBy deposit system — sites with any SWIR coverage:")
        grp = (
            by_site.groupby("category")["has_swir"]
            .agg(["sum", "count"])
            .sort_values("sum", ascending=False)
        )
        for cat, row in grp.iterrows():
            print(f"  {str(cat):<30} {int(row['sum']):>3} of {int(row['count']):>3}")

    print("\nThe SWIR-gap test cohort — epithermal + porphyry nulls with coverage,")
    print("ranked by usable coverage (these are the sites that decide whether")
    print("missing SWIR explains the below-chance medians):")
    if "category" in by_site.columns:
        cohort = by_site[
            by_site["category"].isin(["Epithermal", "Porphyry Cu-Mo-Au"])
            & by_site["has_swir"]
        ].sort_values("swir_coverage_frac", ascending=False)
        print(
            f"\n  {'site_id':<26}{'state':<7}{'category':<20}"
            f"{'cover':>7}{'sensor':>16}{'granules':>10}"
        )
        print("  " + "-" * 84)
        for _, r in cohort.iterrows():
            print(
                f"  {r['site_id']:<26}{str(r['state']):<7}{str(r['category']):<20}"
                f"{r['swir_coverage_frac'] * 100:>6.1f}%{str(r['swir_best_sensor']):>16}"
                f"{int(r['swir_granules']):>10}"
            )
        print(f"\n  {len(cohort)} of the 28 epithermal/porphyry sites are testable.")

    print("\nRobust sites (the positive result) — coverage for alteration-specificity:")
    if "tier" in by_site.columns:
        rob = by_site[by_site["tier"] == "robust"].sort_values(
            "swir_coverage_frac", ascending=False
        )
        for _, r in rob.iterrows():
            marker = "" if r["has_swir"] else "   <- no SWIR coverage"
            print(
                f"  {r['site_id']:<26}{str(r['state']):<5}"
                f"{r['swir_coverage_frac'] * 100:>6.1f}%  "
                f"MASTER {r['master_coverage_frac'] * 100:>5.1f}%{marker}"
            )


if __name__ == "__main__":
    main()
