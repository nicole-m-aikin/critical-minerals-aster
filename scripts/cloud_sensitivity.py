#!/usr/bin/env python
"""
Cloud-contamination sensitivity of the headline significance results.

results/granule_cloud_audit.csv (scripts/audit_granule_cloud_cover.py) shows
most site mosaics blend cloudy ASTER scenes: the original selection ranked
granules on coverage and band count only. Cloud tops are cold, so they depress
TIR radiance and distort band ratios, and a percentile classifier can turn a
cloud edge into a "strong anomaly". Re-downloading clear-sky mosaics for 109
sites is the real fix; this script asks, with the data already on disk,
whether cloud could be driving the reported signal:

1. Association — Spearman correlation of site blended cloud % (and night-scene
   fraction) with enrichment across all sites; cloud % of FDR-significant vs
   non-significant sites (Mann–Whitney).
2. Clear-subset headline — BH-FDR re-applied within sites whose blended bbox
   cloud is below each cutoff (10 / 30 %), and the robust count within them.

Outputs
-------
    results/cloud_sensitivity_sites.csv     per-site cloud, night fraction, results
    results/cloud_sensitivity_summary.csv   association tests + clear-subset counts

Usage
-----
    conda run -n aster-minerals python scripts/cloud_sensitivity.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS = REPO_ROOT / "results"
ALPHA = 0.05
CUTOFFS = (10.0, 30.0)


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(ranked, 1.0)
    return q


def site_table() -> pd.DataFrame:
    audit = pd.read_csv(RESULTS / "granule_cloud_audit.csv")
    per_site = audit.groupby("site_id").agg(
        blend_cloud_pct=("site_blend_cloud_pct", "first"),
        site_flag=("site_flag", "first"),
        n_granules=("granule_id", "nunique"),
        night_fraction=("day_night", lambda s: float((s.str.lower() == "night").mean())),
    )
    fdr = pd.read_csv(RESULTS / "phase4_fdr_corrected_significance.csv").set_index("site_id")
    p5 = pd.read_csv(RESULTS / "phase5_clustering_sensitivity.csv")
    robust_all = p5.pivot(index="site_id", columns="radius_m", values="sig_05").all(axis=1)
    df = fdr[["n_crit", "enrichment", "p_binomial", "sig_fdr_binomial", "sig_fdr_mannwhitney"]].join(per_site)
    df["robust"] = df["sig_fdr_binomial"] & robust_all.reindex(df.index).fillna(False).astype(bool)
    return df


def main() -> None:
    df = site_table()
    missing = df["blend_cloud_pct"].isna().sum()
    df.to_csv(RESULTS / "cloud_sensitivity_sites.csv")

    rows: list[dict] = []
    ok = df.dropna(subset=["blend_cloud_pct", "enrichment"])
    rho, p_rho = spearmanr(ok["blend_cloud_pct"], ok["enrichment"])
    rows.append({"test": "spearman_cloud_vs_enrichment", "n": len(ok), "stat": rho, "p": p_rho})
    rho_n, p_n = spearmanr(ok["night_fraction"], ok["enrichment"])
    rows.append({"test": "spearman_night_fraction_vs_enrichment", "n": len(ok), "stat": rho_n, "p": p_n})
    sig = ok[ok["sig_fdr_binomial"]]["blend_cloud_pct"]
    non = ok[~ok["sig_fdr_binomial"]]["blend_cloud_pct"]
    u, p_u = mannwhitneyu(sig, non, alternative="two-sided")
    rows.append({
        "test": "cloud_pct_fdr_sig_vs_not (median sig / not)",
        "n": len(ok), "stat": float(u), "p": p_u,
        "note": f"{sig.median():.1f} / {non.median():.1f}",
    })

    print(f"Sites with cloud audit: {len(ok)} (missing {missing})")
    print(f"Spearman cloud% vs enrichment: rho={rho:.3f} p={p_rho:.3g}")
    print(f"Spearman night-fraction vs enrichment: rho={rho_n:.3f} p={p_n:.3g}")
    print(f"Median cloud%: FDR-sig {sig.median():.1f} vs not {non.median():.1f} (MWU p={p_u:.3g})")

    for cut in CUTOFFS:
        sub = ok[ok["blend_cloud_pct"] < cut].copy()
        sub["q_sub"] = bh_qvalues(sub["p_binomial"].to_numpy())
        sig_sub = sub[sub["q_sub"] < ALPHA]
        robust_sub = sorted(set(sig_sub.index) & set(ok.index[ok["robust"]]))
        rows.append({
            "test": f"clear_subset_blend_cloud_lt_{int(cut)}pct",
            "n": len(sub),
            "stat": len(sig_sub),
            "p": np.nan,
            "note": f"FDR-sig within subset: {sorted(sig_sub.index)}; robust: {robust_sub}",
        })
        print(f"\nBlended cloud < {cut:.0f}%: {len(sub)} sites, {len(sig_sub)} FDR-significant "
              f"within subset, {len(robust_sub)} of the robust set")
        print(f"  FDR-sig: {sorted(sig_sub.index)}")

    pd.DataFrame(rows).to_csv(RESULTS / "cloud_sensitivity_summary.csv", index=False)
    robust = ok[ok["robust"]].sort_values("blend_cloud_pct")
    print("\nRobust sites by blended cloud:")
    print(robust[["blend_cloud_pct", "night_fraction", "n_granules", "enrichment"]].round(2).to_string())


if __name__ == "__main__":
    main()
