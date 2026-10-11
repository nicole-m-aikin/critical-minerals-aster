#!/usr/bin/env python
"""
Phase 4 — Benjamini-Hochberg FDR correction across the 45 primary site tests.

Input p-values: results/site_specific_null_significance.csv's p_binomial
column — the primary test (site-specific-null one-sided exact binomial,
Phase 2). Also applies the same correction to the continuous-score
Mann-Whitney U p-values from Phase 3 (results/
phase3_monte_carlo_and_continuous_score.csv), since Phase 3 showed several
sites carry a real continuous-score signal independent of the binary test
and that claim needs the same multiple-testing discipline before it can be
reported.

Method: statsmodels.stats.multitest.multipletests(pvals, alpha=0.05,
method="fdr_bh") -- Benjamini-Hochberg. See docs/results.md for the full explanation of what a q-value means and why BH (not
Bonferroni) is the right choice for a 45-site survey.

Since the bedrock-vs-cover stratified null became primary, p_binomial
(Phase 2) and p_mannwhitney (Phase 3) are the STRATIFIED p-values. The
previous geometric-null p-values are corrected the same way, as a separate
family, into the *_geometric columns for comparison only.

Outputs
-------
    results/phase4_fdr_corrected_significance.csv
        site_id, p_binomial, q_binomial, sig_fdr_binomial,
        p_mannwhitney, q_mannwhitney, sig_fdr_mannwhitney,
        p_binomial_geometric, q_binomial_geometric, sig_fdr_binomial_geometric,
        p_mannwhitney_geometric, q_mannwhitney_geometric, sig_fdr_mannwhitney_geometric

Usage:
    conda run -n aster-minerals python scripts/phase4_fdr_correction.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).parent.parent
ALPHA = 0.05


def bh(df: pd.DataFrame, p_col: str, q_col: str, sig_col: str) -> pd.Series:
    """BH-FDR over the non-missing p-values of ``p_col``; returns the validity mask.

    Sites without a p-value (e.g. MWU with < 3 scored deposit pixels) are
    excluded from the family rather than imputed — they were never valid tests.
    """
    mask = df[p_col].notna()
    df[q_col] = np.nan
    df[sig_col] = False
    if mask.sum() > 0:
        reject, q, _, _ = multipletests(df.loc[mask, p_col].values, alpha=ALPHA, method="fdr_bh")
        df.loc[mask, q_col] = np.round(q, 4)
        df.loc[mask, sig_col] = reject
    df[sig_col] = df[sig_col].astype(bool)
    return mask


def main() -> None:
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    p3 = pd.read_csv(REPO_ROOT / "results" / "phase3_monte_carlo_and_continuous_score.csv")

    df = ss[["site_id", "site_name", "n_crit", "hits_crit", "hit_rate_crit_pct",
             "p0_null", "enrichment", "p_binomial",
             "p0_geometric", "enrichment_geometric", "p_binomial_geometric"]].merge(
        p3[["site_id", "p_monte_carlo", "p_mannwhitney", "auc_continuous_score",
            "p_monte_carlo_geometric", "p_mannwhitney_geometric", "auc_continuous_score_geometric"]],
        on="site_id", how="left",
    )

    # --- Primary correction: stratified binomial (Poisson-binomial) p-values ---
    bh(df, "p_binomial", "q_binomial", "sig_fdr_binomial")
    # --- Secondary correction: stratified continuous-score Mann-Whitney p-values ---
    mwu_mask = bh(df, "p_mannwhitney", "q_mannwhitney", "sig_fdr_mannwhitney")
    # --- Comparison only: the same corrections under the geometric null ---
    bh(df, "p_binomial_geometric", "q_binomial_geometric", "sig_fdr_binomial_geometric")
    bh(df, "p_mannwhitney_geometric", "q_mannwhitney_geometric", "sig_fdr_mannwhitney_geometric")

    df = df.sort_values("p_binomial")
    out_path = REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv"
    df.to_csv(out_path, index=False)

    n_raw_binom = int((df["p_binomial"] < ALPHA).sum())
    n_fdr_binom = int(df["sig_fdr_binomial"].sum())
    n_raw_mwu = int((df.loc[mwu_mask, "p_mannwhitney"] < ALPHA).sum())
    n_fdr_mwu = int(df["sig_fdr_mannwhitney"].sum())

    print(f"Saved {out_path}\n")
    n_sites = len(df)
    print(f"BINOMIAL (bedrock/cover-stratified null), {n_sites} tests, alpha={ALPHA}, BH-FDR:")
    print(f"  raw p < {ALPHA}: {n_raw_binom}/{n_sites}")
    print(f"  FDR-significant (q < {ALPHA}): {n_fdr_binom}/{n_sites}\n")

    print(f"CONTINUOUS-SCORE MANN-WHITNEY U, {int(mwu_mask.sum())} valid tests, alpha={ALPHA}, BH-FDR:")
    print(f"  raw p < {ALPHA}: {n_raw_mwu}/{int(mwu_mask.sum())}")
    print(f"  FDR-significant (q < {ALPHA}): {n_fdr_mwu}/{int(mwu_mask.sum())}\n")

    cols = ["site_id", "n_crit", "hit_rate_crit_pct", "enrichment",
            "p_binomial", "q_binomial", "sig_fdr_binomial",
            "p_mannwhitney", "q_mannwhitney", "sig_fdr_mannwhitney"]
    print(df[cols].to_string(index=False))

    print("\n=== Sites significant on EITHER corrected test ===")
    either = df[df["sig_fdr_binomial"] | df["sig_fdr_mannwhitney"]]
    print(either[["site_id", "sig_fdr_binomial", "sig_fdr_mannwhitney"]].to_string(index=False))

    print("\n=== FDR-significant counts: geometric vs stratified null ===")
    for label, b, m in (("geometric", "sig_fdr_binomial_geometric", "sig_fdr_mannwhitney_geometric"),
                        ("stratified (primary)", "sig_fdr_binomial", "sig_fdr_mannwhitney")):
        print(f"  {label:22s} binomial {int(df[b].sum()):3d}   MWU {int(df[m].sum()):3d}   "
              f"either {int((df[b] | df[m]).sum()):3d}   (of {n_sites})")
    lost = sorted(df.loc[df["sig_fdr_binomial_geometric"] & ~df["sig_fdr_binomial"], "site_id"])
    gained = sorted(df.loc[~df["sig_fdr_binomial_geometric"] & df["sig_fdr_binomial"], "site_id"])
    print(f"  binomial: lost under stratification {lost}")
    print(f"  binomial: gained under stratification {gained}")


if __name__ == "__main__":
    main()
