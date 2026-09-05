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

Outputs
-------
    results/phase4_fdr_corrected_significance.csv
        site_id, p_binomial, q_binomial, sig_fdr_binomial,
        p_mannwhitney, q_mannwhitney, sig_fdr_mannwhitney

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


def main() -> None:
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    p3 = pd.read_csv(REPO_ROOT / "results" / "phase3_monte_carlo_and_continuous_score.csv")

    df = ss[["site_id", "site_name", "n_crit", "hits_crit", "hit_rate_crit_pct",
             "p0_null", "enrichment", "p_binomial"]].merge(
        p3[["site_id", "p_monte_carlo", "p_mannwhitney", "auc_continuous_score"]],
        on="site_id", how="left",
    )

    # --- Primary correction: binomial p-values (site-specific null) ---
    reject_b, q_b, _, _ = multipletests(df["p_binomial"].values, alpha=ALPHA, method="fdr_bh")
    df["q_binomial"] = np.round(q_b, 4)
    df["sig_fdr_binomial"] = reject_b

    # --- Secondary correction: continuous-score Mann-Whitney p-values ---
    # Two sites (elk_creek, green_river) have no MWU p-value (n_crit too
    # small / zero) -- exclude them from this correction's test family
    # rather than imputing, since they were never valid tests to begin with.
    mwu_mask = df["p_mannwhitney"].notna()
    df["q_mannwhitney"] = np.nan
    df["sig_fdr_mannwhitney"] = False
    if mwu_mask.sum() > 0:
        reject_m, q_m, _, _ = multipletests(
            df.loc[mwu_mask, "p_mannwhitney"].values, alpha=ALPHA, method="fdr_bh"
        )
        df.loc[mwu_mask, "q_mannwhitney"] = np.round(q_m, 4)
        df.loc[mwu_mask, "sig_fdr_mannwhitney"] = reject_m

    df = df.sort_values("p_binomial")
    out_path = REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv"
    df.to_csv(out_path, index=False)

    n_raw_binom = int((df["p_binomial"] < ALPHA).sum())
    n_fdr_binom = int(df["sig_fdr_binomial"].sum())
    n_raw_mwu = int((df.loc[mwu_mask, "p_mannwhitney"] < ALPHA).sum())
    n_fdr_mwu = int(df["sig_fdr_mannwhitney"].sum())

    print(f"Saved {out_path}\n")
    n_sites = len(df)
    print(f"BINOMIAL (site-specific null), {n_sites} tests, alpha={ALPHA}, BH-FDR:")
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


if __name__ == "__main__":
    main()
