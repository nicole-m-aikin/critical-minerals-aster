#!/usr/bin/env python
"""
Figure 6 -- robustness and discovery-bias summary.

Two panels answering two different questions:

Left: how many sites survive each successive correction, in order applied?
This is the single figure that summarizes the entire corrected pipeline's
net effect on the headline count -- pooled null (original) -> site-specific
null (Phase 2) -> BH-FDR (Phase 4, binomial primary test) -> spatial-
clustering-robust (Phase 5, the strictest cumulative bar). Each bar is a
strict subset of the one before it by construction (this is not 4
independent counts -- it is one shrinking set), which the connecting
funnel shape is meant to communicate visually.

Right: the corrected (Phase 7) discovery-era result, unconditional across
all 45 sites -- observed vs. null hit rate for pre-/post-1950 cohorts, the
result that overturned the original paper's discovery-bias claim.

Usage:
    conda run -n aster-minerals python scripts/fig06_robustness_summary.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent


def _funnel_counts() -> tuple[list[int], int]:
    """Derive the four correction-stage counts from the result CSVs (no hardcoding)."""
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    p4 = pd.read_csv(REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv")
    p5 = pd.read_csv(REPO_ROOT / "results" / "phase5_clustering_sensitivity.csv")
    n_total = len(ss)
    n_pooled = int((ss["p_binomial_pooled"] < 0.05).sum())
    n_raw = int((ss["p_binomial"] < 0.05).sum())
    fdr_sites = set(p4.loc[p4["sig_fdr_binomial"], "site_id"])
    piv = p5.pivot_table(index="site_id", columns="radius_m", values="sig_05")
    robust_all = set(piv[(piv == True).all(axis=1)].index)  # noqa: E712
    n_final = len(fdr_sites & robust_all)
    return [n_pooled, n_raw, len(fdr_sites), n_final], n_total


def main() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5))
    counts, n_total = _funnel_counts()

    # --- Left panel: correction funnel ---
    ax = axes[0]
    labels = [
        "Pooled null\n(raw p<0.05)",
        "Site-specific\nnull (Phase 2)",
        "+ BH-FDR\n(Phase 4, primary)",
        "+ spatial-clustering\nrobust (Phase 5)",
    ]
    colors = ["#BBBBBB", "#E69F00", "#D55E00", "#8B0000"]

    bars = ax.bar(labels, counts, color=colors, width=0.6, zorder=3)
    for bar, c in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, c + 0.3, str(c),
                ha="center", fontsize=13, fontweight="bold")
    ax.set_ylabel("Number of sites called significant")
    ax.set_ylim(0, max(counts) + 6)
    ax.set_title(f"Sites significant, at each correction stage ({n_total}-site survey)\n"
                 "(FDR ⊂ site-specific-raw; clustering-robust ⊂ FDR. Pooled vs. site-specific "
                 "are different sets, not nested)", fontsize=9.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#eeeeee", linewidth=0.7, zorder=0)
    ax.tick_params(axis="x", labelsize=8.5)

    # --- Right panel: discovery-bias, corrected, unconditional ---
    ax = axes[1]
    pooled = pd.read_csv(REPO_ROOT / "results" / "phase7_discovery_bias_pooled.csv")
    primary = pooled[pooled["framing"].str.contains("PRIMARY")]

    x = np.arange(2)
    width = 0.32
    observed = primary["hr_pct"].values
    null = primary["null_hr_pct"].values
    cohort_labels = primary["cohort"].values

    ax.bar(x - width / 2, observed, width, label="Observed hit rate", color="#0072B2", zorder=3)
    ax.bar(x + width / 2, null, width, label="Site-specific null (pooled)", color="#999999", zorder=3)
    for i, row in enumerate(primary.itertuples()):
        ax.text(i, max(row.hr_pct, row.null_hr_pct) + 1.5,
                f"{row.enrichment:.2f}×\np={row.pooled_p:.2f}", ha="center", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(cohort_labels)
    ax.set_ylabel("Hit rate (%)")
    ax.set_ylim(0, max(list(observed) + list(null)) + 6)
    ax.set_title(
        f"Discovery-era split, corrected (Phase 7)\nunconditional across {int(primary['n_sites'].iloc[0])} dated-record sites — no signal in either cohort",
        fontsize=10.5,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#eeeeee", linewidth=0.7, zorder=0)
    ax.legend(loc="upper right", frameon=False, fontsize=8.5)

    fig.suptitle("Figure 6 — robustness and discovery-bias summary", fontsize=13, y=1.02)
    fig.tight_layout()
    out_path = REPO_ROOT / "figures" / "fig06_robustness_summary.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
