#!/usr/bin/env python
"""
Figure 4 -- site-level observed hit rate vs. site-specific null, all 45 sites.

Scientific question: across the full survey, how does each site's observed
critical-mineral hit rate compare to what its own zone geometry predicts by
chance? This is the single figure that most directly visualizes what Phase
2 changed (replacing one pooled null with 45 site-specific ones) and what
Phase 4 found (which of the resulting enrichments survive multiple-testing
correction).

Form: scatter, x = site-specific null (%), y = observed hit rate (%), plus
a 1:1 reference line (y=x, "no enrichment"). Points above the line are
enriched; below are anti-correlated. Color = FDR significance tier
(binomial FDR-significant / continuous-only FDR-significant / not
significant). Marker size ~ n_crit (sample size), so a reader can see at a
glance which points carry more statistical weight.

Usage:
    conda run -n aster-minerals python scripts/fig04_enrichment_vs_null.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent


def main() -> None:
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    fdr = pd.read_csv(REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv")
    clustering = pd.read_csv(REPO_ROOT / "results" / "phase5_clustering_sensitivity.csv")
    clustering_1000 = clustering[clustering["radius_m"] == 1000].set_index("site_id")["sig_05"]

    df = ss.merge(fdr[["site_id", "sig_fdr_binomial", "sig_fdr_mannwhitney"]], on="site_id", how="left")
    df["clustering_robust"] = df["site_id"].map(clustering_1000).fillna(False) & df["sig_fdr_binomial"].fillna(False)

    def tier(row):
        if row["clustering_robust"]:
            return "Binomial FDR-sig. + clustering-robust"
        if row["sig_fdr_binomial"]:
            return "Binomial FDR-significant"
        if row["sig_fdr_mannwhitney"]:
            return "Continuous-score FDR-significant only"
        return "Not FDR-significant"

    df["tier"] = df.apply(tier, axis=1)

    colors = {
        "Binomial FDR-sig. + clustering-robust": "#D55E00",
        "Binomial FDR-significant": "#E69F00",
        "Continuous-score FDR-significant only": "#0072B2",
        "Not FDR-significant": "#BBBBBB",
    }
    order = list(colors)

    fig, ax = plt.subplots(figsize=(7.5, 7))
    # x-axis is the site-specific null p0 (zone area / footprint area); it never
    # exceeds ~18% for any site, so cap the view at 25% rather than letting the
    # (much larger) hit-rate range stretch it.
    x_max = 25.0
    y_max = df["hit_rate_crit_pct"].max() * 1.08
    ax.plot([0, x_max], [0, x_max], color="#555555", linewidth=1.2, linestyle="--", zorder=1,
            label="No enrichment (y = x)")

    size_scale = 900 / df["n_crit"].max()
    for t in order:
        sub = df[df["tier"] == t]
        if sub.empty:
            continue
        ax.scatter(sub["p0_null"] * 100, sub["hit_rate_crit_pct"],
                   s=np.clip(sub["n_crit"] * size_scale, 25, None),
                   color=colors[t], alpha=0.85, edgecolors="white", linewidths=0.6,
                   label=t, zorder=3)

    for _, r in df[df["tier"] != "Not FDR-significant"].iterrows():
        ax.annotate(r["site_id"], (r["p0_null"] * 100, r["hit_rate_crit_pct"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=8)

    ax.set_xlim(0, x_max)
    ax.set_ylim(-2, y_max)
    ax.set_xlabel("Site-specific null p₀ (%) — zone area / footprint area")
    ax.set_ylabel("Observed critical-mineral hit rate (%)")
    ax.set_title(
        f"Figure 4 — observed hit rate vs. site-specific null, all {len(df)} sites\n"
        "Point size ∝ n (critical-mineral deposits); above the line = enriched",
        fontsize=11,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color="#eeeeee", linewidth=0.7, zorder=0)
    ax.legend(loc="upper left", frameon=False, fontsize=8.5)

    fig.tight_layout()
    out_path = REPO_ROOT / "figures" / "fig04_enrichment_vs_null.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
