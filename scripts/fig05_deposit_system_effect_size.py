#!/usr/bin/env python
"""
Figure 5 -- deposit-system effect-size summary.

Scientific question: which deposit-system categories, in aggregate, show
enrichment above chance, and how uncertain is that estimate given how few
sites populate most categories?

Form: horizontal dot-and-whisker plot, one row per category, sorted by
median enrichment. Dot = median enrichment; whiskers = 95% bootstrap CI
(Phase 6, method="percentile"); a vertical reference line at 1.0x (no
enrichment). Categories with insufficient n for a bootstrap CI (Uranium/
Energy, Volcanogenic/VMS, Mafic/Ultramafic -- all 1-2 sites) show a dot
with no whiskers and an explicit "n=1" / "n=2" label rather than a
misleadingly precise interval.

Usage:
    conda run -n aster-minerals python scripts/fig05_deposit_system_effect_size.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent


def main() -> None:
    df = pd.read_csv(REPO_ROOT / "results" / "phase6_deposit_system_categories.csv")
    df = df.sort_values("median_enrichment", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    y = range(len(df))

    has_ci = df["ci_note"] == "ok"
    xerr_low = (df["median_enrichment"] - df["bootstrap_ci_low"]).where(has_ci, 0)
    xerr_high = (df["bootstrap_ci_high"] - df["median_enrichment"]).where(has_ci, 0)

    ax.errorbar(
        df.loc[has_ci, "median_enrichment"], [i for i in y if has_ci[i]],
        xerr=[xerr_low[has_ci], xerr_high[has_ci]],
        fmt="o", markersize=9, color="#0072B2", ecolor="#0072B2",
        elinewidth=1.6, capsize=4, zorder=3,
    )
    ax.scatter(df.loc[~has_ci, "median_enrichment"], [i for i in y if not has_ci[i]],
               marker="o", s=80, color="#999999", zorder=3)

    ax.axvline(1.0, color="#555555", linewidth=1.2, linestyle="--", zorder=1)
    ax.annotate("no enrichment (1.0×)", xy=(1.0, 1.0), xycoords=("data", "axes fraction"),
                xytext=(4, -2), textcoords="offset points", fontsize=8, color="#555555",
                va="top", ha="left")

    labels = []
    for _, r in df.iterrows():
        n_lbl = f"n={r['n_sites']}"
        fdr_lbl = f"{r['n_fdr_significant']}/{r['n_sites']} FDR-sig."
        labels.append(f"{r['category']}  ({n_lbl}, {fdr_lbl})")
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)

    ax.set_xlabel("Median enrichment (observed / site-specific null hit rate)")
    ax.set_title(
        f"Figure 5 — deposit-system effect size, {int(df['n_sites'].sum())}-site survey\n"
        "Dot = median enrichment; whiskers = 95% bootstrap CI (percentile method, 10,000 resamples)",
        fontsize=10.5,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#eeeeee", linewidth=0.7, zorder=0)

    for i, ok in enumerate(~has_ci):
        if ok:
            ax.annotate("insufficient n\nfor bootstrap CI", (df.loc[i, "median_enrichment"], i),
                        textcoords="offset points", xytext=(10, -2), fontsize=7.5, color="#777777")

    fig.tight_layout()
    out_path = REPO_ROOT / "figures" / "fig05_deposit_system_effect_size.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
