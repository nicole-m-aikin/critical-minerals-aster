#!/usr/bin/env python
"""
Phase 5 figure — significance stability across MRDS clustering radii.

Reads results/phase5_clustering_sensitivity.csv and plots p_binomial (log
scale) vs. clustering radius for a contrast set drawn from the current
FDR-significant sites: the three declustering-robust sites with the lowest
worst-case p, and the two non-robust sites with the lowest raw p (the ones
whose significance most visibly collapses once MRDS point clusters are
deduplicated).

Usage:
    conda run -n aster-minerals python scripts/phase5_figure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from critical_minerals_aster.significance import significant_site_sets

REPO_ROOT = Path(__file__).parent.parent

# Okabe-Ito, colorblind-safe: blues/green for robust sites, warm hues for collapsing ones.
ROBUST_COLORS = ["#0072B2", "#009E73", "#56B4E9"]
COLLAPSE_COLORS = ["#D55E00", "#CC79A7"]


def pick_sites(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    fdr_sites, robust = significant_site_sets(REPO_ROOT)
    piv = df.pivot_table(index="site_id", columns="radius_m", values="p_binomial")
    piv = piv[piv.index.isin(fdr_sites)]
    worst = piv.max(axis=1)
    raw = piv[piv.columns.min()]
    robust_ids = worst[worst.index.isin(robust)].sort_values().index[: len(ROBUST_COLORS)].tolist()
    collapse_ids = raw[~raw.index.isin(robust)].sort_values().index[: len(COLLAPSE_COLORS)].tolist()
    return robust_ids, collapse_ids


def main() -> None:
    df = pd.read_csv(REPO_ROOT / "results" / "phase5_clustering_sensitivity.csv")
    robust_ids, collapse_ids = pick_sites(df)
    names = df.drop_duplicates("site_id").set_index("site_id").get("site_name")

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    series = [(s, c, "robust") for s, c in zip(robust_ids, ROBUST_COLORS)]
    series += [(s, c, "not robust") for s, c in zip(collapse_ids, COLLAPSE_COLORS)]
    for site_id, colour, tag in series:
        sub = df[df["site_id"] == site_id].sort_values("radius_m")
        # Floor at 1e-6 for log display; p_binomial can underflow to 0.0 at radius 0.
        y = sub["p_binomial"].clip(lower=1e-6)
        name = names.get(site_id, site_id) if names is not None else site_id
        ax.plot(
            sub["radius_m"], y,
            marker="o", markersize=6, linewidth=2,
            color=colour, label=f"{name} ({tag})",
        )

    ax.axhline(0.05, color="#555555", linewidth=1, linestyle="--", zorder=0)
    ax.text(1000, 0.05, "  α = 0.05", va="bottom", ha="right", fontsize=9, color="#555555")

    ax.set_yscale("log")
    ax.set_xticks([0, 250, 500, 1000])
    ax.set_xlabel("DBSCAN clustering radius (m)")
    ax.set_ylabel("p-value (bedrock/cover-stratified null, log scale)")
    ax.set_title(
        "Significance stability under MRDS spatial-clustering correction\n"
        "(medoid deduplication, min_samples=2)",
        fontsize=11,
    )
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", which="major", color="#e5e5e5", linewidth=0.8, zorder=-1)

    fig.tight_layout()
    out_path = REPO_ROOT / "figures" / "phase5_clustering_stability.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    print(f"Saved {out_path}: robust {robust_ids}, not robust {collapse_ids}")


if __name__ == "__main__":
    main()
