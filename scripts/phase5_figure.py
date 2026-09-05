#!/usr/bin/env python
"""
Phase 5 figure — significance stability across MRDS clustering radii.

Reads results/phase5_clustering_sensitivity.csv and plots p_binomial (log
scale) vs. clustering radius for the 5 sites that were significant at
radius=0 (raw, no deduplication) at p<0.05, so a reader can see directly
whether spatial-clustering correction changes which sites the paper can
defend.

Usage:
    conda run -n aster-minerals python scripts/phase5_figure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent

# Colorblind-safe categorical set (Okabe-Ito), assigned by fixed site identity
# -- never re-cycled if the site list changes.
# Post-109-site expansion: contrast the declustering-ROBUST set (bisbee + two
# new NM skarn replication sites) against sites whose raw significance COLLAPSES
# once MRDS point clusters are deduplicated (eureka -- an original headline site;
# tri_state -- an MVT district, raw p~0 -> p~1 by 500 m).
SITE_COLORS = {
    "bisbee": "#0072B2",             # blue
    "magdalena_kelly": "#009E73",    # bluish green
    "organ_mountains": "#56B4E9",    # sky blue
    "eureka": "#D55E00",             # vermillion
    "tri_state": "#CC79A7",          # reddish purple
}
SITE_LABELS = {
    "bisbee": "Bisbee, AZ (robust)",
    "magdalena_kelly": "Magdalena-Kelly, NM (robust, new)",
    "organ_mountains": "Organ Mtns, NM (robust, new)",
    "eureka": "Eureka, NV (collapses)",
    "tri_state": "Tri-State MVT (collapses)",
}


def main() -> None:
    df = pd.read_csv(REPO_ROOT / "results" / "phase5_clustering_sensitivity.csv")
    df = df[df["site_id"].isin(SITE_COLORS)]

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    for site_id in SITE_COLORS:
        sub = df[df["site_id"] == site_id].sort_values("radius_m")
        # Floor at 1e-6 for log display; p_binomial can be exactly 0.0 for
        # Bisbee/Eureka at radius=0 (below float precision, not a real zero).
        y = sub["p_binomial"].clip(lower=1e-6)
        ax.plot(
            sub["radius_m"], y,
            marker="o", markersize=6, linewidth=2,
            color=SITE_COLORS[site_id], label=SITE_LABELS[site_id],
        )

    ax.axhline(0.05, color="#555555", linewidth=1, linestyle="--", zorder=0)
    ax.text(1000, 0.05, "  α = 0.05", va="bottom", ha="right", fontsize=9, color="#555555")

    ax.set_yscale("log")
    ax.set_xticks([0, 250, 500, 1000])
    ax.set_xlabel("DBSCAN clustering radius (m)")
    ax.set_ylabel("Binomial p-value (site-specific null, log scale)")
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
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
