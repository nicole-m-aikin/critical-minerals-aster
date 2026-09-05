#!/usr/bin/env python
"""
Figure 1 -- map of all study sites over a US state-boundary line map.

Scientific question: where are the sites, and which ones carry the
final corrected significance result?

Encoding (all three dimensions carried by the marker itself -- no
overplotted rings):
  * fill colour  = deposit-system category (Phase 6)
  * heavy black outline = FDR-significant (Phase 4, binomial or continuous)
  * marker shape = star -> also robust to spatial-clustering correction
                   (Phase 5, sig at every DBSCAN radius); circle -> not

The state outlines come from a small committed GeoJSON
(data/basemap/us_states.geojson, ~90 KB) so the figure stays
offline-reproducible -- no tile fetch.

Usage:
    conda run -n aster-minerals python scripts/fig01_site_map.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from matplotlib.lines import Line2D

from critical_minerals_aster.config import load_site_by_id

REPO_ROOT = Path(__file__).parent.parent

CATEGORY_COLORS = {
    "Porphyry Cu-Mo-Au": "#0072B2",
    "Skarn/Carbonate Replacement": "#D55E00",
    "Epithermal": "#E69F00",
    "Sediment-hosted": "#009E73",
    "Alkaline/Carbonatite": "#CC79A7",
    "Uranium/Energy": "#56B4E9",
    "Mafic/Ultramafic": "#8B4513",
    "Volcanogenic/VMS": "#999999",
    "Mixed/Other": "#666666",
}


def main() -> None:
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        site_ids = yaml.safe_load(f)["sites"]

    cats = pd.read_csv(REPO_ROOT / "results" / "deposit_system_categories.csv")
    fdr = pd.read_csv(REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv")
    clustering = pd.read_csv(REPO_ROOT / "results" / "phase5_clustering_sensitivity.csv")
    piv = clustering.pivot_table(index="site_id", columns="radius_m", values="sig_05")
    robust_all_radii = piv[(piv == True).all(axis=1)].index  # noqa: E712

    rows = []
    for site_id in site_ids:
        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        lon = (site.bbox_wgs84[0] + site.bbox_wgs84[2]) / 2
        lat = (site.bbox_wgs84[1] + site.bbox_wgs84[3]) / 2
        rows.append({"site_id": site_id, "lon": lon, "lat": lat, "name": site.name})
    df = pd.DataFrame(rows).merge(cats[["site_id", "category"]], on="site_id", how="left")
    df = df.merge(fdr[["site_id", "sig_fdr_binomial", "sig_fdr_mannwhitney"]], on="site_id", how="left")
    df["sig_fdr_either"] = df["sig_fdr_binomial"].fillna(False) | df["sig_fdr_mannwhitney"].fillna(False)
    df["clustering_robust"] = (
        df["site_id"].isin(robust_all_radii) & df["sig_fdr_binomial"].fillna(False)
    )
    df["color"] = df["category"].map(CATEGORY_COLORS).fillna("#333333")

    fig, ax = plt.subplots(figsize=(14, 8))

    # --- US state-boundary line map ---
    states = gpd.read_file(REPO_ROOT / "data" / "basemap" / "us_states.geojson")
    states.boundary.plot(ax=ax, color="#bcbcbc", linewidth=0.6, zorder=1)

    # --- sites: shape = clustering-robustness, outline = FDR, fill = category ---
    for is_robust, marker, size in [(False, "o", 70), (True, "*", 430)]:
        for is_sig, edge_c, edge_w in [(False, "#b0b0b0", 0.5), (True, "black", 1.8)]:
            sub = df[(df["clustering_robust"] == is_robust) & (df["sig_fdr_either"] == is_sig)]
            if sub.empty:
                continue
            ax.scatter(sub["lon"], sub["lat"], s=size, marker=marker,
                       c=sub["color"], edgecolors=edge_c, linewidths=edge_w, zorder=4)

    # Label the clustering-robust (star) sites only; the SE-Arizona skarn belt
    # (bisbee/tombstone/courtland_gleeson/christmas) is a ~1-degree cluster -> label once.
    az_cluster = {"bisbee", "tombstone", "courtland_gleeson", "christmas"}
    for _, r in df[df["clustering_robust"]].iterrows():
        if r["site_id"] in az_cluster:
            continue
        ax.annotate(r["site_id"], (r["lon"], r["lat"]), textcoords="offset points",
                    xytext=(6, 4), fontsize=7.5, color="#222222", zorder=6)
    az_pts = df[df["site_id"].isin(az_cluster) & df["clustering_robust"]]
    if not az_pts.empty:
        ax.annotate("SE-AZ skarn belt (4):\nbisbee, tombstone,\ncourtland_gleeson,\nchristmas",
                    (az_pts["lon"].min(), az_pts["lat"].mean()), textcoords="offset points",
                    xytext=(-96, -30), fontsize=7, color="#222222",
                    arrowprops=dict(arrowstyle="-", color="#888888", linewidth=0.8), zorder=6)

    ax.set_xlim(-125.5, -66)
    ax.set_ylim(24, 49.5)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect(1.30)  # ~1/cos(39 N)

    n_total = len(df)
    n_sig = int(df["sig_fdr_either"].sum())
    n_rob = int(df["clustering_robust"].sum())
    ax.set_title(
        f"Figure 1 — {n_total} critical-mineral study sites\n"
        f"fill = deposit system;  heavy outline = FDR-significant (n={n_sig});  "
        f"star = also spatial-clustering-robust (n={n_rob})",
        fontsize=11,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    cat_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=c,
               markeredgecolor="#b0b0b0", label=cat.replace("Carbonate Replacement", "Carb. Repl."))
        for cat, c in CATEGORY_COLORS.items()
    ]
    shape_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor="#cccccc",
               markeredgecolor="black", markeredgewidth=1.8, label="FDR-significant"),
        Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor="#cccccc",
               markeredgecolor="#b0b0b0", markeredgewidth=0.5, label="not FDR-significant"),
        Line2D([0], [0], marker="*", linestyle="", markersize=15, markerfacecolor="#cccccc",
               markeredgecolor="black", markeredgewidth=1.2, label="spatial-clustering-robust"),
    ]
    leg1 = ax.legend(handles=cat_handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                     frameon=False, fontsize=8, title="Deposit system", title_fontsize=9)
    ax.add_artist(leg1)
    ax.legend(handles=shape_handles, loc="upper left", bbox_to_anchor=(1.01, 0.42),
              frameon=False, fontsize=8, title="Significance", title_fontsize=9)

    fig.tight_layout()
    out_path = REPO_ROOT / "figures" / "fig01_site_map.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
