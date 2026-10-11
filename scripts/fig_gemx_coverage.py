#!/usr/bin/env python
"""GEMx / AVIRIS coverage figures for the 109-site survey.

Scientific question: which study sites actually sit under an imaging-
spectrometer flight line, and which of those can test the SWIR detection
gap named in docs/results.md limitation 2?

Two figures, one script:

  figures/gemx_coverage_map.png
      Geographic. Fill = SWIR bbox-coverage tier; shape = scientific
      role (robust / epithermal-or-porphyry / other). Offline state
      boundaries, same basemap as Figure 1.

  figures/gemx_coverage_by_system.png
      Deposit-system summary: what fraction of each category sits
      under GEMx-era AVIRIS coverage.

Reads results/gemx_coverage_by_site.csv (from
scripts/gemx_coverage_report.py). No network access.

Usage:
    conda run -n aster-minerals python scripts/fig_gemx_coverage.py
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
from matplotlib.patches import Patch

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

TIER_COLORS = {
    "good": "#117733",
    "partial": "#DDCC77",
    "marginal": "#EE7733",
    "none": "#BBBBBB",
}

SWIR_GAP_CATEGORIES = {"Epithermal", "Porphyry Cu-Mo-Au"}
AZ_CLUSTER = {"bisbee", "tombstone", "courtland_gleeson", "christmas"}


def load_coverage() -> pd.DataFrame:
    path = REPO_ROOT / "results" / "gemx_coverage_by_site.csv"
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} is missing. Run scripts/gemx_coverage_report.py first."
        )
    cov = pd.read_csv(path)
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        site_ids = yaml.safe_load(f)["sites"]
    rows = []
    for site_id in site_ids:
        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        lon = (site.bbox_wgs84[0] + site.bbox_wgs84[2]) / 2
        lat = (site.bbox_wgs84[1] + site.bbox_wgs84[3]) / 2
        rows.append({"site_id": site_id, "lon": lon, "lat": lat})
    df = pd.DataFrame(rows).merge(cov, on="site_id", how="left")
    df["swir_coverage_frac"] = df["swir_coverage_frac"].fillna(0.0)
    df["swir_tier"] = df["swir_tier"].fillna("none")
    df["has_swir"] = df["has_swir"].fillna(False).astype(bool)
    df["is_robust"] = df["tier"] == "robust"
    df["is_swir_gap"] = df["category"].isin(SWIR_GAP_CATEGORIES)
    df["role"] = "other"
    df.loc[df["is_swir_gap"], "role"] = "swir_gap"
    df.loc[df["is_robust"], "role"] = "robust"
    return df


def _scatter_roles(ax, df: pd.DataFrame, size_scale: float = 1.0) -> None:
    """Plot sites: fill = coverage tier, shape = scientific role."""
    specs = [
        ("other", "o", 70 * size_scale, 0.6),
        ("swir_gap", "D", 90 * size_scale, 0.9),
        ("robust", "*", 430 * size_scale, 1.2),
    ]
    for role, marker, size, edge_w in specs:
        sub = df[df["role"] == role]
        if sub.empty:
            continue
        ax.scatter(
            sub["lon"], sub["lat"],
            s=size, marker=marker,
            c=sub["swir_tier"].map(TIER_COLORS).fillna("#BBBBBB"),
            edgecolors="#222222", linewidths=edge_w, zorder=4,
        )


def save_coverage_map(df: pd.DataFrame) -> Path:
    states = gpd.read_file(REPO_ROOT / "data" / "basemap" / "us_states.geojson")

    fig, (ax, ax_sw) = plt.subplots(
        1, 2, figsize=(15.5, 7.2),
        gridspec_kw={"width_ratios": [1.55, 1], "wspace": 0.06},
    )
    states.boundary.plot(ax=ax, color="#bcbcbc", linewidth=0.6, zorder=1)
    _scatter_roles(ax, df)

    # Label uncovered robust sites (the ones this data cannot speak to).
    uncovered_robust = df[df["is_robust"] & ~df["has_swir"]]
    pretty = {
        "holden": "Holden",
        "crooks_gap": "Crooks Gap",
        "karnes_county": "Karnes County",
    }
    offsets = {
        "holden": (6, 5),
        "crooks_gap": (6, 5),
        "karnes_county": (6, -10),
    }
    for _, r in uncovered_robust.iterrows():
        dx, dy = offsets.get(r["site_id"], (6, 4))
        ax.annotate(
            pretty.get(r["site_id"], r["site_id"]), (r["lon"], r["lat"]),
            textcoords="offset points", xytext=(dx, dy),
            fontsize=7.5, color="#222222", zorder=6,
        )

    ax.set_xlim(-125.5, -66)
    ax.set_ylim(24, 49.5)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect(1.30)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("CONUS", fontsize=10)
    ax.annotate(
        "N", xy=(-68.5, 47.6), xytext=(-68.5, 45.4),
        ha="center", va="bottom", fontsize=9, color="#333333",
        arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.1),
    )

    states.boundary.plot(ax=ax_sw, color="#bcbcbc", linewidth=0.6, zorder=1)
    _scatter_roles(ax_sw, df)
    az_pts = df[df["site_id"].isin(AZ_CLUSTER) & df["is_robust"]]
    if not az_pts.empty:
        ax_sw.annotate(
            "SE-AZ skarn belt (4)",
            (az_pts["lon"].mean(), az_pts["lat"].mean()),
            textcoords="offset points", xytext=(14, -28),
            fontsize=7, color="#222222",
            arrowprops=dict(arrowstyle="-", color="#888888", linewidth=0.8),
            zorder=6,
        )
    ax_sw.set_xlim(-125.0, -105.2)
    ax_sw.set_ylim(31.0, 43.5)
    ax_sw.set_xlabel("Longitude")
    ax_sw.set_ylabel("")
    ax_sw.set_aspect(1.30)
    ax_sw.spines["top"].set_visible(False)
    ax_sw.spines["right"].set_visible(False)
    ax_sw.set_title("Southwest (GEMx core)", fontsize=10)

    n = len(df)
    n_good = int((df["swir_tier"] == "good").sum())
    n_gap = int(df["is_swir_gap"].sum())
    n_gap_cov = int((df["is_swir_gap"] & df["has_swir"]).sum())
    n_rob = int(df["is_robust"].sum())
    n_rob_cov = int((df["is_robust"] & df["has_swir"]).sum())
    fig.suptitle(
        f"GEMx / AVIRIS SWIR coverage of {n} study sites\n"
        f"fill = bbox coverage  ·  {n_good} sites ≥75%  ·  "
        f"SWIR-gap cohort (epithermal + porphyry) {n_gap_cov}/{n_gap} covered  ·  "
        f"robust {n_rob_cov}/{n_rob} covered",
        fontsize=11, y=1.02,
    )

    tier_handles = [
        Patch(facecolor=TIER_COLORS["good"], edgecolor="#222222",
              label="good (≥75% of bbox)"),
        Patch(facecolor=TIER_COLORS["partial"], edgecolor="#222222",
              label="partial (25–75%)"),
        Patch(facecolor=TIER_COLORS["none"], edgecolor="#222222",
              label="none"),
    ]
    role_handles = [
        Line2D([0], [0], marker="*", linestyle="", markersize=15,
               markerfacecolor="#cccccc", markeredgecolor="#222222",
               label="fully robust (TIR result)"),
        Line2D([0], [0], marker="D", linestyle="", markersize=8,
               markerfacecolor="#cccccc", markeredgecolor="#222222",
               label="epithermal / porphyry (SWIR-gap test)"),
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               markerfacecolor="#cccccc", markeredgecolor="#222222",
               label="other deposit system"),
    ]
    fig.legend(
        handles=tier_handles + role_handles,
        loc="lower center", ncol=3, frameon=False, fontsize=8,
        bbox_to_anchor=(0.5, -0.08),
    )

    fig.subplots_adjust(top=0.86, bottom=0.14, wspace=0.12)
    out = REPO_ROOT / "figures" / "gemx_coverage_map.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")
    return out


def save_coverage_ranking(df: pd.DataFrame) -> Path:
    by_cat = (
        df.groupby("category", dropna=False)
        .agg(n=("site_id", "size"), n_covered=("has_swir", "sum"))
        .reset_index()
    )
    by_cat["pct"] = by_cat["n_covered"] / by_cat["n"] * 100
    by_cat = by_cat.sort_values("pct", ascending=True).reset_index(drop=True)
    by_cat["label"] = by_cat.apply(
        lambda r: f"{r['category']}  ({int(r['n_covered'])}/{int(r['n'])})",
        axis=1,
    )

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.barh(
        by_cat["label"], by_cat["pct"],
        color=by_cat["category"].map(CATEGORY_COLORS).fillna("#333333"),
        edgecolor="white", linewidth=0.5, height=0.7,
    )
    for i, r in by_cat.iterrows():
        ax.annotate(
            f"{r['pct']:.0f}%",
            (r["pct"], i),
            textcoords="offset points", xytext=(5, 0),
            va="center", fontsize=8, color="#333333",
        )
    ax.set_xlim(0, 115)
    ax.set_xlabel("% of sites with any GEMx-era SWIR coverage")
    n = len(df)
    n_cov = int(df["has_swir"].sum())
    ax.set_title(
        f"GEMx / AVIRIS coverage by deposit system, {n} sites "
        f"({n_cov} covered)\n"
        "A site is covered if any of its bbox is under AVIRIS-C or "
        "AVIRIS-5 flight lines from 2023 onward\n"
        "Coverage of the box, not a mineral identification",
        fontsize=11,
    )
    ax.tick_params(axis="y", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#eeeeee", linewidth=0.7, zorder=0)

    fig.tight_layout()
    out = REPO_ROOT / "figures" / "gemx_coverage_by_system.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")
    return out


def main() -> None:
    df = load_coverage()
    save_coverage_map(df)
    save_coverage_ranking(df)


if __name__ == "__main__":
    main()
