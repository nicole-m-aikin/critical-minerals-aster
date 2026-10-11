#!/usr/bin/env python
"""
Lithology-conditioned null: do anomaly zones and deposits overlap beyond what
shared host lithology alone predicts?

The primary site-specific null (p0 = zone area / valid TIR footprint area)
assumes a deposit is equally likely anywhere in the footprint. If strong
anomaly zones concentrate in, say, carbonate units and deposits also
concentrate in carbonate units, the overlap can be significant without ASTER
detecting any hydrothermal alteration (docs/publication_audit.md, reviewer
risk "lithology-conditioned randomization"). This script conditions the null
on USGS SGMC geologic-map lithology:

1. Lithology-stratified null. For each coarse lithology class L,
   p0_L = zone area in (L ∩ footprint) / area of (L ∩ footprint).
   Expected hits = Σ_L n_L · p0_L (n_L = critical-mineral deposits in L);
   exact one-sided Poisson-binomial p-value → p_lith; BH-FDR across all sites
   run → q_lith.
2. Host-unit-restricted test. Footprint, zones and deposits restricted to the
   classes hosting ≥ 1 deposit; binomial test on the restricted p0.
3. Descriptive zone-area vs deposit fractions by lithology class.
   Sensitivity: the same stratified test at SGMC map-unit level (SGMC_LABEL
   per state) — much finer strata, so conservative.

Deposit set, zones and footprint mirror scripts/regenerate_site_summaries.py
(validated MRDS → raster bbox → TIR footprint clip) restricted to critical
minerals as in scripts/site_specific_null_significance.py. Unconditioned
n / hits / p0 are cross-checked against results/phase4_fdr_corrected_significance.csv.

SGMC polygons are fetched from the ScienceBase ArcGIS FeatureServer layer 3
(SGMC_Geology), bbox-bounded and paged, and cached to
data/lithology/{site_id}_sgmc_units.geojson (gitignored; re-runs are offline).

Outputs
-------
    results/lithology_matched_null.csv            — one row per site
    results/lithology_zone_deposit_fractions.csv  — long format, site × class
    figures/verification_lithology_null.png       — 20 FDR sites, unconditioned vs conditioned

Usage:
    conda run -n aster-minerals python scripts/lithology_matched_null.py
    conda run -n aster-minerals python scripts/lithology_matched_null.py --sites bisbee holden
    conda run -n aster-minerals python scripts/lithology_matched_null.py --refresh   # re-fetch SGMC
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from shapely.ops import unary_union
from shapely.validation import make_valid
from statsmodels.stats.multitest import multipletests

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.lithology import (
    SGMC_SOURCE,
    UNMAPPED,
    add_lith_class,
    assign_strata,
    build_strata,
    host_restricted_test,
    load_sgmc_units,
    stratified_test,
    stratum_table,
)
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    is_critical_mineral,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.significance import run_binomial, significant_site_sets

REPO_ROOT = Path(__file__).parent.parent
ALPHA = 0.05

# FDR-significant on the primary binomial test and significant at every
# declustering radius, read from the current phase4 / phase5 results.
ROBUST_SITES = tuple(sorted(significant_site_sets(REPO_ROOT)[1]))
THEORY_VIOLATING = ("crooks_gap", "karnes_county", "holden")


def load_site_ids() -> list[str]:
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        return yaml.safe_load(f)["sites"]


def load_site_inputs(site_id: str):
    """Zones, footprint (zones CRS), raster bbox and critical deposits with hit flags."""
    site = load_site_by_id(site_id, REPO_ROOT / "sites")
    paths = site_paths_for(site, REPO_ROOT)
    prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
    if not paths.strong_zones_geojson.exists() or not prov_path.exists():
        raise FileNotFoundError("missing zones or provenance")
    zones = gpd.read_file(paths.strong_zones_geojson)
    if zones.empty:
        raise ValueError("empty zones")
    with open(prov_path) as f:
        prov = json.load(f)
    raster_bbox = tuple(prov["raster_bbox_wgs84"]) if prov.get("raster_bbox_wgs84") else site.bbox_wgs84
    if prov.get("tir_footprint_wgs84") is None:
        raise ValueError("no TIR footprint in provenance")
    footprint = gpd.GeoDataFrame.from_features(
        prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
    ).to_crs(zones.crs)

    local = filter_mrds_bbox(read_mrds_national(paths), raster_bbox)
    deposits = gpd.clip(mrds_to_points_gdf(local, zones.crs), footprint)
    deposits = reclassify_mrds_earth_mri(deposits)
    deposits = deposits[deposits["earth_mri_category"].apply(is_critical_mineral)].copy()
    _, hits, _ = spatial_join_deposits_zones(deposits, zones)
    deposits["hit"] = deposits.index.isin(hits.index.unique())
    return site, zones, footprint, raster_bbox, deposits


def analyse_site(site_id: str, refresh: bool) -> tuple[dict, list[dict]]:
    site, zones, footprint, raster_bbox, deposits = load_site_inputs(site_id)
    zone_union = make_valid(unary_union(zones.geometry.values))
    fp_geom = make_valid(unary_union(footprint.geometry.values))

    n = len(deposits)
    hits = int(deposits["hit"].sum())
    # Same denominator as significance.coverage_fraction_footprint.
    p0 = min(float(zone_union.area / footprint.iloc[0].geometry.area), 1.0)
    p_unc, exp_unc = run_binomial(hits, n, p0)

    units = load_sgmc_units(site_id, raster_bbox, REPO_ROOT, refresh=refresh)
    units = add_lith_class(units.to_crs(zones.crs)) if len(units) else units.assign(lith_class=[])
    strata = build_strata(units, fp_geom, "lith_class")
    labels = assign_strata(deposits, strata, "lith_class")
    tab = stratum_table(strata, zone_union, labels, deposits["hit"], "lith_class")
    strat = stratified_test(tab)
    host = host_restricted_test(tab, "lith_class")

    unit_res = {"p": float("nan"), "enrichment": float("nan"), "n_strata": 0}
    if len(units):
        units["unit"] = units["STATE"].astype(str) + ":" + units["SGMC_LABEL"].astype(str)
        u_strata = build_strata(units, fp_geom, "unit")
        u_tab = stratum_table(u_strata, zone_union, assign_strata(deposits, u_strata, "unit"),
                              deposits["hit"], "unit")
        u = stratified_test(u_tab)
        unit_res = {"p": u["p"], "enrichment": u["enrichment"], "n_strata": int((u_tab["area"] > 0).sum())}

    fp_area = float(tab["area"].sum())
    zone_area_fp = float(tab["zone_area"].sum())
    host_mask = tab["n"] > 0
    dep_tab = tab[host_mask].sort_values("n", ascending=False)
    dominant = dep_tab.iloc[0] if len(dep_tab) else None

    row = {
        "site_id": site_id,
        "site_name": site.name,
        "n": n,
        "hits": hits,
        "p0_unconditioned": round(p0, 4),
        "expected_unconditioned": round(exp_unc, 2),
        "enrichment_unconditioned": round((hits / n) / p0, 3) if n and p0 > 0 else np.nan,
        "p_unconditioned": p_unc,
        "expected_lith": round(strat["expected"], 2),
        "enrichment_lith": round(strat["enrichment"], 3) if n else np.nan,
        "p_lith": strat["p"],
        "p0_host": round(host["p0"], 4),
        "expected_host": round(host["expected"], 2),
        "enrichment_host": round(host["enrichment"], 3) if n else np.nan,
        "p_host": host["p"],
        "host_classes": ";".join(host["host_classes"]),
        "n_host_classes": len(host["host_classes"]),
        "host_area_frac": round(host["host_area_frac"], 4),
        "zone_share_in_host": round(float(tab.loc[host_mask, "zone_area"].sum()) / zone_area_fp, 4)
        if zone_area_fp > 0 else np.nan,
        "n_lith_classes": int(((tab["area"] > 0) & (tab["lith_class"] != UNMAPPED)).sum()),
        "dominant_host_lith": dominant["lith_class"] if dominant is not None else "",
        "dominant_host_share": round(dominant["n"] / n, 3) if dominant is not None and n else np.nan,
        "dominant_host_p0": round(float(dominant["p0"]), 4) if dominant is not None else np.nan,
        "unmapped_area_frac": round(float(tab.loc[tab["lith_class"] == UNMAPPED, "area"].sum()) / fp_area, 4)
        if fp_area > 0 else np.nan,
        "n_sgmc_units": len(units),
        "enrichment_unit": round(unit_res["enrichment"], 3) if n else np.nan,
        "p_unit": unit_res["p"],
        "n_unit_strata": unit_res["n_strata"],
        "data_source": SGMC_SOURCE,
    }

    frac_rows = []
    for _, r in tab.iterrows():
        frac_rows.append({
            "site_id": site_id,
            "lith_class": r["lith_class"],
            "area_km2": round(r["area"] / 1e6, 3),
            "zone_area_km2": round(r["zone_area"] / 1e6, 3),
            "p0_class": round(r["p0"], 4),
            "footprint_frac": round(r["area"] / fp_area, 4) if fp_area > 0 else np.nan,
            "zone_area_frac": round(r["zone_area"] / zone_area_fp, 4) if zone_area_fp > 0 else np.nan,
            "n_deposits": int(r["n"]),
            "deposit_frac": round(r["n"] / n, 4) if n else np.nan,
            "hits": int(r["hits"]),
            "hit_rate_class": round(r["hits"] / r["n"], 4) if r["n"] else np.nan,
        })
    return row, frac_rows


def add_fdr(df: pd.DataFrame, p_col: str, q_col: str, sig_col: str) -> None:
    valid = df[p_col].notna()
    df[q_col] = np.nan
    df[sig_col] = False
    if valid.any():
        reject, q, _, _ = multipletests(df.loc[valid, p_col].values, alpha=ALPHA, method="fdr_bh")
        df.loc[valid, q_col] = q
        df.loc[valid, sig_col] = reject


def make_figure(df: pd.DataFrame, out_path: Path) -> None:
    fdr = df[df["sig_fdr_ref"]].copy()
    fdr["robust"] = fdr["site_id"].isin(ROBUST_SITES)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7.2), gridspec_kw={"width_ratios": [1.15, 1]})

    lim = max(fdr["enrichment_unconditioned"].max(), fdr["enrichment_lith"].max()) * 1.1
    ax1.plot([0, lim], [0, lim], color="0.6", lw=1, ls="--", zorder=0, label="1:1 (no lithology effect)")
    ax1.axhline(1.0, color="0.8", lw=0.8, zorder=0)
    for robust, marker, colour, label in ((True, "*", "#c0392b", "Robust site (FDR + clustering)"),
                                          (False, "o", "#2c7fb8", "FDR-only site")):
        for sig in (True, False):
            sub = fdr[(fdr["robust"] == robust) & (fdr["sig_fdr_lith"] == sig)]
            if sub.empty:
                continue
            ax1.scatter(
                sub["enrichment_unconditioned"], sub["enrichment_lith"],
                marker=marker, s=190 if marker == "*" else 70,
                facecolors=colour if sig else "white", edgecolors=colour, linewidths=1.5,
                label=f"{label}, {'FDR-sig' if sig else 'not sig'} under lithology null",
                zorder=3,
            )
    # Points in the dense low-enrichment cluster get labels in a callout column (ordered by y so
    # leader lines don't cross); the rest are labelled in place.
    crowded = fdr[(fdr["enrichment_lith"] < 1.6) & (fdr["enrichment_unconditioned"] < 2.9)]
    crowded = crowded.sort_values("enrichment_lith", ascending=False)
    callout_x = max(3.15, float(crowded["enrichment_unconditioned"].max()) + 0.5) if len(crowded) else 0
    callout_y = dict(zip(crowded["site_id"], np.linspace(1.55, 0.25, len(crowded)))) if len(crowded) else {}
    for _, r in fdr.iterrows():
        special = r["site_id"] in THEORY_VIOLATING
        style = dict(fontsize=9.5 if special else 7.5, fontweight="bold" if special else "normal",
                     color="#7b1fa2" if special else "0.25")
        xy = (r["enrichment_unconditioned"], r["enrichment_lith"])
        if r["site_id"] in callout_y:
            ax1.annotate(r["site_id"], xy, xytext=(callout_x, callout_y[r["site_id"]]), ha="left", va="center",
                         arrowprops=dict(arrowstyle="-", color="0.6", lw=0.6, shrinkA=0, shrinkB=4), **style)
        else:
            ax1.annotate(r["site_id"], xy, xytext=(5, 4), textcoords="offset points", ha="left", **style)
    ax1.set_xlim(0.8, lim)
    ax1.set_ylim(0.0, lim)
    ax1.set_xlabel("Enrichment vs site-specific geometric null (observed / expected)")
    ax1.set_ylabel("Enrichment vs lithology-stratified null (SGMC classes)")
    ax1.set_title("Enrichment before vs after conditioning on host lithology")
    ax1.legend(loc="upper left", fontsize=8, frameon=True)
    ax1.grid(alpha=0.25)

    order = fdr.sort_values("q_lith", ascending=False)
    y = np.arange(len(order))
    floor = 1e-20
    q_unc = -np.log10(order["q_unconditioned"].clip(lower=floor))
    q_lith = -np.log10(order["q_lith"].clip(lower=floor))
    colours = ["#c0392b" if s in ROBUST_SITES else "#2c7fb8" for s in order["site_id"]]
    ax2.barh(y + 0.2, q_unc, height=0.4, color="0.78")
    ax2.barh(y - 0.2, q_lith, height=0.4, color=colours)
    ax2.axvline(-np.log10(ALPHA), color="k", lw=1, ls="--")
    ax2.set_yticks(y)
    ax2.set_yticklabels(order["site_id"], fontsize=8.5)
    ax2.set_ylim(-0.6, len(order) - 0.4)
    for tick, sid in zip(ax2.get_yticklabels(), order["site_id"]):
        if sid in THEORY_VIOLATING:
            tick.set_color("#7b1fa2")
            tick.set_fontweight("bold")
    ax2.set_xlabel(f"−log10 BH q-value (across {len(df)} sites; capped at q = {floor:g})")
    ax2.set_title("FDR significance: geometric vs lithology-conditioned")
    handles = [
        plt.Rectangle((0, 0), 1, 1, color="0.78"),
        plt.Rectangle((0, 0), 1, 1, color="#c0392b"),
        plt.Rectangle((0, 0), 1, 1, color="#2c7fb8"),
        plt.Line2D([0], [0], color="k", lw=1, ls="--"),
    ]
    labels = [
        "Geometric null (this run)",
        "Lithology-stratified null — robust site",
        "Lithology-stratified null — FDR-only site",
        "q = 0.05",
    ]
    ax2.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, fontsize=8,
               frameon=False)
    ax2.grid(axis="x", alpha=0.25)

    fig.suptitle(
        f"Lithology-conditioned null (USGS SGMC, ~1:500k) — {int(df['sig_fdr_ref'].sum())} primary-FDR sites",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--sites", nargs="*", help="subset of site ids (default: all in sites/index.yaml)")
    ap.add_argument("--refresh", action="store_true", help="re-fetch SGMC polygons, ignoring cache")
    ap.add_argument("--figure-only", action="store_true",
                    help="redraw the figure from results/lithology_matched_null.csv")
    args = ap.parse_args()

    fig_path = REPO_ROOT / "figures" / "verification_lithology_null.png"
    if args.figure_only:
        make_figure(pd.read_csv(REPO_ROOT / "results" / "lithology_matched_null.csv"), fig_path)
        print(f"Saved {fig_path}")
        return

    site_ids = args.sites or load_site_ids()
    ref = pd.read_csv(REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv").set_index("site_id")

    rows, frac_rows, failed = [], [], []
    for site_id in site_ids:
        t0 = time.time()
        print(f"{site_id:24s}", end=" ", flush=True)
        try:
            row, fr = analyse_site(site_id, args.refresh)
        except Exception as exc:
            print(f"SKIP — {exc!r}")
            failed.append((site_id, repr(exc)))
            continue
        rows.append(row)
        frac_rows.extend(fr)
        print(
            f"n={row['n']:4d} hits={row['hits']:3d} p0={row['p0_unconditioned']:.4f} "
            f"E_unc={row['enrichment_unconditioned']:.2f} E_lith={row['enrichment_lith']:.2f} "
            f"p_lith={row['p_lith']:.3g} p_host={row['p_host']:.3g} "
            f"host={row['dominant_host_lith']} ({time.time() - t0:.1f}s)"
        )

    df = pd.DataFrame(rows)
    df["n_ref"] = df["site_id"].map(ref["n_crit"])
    df["hits_ref"] = df["site_id"].map(ref["hits_crit"])
    df["p0_ref"] = df["site_id"].map(ref["p0_geometric"] if "p0_geometric" in ref.columns else ref["p0_null"])
    df["sig_fdr_ref"] = df["site_id"].map(ref["sig_fdr_binomial"]).fillna(False).astype(bool)
    df["robust"] = df["site_id"].isin(ROBUST_SITES)
    add_fdr(df, "p_unconditioned", "q_unconditioned", "sig_fdr_unconditioned")
    add_fdr(df, "p_lith", "q_lith", "sig_fdr_lith")
    add_fdr(df, "p_host", "q_host", "sig_fdr_host")
    add_fdr(df, "p_unit", "q_unit", "sig_fdr_unit")

    lead = [
        "site_id", "site_name", "robust", "sig_fdr_ref",
        "n", "hits", "n_ref", "hits_ref", "p0_unconditioned", "p0_ref",
        "expected_unconditioned", "enrichment_unconditioned", "p_unconditioned",
        "q_unconditioned", "sig_fdr_unconditioned",
        "expected_lith", "enrichment_lith", "p_lith", "q_lith", "sig_fdr_lith",
        "p0_host", "expected_host", "enrichment_host", "p_host", "q_host", "sig_fdr_host",
    ]
    df = df[lead + [c for c in df.columns if c not in lead]].sort_values("p_lith")

    out_csv = REPO_ROOT / "results" / "lithology_matched_null.csv"
    df.to_csv(out_csv, index=False)
    fr_csv = REPO_ROOT / "results" / "lithology_zone_deposit_fractions.csv"
    pd.DataFrame(frac_rows).to_csv(fr_csv, index=False)
    print(f"\nSaved {out_csv}\nSaved {fr_csv}")

    mism = df[(df["n"] != df["n_ref"]) | (df["hits"] != df["hits_ref"]) | ((df["p0_unconditioned"] - df["p0_ref"]).abs() > 5e-4)]
    print(f"\nUnconditioned cross-check vs phase4 file: {len(df) - len(mism)}/{len(df)} exact on n/hits/p0")
    if not mism.empty:
        print(mism[["site_id", "n", "n_ref", "hits", "hits_ref", "p0_unconditioned", "p0_ref"]].to_string(index=False))

    fdr = df[df["sig_fdr_ref"]]
    cols = ["site_id", "robust", "n", "hits", "enrichment_unconditioned", "p_unconditioned",
            "enrichment_lith", "p_lith", "q_lith", "sig_fdr_lith", "enrichment_host", "p_host",
            "dominant_host_lith", "enrichment_unit", "p_unit", "q_unit"]
    print(f"\n{len(fdr)} primary-FDR sites:")
    print(fdr[cols].to_string(index=False))
    print(f"\nFDR-significant under lithology null: {int(df['sig_fdr_lith'].sum())}/{len(df)} "
          f"(geometric null this run: {int(df['sig_fdr_unconditioned'].sum())})")
    print(f"Robust sites surviving lithology null: "
          f"{sorted(df.loc[df['robust'] & df['sig_fdr_lith'], 'site_id'])}")
    if failed:
        print(f"\nFailed/skipped ({len(failed)}): {failed}")

    if not fdr.empty:
        make_figure(df, fig_path)
        print(f"Saved {fig_path}")


if __name__ == "__main__":
    main()
