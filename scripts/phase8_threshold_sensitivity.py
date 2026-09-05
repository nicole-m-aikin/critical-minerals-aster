#!/usr/bin/env python
"""
Phase 8a -- threshold-perturbation and anomaly-score-cutoff sensitivity.

Two parameters in the classification pipeline are empirical, not physical
constants (see src/critical_minerals_aster/config.py's ClassificationParams
and docs/results.md):

  1. The percentile thresholds (low_pct=70, high_pct=90) that turn each
     continuous band ratio into a 3-class map (background/moderate/strong).
  2. The combined-score cutoff (strong_score_min=3 of a possible 0-6) that
     turns the 3 per-ratio class maps into a single "strong anomaly" mask.

Neither is literature-derived or theoretically required -- they were
initial, reasonable-looking choices (a top-30%/top-10% split; "majority of
the three ratios must agree") that this script tests for sensitivity.
Unlike Phase 3's continuous-score test (which removes thresholds entirely),
this asks a narrower, more direct question: if the SAME kind of binary
decision rule is kept but its two knobs are turned, does the significance
conclusion change?

Method: for each scenario, reclassify fresh from the cached ratio mosaic
(same infrastructure as scripts/phase3_monte_carlo_and_continuous_score.py),
recompute the pixel-level coverage fraction p0 = strong_pixels / valid_pixels
for THAT scenario specifically (the null must be recomputed per scenario --
a looser threshold produces a larger zone and therefore a larger p0), sample
each critical-mineral deposit's pixel against the new mask, and rerun the
one-sided binomial test.

Sweep A (percentile thresholds, score_min fixed at nominal 3):
    (65, 85) loose -- (70, 90) nominal -- (75, 95) tight
Sweep B (score cutoff, percentiles fixed at nominal 70/90):
    score_min in {2, 3, 4} of 6

Outputs
-------
    results/phase8_threshold_sensitivity.csv
    results/phase8_score_cutoff_sensitivity.csv
    figures/phase8_threshold_sensitivity.png

Usage:
    conda run -n aster-minerals python scripts/phase8_threshold_sensitivity.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.classification import classify_percentiles, combined_score
from critical_minerals_aster.spectral import load_ratio_mosaic, clip_bands_to_bbox
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import filter_mrds_bbox, mrds_to_points_gdf, reclassify_mrds_earth_mri
from critical_minerals_aster.significance import run_binomial

REPO_ROOT = Path(__file__).parent.parent
PCT_SCENARIOS = [(65.0, 85.0, "loose"), (70.0, 90.0, "nominal"), (75.0, 95.0, "tight")]
SCORE_SCENARIOS = [2, 3, 4]

# Post-109-site-expansion: a 9-site readable subset (matplotlib's 10-colour
# cycle wraps beyond this) -- original robust skarn + 3 new declustering-robust
# skarn replication sites + 2 robust "surprises" + eureka as a now-fragile
# contrast. The full per-site scenario grid is in results/phase8_threshold_sensitivity.csv.
KEY_SITES = ["bisbee", "tombstone", "courtland_gleeson",
             "magdalena_kelly", "organ_mountains", "eagle_mountain",
             "crooks_gap", "holden", "eureka"]


def load_deposit_pixels(paths, site, raster_bbox, footprint, transform, crs, shape):
    mrds = read_mrds_national(paths)
    local = filter_mrds_bbox(mrds, raster_bbox)
    deposits = mrds_to_points_gdf(local, crs)
    if footprint is not None and len(footprint):
        try:
            deposits = gpd.clip(deposits, footprint)
        except Exception:
            pass
    deposits = reclassify_mrds_earth_mri(deposits)
    deposits = deposits[deposits["earth_mri_category"] != "Non-Critical"]
    if deposits.empty:
        return None, None
    rows_px, cols_px = rasterio.transform.rowcol(
        transform, deposits.geometry.x.values, deposits.geometry.y.values
    )
    rows_px, cols_px = np.asarray(rows_px), np.asarray(cols_px)
    in_bounds = (rows_px >= 0) & (rows_px < shape[0]) & (cols_px >= 0) & (cols_px < shape[1])
    return rows_px[in_bounds], cols_px[in_bounds]


def main() -> None:
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    pct_rows, score_rows = [], []

    for _, ss_row in ss.iterrows():
        site_id = ss_row["site_id"]
        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        paths = site_paths_for(site, REPO_ROOT)
        prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
        if not prov_path.exists():
            continue
        with open(prov_path) as f:
            prov = json.load(f)
        raster_bbox = tuple(prov["raster_bbox_wgs84"])

        prefix = f"{site_id}_mosaic"
        ratio_paths = {n: paths.aster_dir / f"{prefix}_ratio_{n}.tif" for n in
                        ("silica", "carbonate", "mafic")}
        if not all(p.is_file() for p in ratio_paths.values()):
            print(f"{site_id:22s} SKIP — no ratio mosaic")
            continue
        silica, carbonate, mafic, _, transform, crs = load_ratio_mosaic(paths.aster_dir, prefix)
        # Clip to the site bbox BEFORE classification -- matches pipeline.py exactly
        # and avoids computing percentiles over the much larger raw mosaic extent
        # (measured at 1.6x-19x the true footprint area). An earlier version of this
        # script and phase3's script both skipped this clip, producing classification
        # results that disagreed with the cached zones by up to several-fold.
        (silica, carbonate, mafic), transform = clip_bands_to_bbox(
            [silica, carbonate, mafic], transform, crs, raster_bbox
        )
        valid = np.isfinite(silica) & np.isfinite(carbonate) & np.isfinite(mafic)

        footprint = None
        if prov.get("tir_footprint_wgs84") is not None:
            try:
                footprint = gpd.GeoDataFrame.from_features(
                    prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
                ).to_crs(crs)
            except Exception:
                footprint = None

        rows_px, cols_px = load_deposit_pixels(paths, site, raster_bbox, footprint,
                                                 transform, crs, silica.shape)
        if rows_px is None or len(rows_px) == 0:
            print(f"{site_id:22s} SKIP — no critical deposits")
            continue

        print(f"{site_id:22s}", end=" ", flush=True)

        # --- Sweep A: percentile thresholds, score_min fixed at 3 ---
        for low_pct, high_pct, label in PCT_SCENARIOS:
            s_cls, _, _ = classify_percentiles(silica, low_pct, high_pct)
            c_cls, _, _ = classify_percentiles(carbonate, low_pct, high_pct)
            m_cls, _, _ = classify_percentiles(mafic, low_pct, high_pct)
            score = combined_score(s_cls, c_cls, m_cls)
            mask = (score >= 3) & valid

            p0 = float(mask.sum()) / float(valid.sum()) if valid.sum() else 0.0
            hit_vals = mask[rows_px, cols_px]
            n, hits = len(hit_vals), int(hit_vals.sum())
            p_val, _ = run_binomial(hits, n, p0)
            enrichment = (hits / n) / p0 if n and p0 else float("nan")

            pct_rows.append({
                "site_id": site_id, "scenario": label, "low_pct": low_pct, "high_pct": high_pct,
                "p0_null": round(p0, 4), "n_crit": n, "hits_crit": hits,
                "hit_rate_pct": round(hits / n * 100, 2) if n else 0.0,
                "enrichment": round(enrichment, 2), "p_binomial": round(p_val, 6),
                "sig_05": p_val < 0.05,
            })

        # --- Sweep B: score cutoff, percentiles fixed at nominal 70/90 ---
        s_cls, _, _ = classify_percentiles(silica, 70.0, 90.0)
        c_cls, _, _ = classify_percentiles(carbonate, 70.0, 90.0)
        m_cls, _, _ = classify_percentiles(mafic, 70.0, 90.0)
        score = combined_score(s_cls, c_cls, m_cls)
        for score_min in SCORE_SCENARIOS:
            mask = (score >= score_min) & valid
            p0 = float(mask.sum()) / float(valid.sum()) if valid.sum() else 0.0
            hit_vals = mask[rows_px, cols_px]
            n, hits = len(hit_vals), int(hit_vals.sum())
            p_val, _ = run_binomial(hits, n, p0)
            enrichment = (hits / n) / p0 if n and p0 else float("nan")

            score_rows.append({
                "site_id": site_id, "score_min": score_min,
                "p0_null": round(p0, 4), "n_crit": n, "hits_crit": hits,
                "hit_rate_pct": round(hits / n * 100, 2) if n else 0.0,
                "enrichment": round(enrichment, 2), "p_binomial": round(p_val, 6),
                "sig_05": p_val < 0.05,
            })

        nominal_p = [r for r in pct_rows if r["site_id"] == site_id and r["scenario"] == "nominal"][0]["p_binomial"]
        print(f"p(nominal)={nominal_p:.4g}")

    pct_df = pd.DataFrame(pct_rows)
    score_df = pd.DataFrame(score_rows)
    pct_df.to_csv(REPO_ROOT / "results" / "phase8_threshold_sensitivity.csv", index=False)
    score_df.to_csv(REPO_ROOT / "results" / "phase8_score_cutoff_sensitivity.csv", index=False)
    print(f"\nSaved results/phase8_threshold_sensitivity.csv and phase8_score_cutoff_sensitivity.csv")

    # --- Stability summary for key sites ---
    print("\n=== Percentile-threshold stability (key sites) ===")
    piv = pct_df[pct_df.site_id.isin(KEY_SITES)].pivot(index="site_id", columns="scenario", values="p_binomial")
    piv = piv[["loose", "nominal", "tight"]]
    print(piv.to_string())

    print("\n=== Score-cutoff stability (key sites) ===")
    piv2 = score_df[score_df.site_id.isin(KEY_SITES)].pivot(index="site_id", columns="score_min", values="p_binomial")
    print(piv2.to_string())

    # --- Figure ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))
    colors = plt.cm.tab10(np.linspace(0, 1, len(KEY_SITES)))
    color_map = dict(zip(KEY_SITES, colors))

    ax = axes[0]
    for site_id in KEY_SITES:
        sub = pct_df[pct_df.site_id == site_id]
        order = {"loose": 0, "nominal": 1, "tight": 2}
        sub = sub.assign(_x=sub["scenario"].map(order)).sort_values("_x")
        y = sub["p_binomial"].clip(lower=1e-6)
        ax.plot(sub["_x"], y, marker="o", linewidth=2, color=color_map[site_id], label=site_id)
    ax.axhline(0.05, color="#555555", linewidth=1, linestyle="--")
    ax.set_yscale("log")
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["loose\n(65/85)", "nominal\n(70/90)", "tight\n(75/95)"])
    ax.set_ylabel("Binomial p-value (log scale)")
    ax.set_title("Percentile-threshold sensitivity")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    for site_id in KEY_SITES:
        sub = score_df[score_df.site_id == site_id].sort_values("score_min")
        y = sub["p_binomial"].clip(lower=1e-6)
        ax.plot(sub["score_min"], y, marker="o", linewidth=2, color=color_map[site_id], label=site_id)
    ax.axhline(0.05, color="#555555", linewidth=1, linestyle="--")
    ax.set_yscale("log")
    ax.set_xticks([2, 3, 4])
    ax.set_xlabel("Combined-score cutoff (of 6)")
    ax.set_title("Anomaly-score-cutoff sensitivity")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False, fontsize=8)

    fig.suptitle("Phase 8 — classification-parameter sensitivity, key sites", fontsize=12)
    fig.tight_layout()
    out_path = REPO_ROOT / "figures" / "phase8_threshold_sensitivity.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
