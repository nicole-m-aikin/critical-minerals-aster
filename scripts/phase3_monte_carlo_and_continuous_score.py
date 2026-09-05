#!/usr/bin/env python
"""
Phase 3 — Monte Carlo spatial null + threshold-free continuous-score check.

Two independent robustness checks on the Phase 2 primary result
(results/site_specific_null_significance.csv), both against the same
site-specific null geometry:

1. MONTE CARLO SPATIAL PERMUTATION (run_permutation, now footprint-aware).
   Preserves: the valid TIR footprint polygon shape, n_crit (deposit count).
   Randomizes: the locations of n_crit synthetic points, drawn uniformly
   from footprint grid cells only (not the rectangular bbox — see the
   docstring added to run_permutation() for why this now matches the
   binomial test's denominator). 10,000 iterations, seed=42.
   This is a distribution-free cross-check of the analytical binomial
   p-value — it does not assume a fixed p and independent Bernoulli trials
   the way binomtest's *formula* does; it directly simulates the null
   point process instead.

2. CONTINUOUS-SCORE MANN-WHITNEY U TEST (new).
   Motivation (docs/results.md — "audit whether the method needs to
   change"): the binary hit-rate test depends on one specific choice of
   classification threshold (70th/90th percentile, combined score >= 3).
   Phase 2 showed the exact denominator this threshold produces varies
   9.0-14.7% across sites for reasons unrelated to deposit clustering
   (spatial correlation among the three band ratios) -- proof the binary
   zone is somewhat threshold-sensitive. This test asks the same question
   without needing any threshold at all: do deposits sit on pixels with a
   HIGHER continuous combined_score (0-6, silica+carbonate+mafic class sum)
   than the rest of the surveyed footprint, using every pixel's exact score
   rather than collapsing it to "in a top-10%-ish zone or not"?
   Uses scipy.stats.mannwhitneyu (Wilcoxon rank-sum), one-sided
   ("greater"), on deposit-pixel scores vs. a random background sample of
   footprint pixels. Effect size reported as the common-language / AUC
   statistic: P(a random deposit pixel scores higher than a random
   background pixel), where 0.5 = no discrimination, 1.0 = perfect
   separation. This requires no percentile threshold at all, so a
   significant result here that agrees with the binary test is strong
   evidence the finding is not a threshold artifact.

Reads combined-score rasters from the cached mosaic classification files
(data/sites/{id}/aster/{id}_mosaic_cls_{silica,carbonate,mafic,valid}.tif)
that the pipeline already wrote when zones were generated — no raster
reprocessing, no network calls.

Outputs
-------
    results/phase3_monte_carlo_and_continuous_score.csv

Usage:
    conda run -n aster-minerals python scripts/phase3_monte_carlo_and_continuous_score.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from scipy.stats import mannwhitneyu

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.classification import classify_percentiles, combined_score
from critical_minerals_aster.spectral import load_ratio_mosaic, clip_bands_to_bbox
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
)
from critical_minerals_aster.significance import run_permutation

REPO_ROOT = Path(__file__).parent.parent
N_ITER = 10_000
SEED = 42
BACKGROUND_SAMPLE = 50_000


def load_combined_score_raster(aster_dir: Path, site_id: str, site, raster_bbox):
    """Classify fresh from the cached ratio mosaic (mosaic_ratio_{silica,
    carbonate,mafic}.tif — present for all 45 sites) using the site's own
    classification thresholds, rather than trusting cached _cls_*.tif files.

    Several sites (jerome, morenci, randsburg, eureka, lordsburg, sierrita,
    ray_mine, bodie, patagonia) had their old per-granule _cls_*.tif files
    deleted when the mosaic-level classification bug was fixed (see
    docs/results.md "Methodology note: classification approach and bug
    fix") and were never regenerated on disk, even though the ratio mosaic
    and the current zones geojson both reflect the corrected method.
    Reclassifying directly from the ratio mosaic guarantees every site uses
    the identical, current, mosaic-level classification consistently,
    rather than depending on which intermediate files happen to still
    exist. Returns None only if the ratio mosaic itself is missing.

    CRITICAL: the ratio mosaic .tif on disk covers a much larger raw extent
    than the site's configured bbox (multi-granule mosaics are built by
    stitching full granule footprints, not just the target bbox — measured
    at 1.6x-19x larger across all 45 sites). pipeline.py always calls
    clip_bands_to_bbox(..., site.bbox_wgs84) before classification
    (pipeline.py's build-zones path) -- an earlier version of this script
    skipped that step, which both inflated the apparent "valid footprint"
    area and computed percentile thresholds over a much larger, wrong
    population than the real pipeline uses, producing classification
    results that disagreed with the cached zones by up to several-fold at
    some sites. See docs/results.md for the full
    story of how this was caught (it looked at first like the CACHED zones
    were stale; reverting a speculative zone-regeneration and checking
    scope revealed the bug was in this script, not the cached data).
    """
    prefix = f"{site_id}_mosaic"
    ratio_paths = {name: aster_dir / f"{prefix}_ratio_{name}.tif" for name in
                    ("silica", "carbonate", "mafic")}
    if not all(p.is_file() for p in ratio_paths.values()):
        return None

    silica, carbonate, mafic, b10, transform, crs = load_ratio_mosaic(aster_dir, prefix)
    (silica, carbonate, mafic), transform = clip_bands_to_bbox(
        [silica, carbonate, mafic], transform, crs, raster_bbox
    )

    cp = site.classification
    s_cls, _, _ = classify_percentiles(silica, cp.low_pct, cp.high_pct)
    c_cls, _, _ = classify_percentiles(carbonate, cp.low_pct, cp.high_pct)
    m_cls, _, _ = classify_percentiles(mafic, cp.low_pct, cp.high_pct)

    score = combined_score(s_cls, c_cls, m_cls)
    valid = np.isfinite(silica) & np.isfinite(carbonate) & np.isfinite(mafic)
    return score, valid, transform, crs


def main() -> None:
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    rows = []

    for _, ss_row in ss.iterrows():
        site_id = ss_row["site_id"]
        print(f"{site_id:22s}", end=" ", flush=True)

        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        paths = site_paths_for(site, REPO_ROOT)

        zones_path = paths.strong_zones_geojson
        prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
        if not zones_path.exists() or not prov_path.exists():
            print("SKIP — no zones/provenance")
            continue
        zones = gpd.read_file(zones_path)
        with open(prov_path) as f:
            prov = json.load(f)
        raster_bbox = tuple(prov["raster_bbox_wgs84"])
        footprint = None
        if prov.get("tir_footprint_wgs84") is not None:
            try:
                footprint = gpd.GeoDataFrame.from_features(
                    prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
                ).to_crs(zones.crs)
            except Exception:
                footprint = None

        n_crit = int(ss_row["n_crit"])
        hits_crit = int(ss_row["hits_crit"])
        p_binomial = float(ss_row["p_binomial"])

        # --- 1. Footprint-aware Monte Carlo spatial permutation ---
        p_perm = run_permutation(
            zones, raster_bbox, n_crit, hits_crit,
            n_iter=N_ITER, seed=SEED, footprint=footprint,
        )
        p_perm_label = f"<{1/N_ITER:.0e}" if p_perm == 0.0 else f"{p_perm:.4f}"

        # --- 2. Continuous-score Mann-Whitney U test ---
        cs = load_combined_score_raster(paths.aster_dir, site_id, site, raster_bbox)
        mwu_p, mwu_auc, n_deposit_px = (float("nan"), float("nan"), 0)
        if cs is not None and n_crit > 0:
            score, valid, transform, crs = cs

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

            if len(deposits):
                rows_px, cols_px = rasterio.transform.rowcol(
                    transform, deposits.geometry.x.values, deposits.geometry.y.values
                )
                rows_px = np.asarray(rows_px)
                cols_px = np.asarray(cols_px)
                in_bounds = (
                    (rows_px >= 0) & (rows_px < score.shape[0])
                    & (cols_px >= 0) & (cols_px < score.shape[1])
                )
                rows_px, cols_px = rows_px[in_bounds], cols_px[in_bounds]
                px_valid = valid[rows_px, cols_px]
                deposit_scores = score[rows_px, cols_px][px_valid].astype(float)
                n_deposit_px = len(deposit_scores)

                background_all = score[valid].astype(float)
                rng = np.random.default_rng(SEED)
                if len(background_all) > BACKGROUND_SAMPLE:
                    background = rng.choice(background_all, size=BACKGROUND_SAMPLE, replace=False)
                else:
                    background = background_all

                if n_deposit_px >= 3 and len(background) >= 3:
                    u_stat, mwu_p = mannwhitneyu(
                        deposit_scores, background, alternative="greater"
                    )
                    mwu_auc = float(u_stat) / (n_deposit_px * len(background))

        rows.append({
            "site_id": site_id,
            "n_crit": n_crit,
            "hits_crit": hits_crit,
            "p_binomial": round(p_binomial, 6),
            "p_monte_carlo": round(p_perm, 4),
            "p_monte_carlo_label": p_perm_label,
            "n_deposit_pixels_scored": n_deposit_px,
            "p_mannwhitney": round(mwu_p, 6) if not np.isnan(mwu_p) else None,
            "auc_continuous_score": round(mwu_auc, 4) if not np.isnan(mwu_auc) else None,
        })
        print(
            f"p_binom={p_binomial:.4g}  p_MC={p_perm_label}  "
            f"p_MWU={'n/a' if np.isnan(mwu_p) else f'{mwu_p:.4g}'}  "
            f"AUC={'n/a' if np.isnan(mwu_auc) else f'{mwu_auc:.3f}'}"
        )

    out = pd.DataFrame(rows).sort_values("p_binomial")
    out_path = REPO_ROOT / "results" / "phase3_monte_carlo_and_continuous_score.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
