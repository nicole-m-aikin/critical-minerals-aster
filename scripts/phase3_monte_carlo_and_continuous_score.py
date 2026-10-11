#!/usr/bin/env python
"""
Phase 3 — Monte Carlo spatial null + threshold-free continuous-score check,
both STRATIFIED by bedrock / cover / unmapped ground (the Phase 2 primary null).

Two independent robustness checks on the Phase 2 primary result
(results/site_specific_null_significance.csv). Strata come from
critical_minerals_aster.lithology.cover_strata (USGS SGMC: cover =
"Unconsolidated…" units, unmapped = no SGMC polygon, bedrock = the rest),
clipped to the valid TIR footprint.

1. MONTE CARLO SPATIAL PERMUTATION (run_permutation_stratified).
   Preserves: the footprint and the number of critical deposits in EACH
   stratum (n_bedrock, n_cover, n_unmapped — from deposit locations).
   Randomizes: each stratum's n_s points, drawn uniformly from grid cells of
   that stratum only. 10,000 iterations, seed=42. This is the simulation
   counterpart of Phase 2's Poisson-binomial p-value: a deposit on alluvium
   is only ever compared with random alluvium locations, so a zone set that
   simply favours outcrop cannot look enriched.
   The previous footprint-wide permutation is kept as p_monte_carlo_geometric.

2. CONTINUOUS-SCORE STRATIFIED MANN-WHITNEY U (stratified_mannwhitney).
   Asks the threshold-free version of the question: do deposits sit on
   pixels with a higher continuous combined_score (0-6) than background?
   Design choice — WITHIN-STRATUM RANKS, not a stratum-weighted pooled
   background: within each stratum s, deposit-pixel scores are ranked only
   against background pixels of the same stratum (U_s); strata combine
   with deposit-count weights w_s = n_s / n into a stratified AUC
   (P(deposit pixel outscores a random background pixel from its own
   stratum)), with a tie-corrected normal approximation (combined_score
   is an integer 0-6, so ties are pervasive). A pooled MWU against a
   background resampled to the deposits' stratum mix has the right null
   mean but still scores bedrock deposits against cover pixels (and vice
   versa); the within-stratum form removes every cross-stratum comparison,
   so a detector that just scores outcrop above valley fill contributes
   nothing. Background per stratum: all valid pixels of that stratum,
   randomly subsampled to 50,000 (seed=42) when larger. Pixels and deposits
   take their stratum from the same rasterized strata grid. The previous
   unstratified test is kept as p_mannwhitney_geometric /
   auc_continuous_score_geometric.

Reads combined-score rasters freshly classified from the cached ratio
mosaic (data/sites/{id}/aster/{id}_mosaic_ratio_{silica,carbonate,mafic}.tif).

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
from shapely.geometry import box
from shapely.ops import unary_union
from shapely.validation import make_valid

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.classification import classify_percentiles, combined_score
from critical_minerals_aster.lithology import (
    COVER_STRATA,
    assign_cover_stratum,
    cover_strata,
    rasterize_cover_strata,
)
from critical_minerals_aster.spectral import load_ratio_mosaic, clip_bands_to_bbox
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
)
from critical_minerals_aster.significance import (
    run_permutation,
    run_permutation_stratified,
    stratified_mannwhitney,
)

REPO_ROOT = Path(__file__).parent.parent
N_ITER = 10_000
SEED = 42
BACKGROUND_SAMPLE = 50_000


def load_combined_score_raster(aster_dir: Path, site_id: str, site, raster_bbox):
    """Classify fresh from the cached ratio mosaic (mosaic_ratio_{silica,
    carbonate,mafic}.tif) using the site's own classification thresholds,
    rather than trusting cached _cls_*.tif files.

    Several sites had their old per-granule _cls_*.tif files deleted when the
    mosaic-level classification bug was fixed (docs/results.md "Methodology
    note: classification approach and bug fix"); reclassifying from the ratio
    mosaic guarantees every site uses the identical, current classification.
    Returns None only if the ratio mosaic itself is missing.

    CRITICAL: the ratio mosaic covers a much larger raw extent than the
    site's bbox (1.6x-19x). pipeline.py clips to the bbox before
    classification, so this does too — skipping the clip both inflates the
    valid footprint and computes percentile thresholds over the wrong
    population (docs/results.md).
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


def load_critical_deposits(paths, raster_bbox, footprint, crs) -> gpd.GeoDataFrame:
    deposits = mrds_to_points_gdf(filter_mrds_bbox(read_mrds_national(paths), raster_bbox), crs)
    if footprint is not None and len(footprint):
        try:
            deposits = gpd.clip(deposits, footprint.to_crs(crs))
        except Exception:
            pass
    deposits = reclassify_mrds_earth_mri(deposits)
    return deposits[deposits["earth_mri_category"] != "Non-Critical"].copy()


def continuous_score_tests(score, valid, transform, deposits, strata, rng):
    """(n_deposit_px, stratified MWU dict, geometric p, geometric AUC)."""
    rows_px, cols_px = rasterio.transform.rowcol(
        transform, deposits.geometry.x.values, deposits.geometry.y.values
    )
    rows_px, cols_px = np.asarray(rows_px), np.asarray(cols_px)
    in_bounds = (rows_px >= 0) & (rows_px < score.shape[0]) & (cols_px >= 0) & (cols_px < score.shape[1])
    rows_px, cols_px = rows_px[in_bounds], cols_px[in_bounds]
    px_valid = valid[rows_px, cols_px]
    rows_px, cols_px = rows_px[px_valid], cols_px[px_valid]
    deposit_scores = score[rows_px, cols_px].astype(float)
    n_deposit_px = len(deposit_scores)

    # --- Geometric (unstratified) MWU: previous primary continuous test ---
    p_geo, auc_geo = float("nan"), float("nan")
    background_all = score[valid].astype(float)
    background = (rng.choice(background_all, size=BACKGROUND_SAMPLE, replace=False)
                  if len(background_all) > BACKGROUND_SAMPLE else background_all)
    if n_deposit_px >= 3 and len(background) >= 3:
        u_stat, p_geo = mannwhitneyu(deposit_scores, background, alternative="greater")
        auc_geo = float(u_stat) / (n_deposit_px * len(background))

    # --- Stratified MWU: deposits ranked only against their own stratum ---
    codes = rasterize_cover_strata(strata, transform, score.shape)
    dep_codes = codes[rows_px, cols_px]
    bg_scores, bg_codes = [], []
    for code in np.unique(dep_codes):
        pool = score[valid & (codes == code)].astype(float)
        if len(pool) > BACKGROUND_SAMPLE:
            pool = rng.choice(pool, size=BACKGROUND_SAMPLE, replace=False)
        bg_scores.append(pool)
        bg_codes.append(np.full(len(pool), code, dtype=np.uint8))
    strat = {"p": float("nan"), "auc": float("nan"), "n_used": 0}
    if n_deposit_px >= 3 and bg_scores:
        strat = stratified_mannwhitney(deposit_scores, dep_codes,
                                       np.concatenate(bg_scores), np.concatenate(bg_codes))
    return n_deposit_px, strat, float(p_geo), auc_geo


def _r(x, nd):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), nd)


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
        if footprint is not None and len(footprint):
            fp_geom = make_valid(unary_union(footprint.geometry.values))
        else:
            fp_geom = gpd.GeoSeries([box(*raster_bbox)], crs="EPSG:4326").to_crs(zones.crs).iloc[0]

        n_crit = int(ss_row["n_crit"])
        hits_crit = int(ss_row["hits_crit"])
        p_binomial = float(ss_row["p_binomial"])

        strata = cover_strata(site_id, fp_geom, zones.crs, bbox_wgs84=raster_bbox)
        deposits = load_critical_deposits(paths, raster_bbox, footprint, zones.crs)
        labels = assign_cover_stratum(deposits, strata) if len(deposits) else pd.Series(dtype=object)
        n_per = {s: int((labels == s).sum()) for s in COVER_STRATA}
        if sum(n_per.values()) != n_crit:
            print(f"[WARN n_strata={sum(n_per.values())} vs n_crit={n_crit}]", end=" ")

        # --- 1. Monte Carlo: stratified (primary) and footprint-wide (geometric) ---
        p_perm = run_permutation_stratified(
            zones, raster_bbox, strata, n_per, hits_crit, n_iter=N_ITER, seed=SEED,
        )
        p_perm_geo = run_permutation(
            zones, raster_bbox, n_crit, hits_crit, n_iter=N_ITER, seed=SEED, footprint=footprint,
        )
        p_perm_label = f"<{1/N_ITER:.0e}" if p_perm == 0.0 else f"{p_perm:.4f}"

        # --- 2. Continuous-score Mann-Whitney: stratified (primary) and geometric ---
        cs = load_combined_score_raster(paths.aster_dir, site_id, site, raster_bbox)
        n_deposit_px = 0
        strat = {"p": float("nan"), "auc": float("nan"), "n_used": 0}
        p_geo, auc_geo = float("nan"), float("nan")
        if cs is not None and len(deposits):
            score, valid, transform, crs = cs
            deps = deposits if deposits.crs == crs else deposits.to_crs(crs)
            strata_r = strata if zones.crs == crs else {
                k: gpd.GeoSeries([g], crs=zones.crs).to_crs(crs).iloc[0] for k, g in strata.items()
            }
            rng = np.random.default_rng(SEED)
            n_deposit_px, strat, p_geo, auc_geo = continuous_score_tests(
                score, valid, transform, deps, strata_r, rng
            )

        rows.append({
            "site_id": site_id,
            "n_crit": n_crit,
            "hits_crit": hits_crit,
            "n_crit_bedrock": n_per["bedrock"],
            "n_crit_cover": n_per["cover"],
            "n_crit_unmapped": n_per["unmapped"],
            "p_binomial": p_binomial,
            "p_monte_carlo": round(p_perm, 4),
            "p_monte_carlo_label": p_perm_label,
            "p_monte_carlo_geometric": round(p_perm_geo, 4),
            "n_deposit_pixels_scored": n_deposit_px,
            "n_deposit_pixels_stratified": int(strat["n_used"]),
            "p_mannwhitney": _r(strat["p"], 8),
            "auc_continuous_score": _r(strat["auc"], 4),
            "p_mannwhitney_geometric": _r(p_geo, 8),
            "auc_continuous_score_geometric": _r(auc_geo, 4),
        })
        print(
            f"p_binom={p_binomial:.4g}  p_MC={p_perm_label} (geo {p_perm_geo:.4g})  "
            f"p_MWU={strat['p']:.4g} AUC={strat['auc']:.3f} (geo AUC={auc_geo:.3f})"
        )

    out = pd.DataFrame(rows).sort_values("p_binomial")
    out_path = REPO_ROOT / "results" / "phase3_monte_carlo_and_continuous_score.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
