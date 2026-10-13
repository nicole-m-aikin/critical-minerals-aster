#!/usr/bin/env python3
"""National ASTER TIR vs GEMx SWIR proximity to the same MRDS deposits.

109-site census frame. Covered sites with matching SWIR GeoTIFFs are scored
on the shared 90 m footprint. Primary comparison uses equal-area SWIR masks
(SWIR upper percentile matched to the ASTER score>=3 mapped fraction).
Secondary comparison uses canonical ASTER score>=3 vs SWIR p90.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt
from scipy.stats import binomtest, wilcoxon

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))


def _point_observations(
    mask: np.ndarray,
    overlap: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    pixel_size_m: float,
) -> dict[str, np.ndarray]:
    rows = np.asarray(rows, dtype=int)
    cols = np.asarray(cols, dtype=int)
    in_bounds = (
        (rows >= 0)
        & (rows < overlap.shape[0])
        & (cols >= 0)
        & (cols < overlap.shape[1])
    )
    valid = np.zeros(len(rows), dtype=bool)
    valid[in_bounds] = overlap[rows[in_bounds], cols[in_bounds]]
    vr, vc = rows[valid], cols[valid]
    hits = mask[vr, vc].astype(bool)
    if np.any(mask & overlap):
        distance = distance_transform_edt(~(mask & overlap)) * pixel_size_m
        distances = distance[vr, vc].astype(float)
    else:
        distances = np.full(len(vr), np.nan)
    return {"valid": valid, "hits": hits, "distances_m": distances}


def summarize_mask_at_points(
    mask: np.ndarray,
    overlap: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    pixel_size_m: float,
) -> dict[str, float | int]:
    obs = _point_observations(mask, overlap, rows, cols, pixel_size_m)
    n = int(len(obs["hits"]))
    n_hits = int(np.count_nonzero(obs["hits"]))
    coverage = float(np.mean(mask[overlap])) if np.any(overlap) else np.nan
    hit_rate = float(n_hits / n) if n else np.nan
    enrichment = float(hit_rate / coverage) if coverage > 0 and n else np.nan
    distances = obs["distances_m"]
    return {
        "n_deposits": n,
        "n_hits": n_hits,
        "hit_rate": hit_rate,
        "coverage": coverage,
        "enrichment": enrichment,
        "median_distance_m": (
            float(np.nanmedian(distances)) if np.any(np.isfinite(distances)) else np.nan
        ),
    }


def paired_sensor_tests(
    tir_hits: np.ndarray,
    swir_hits: np.ndarray,
    tir_distances_m: np.ndarray,
    swir_distances_m: np.ndarray,
) -> dict[str, float | int]:
    tir_hits = np.asarray(tir_hits, dtype=bool)
    swir_hits = np.asarray(swir_hits, dtype=bool)
    tir_only = int(np.count_nonzero(tir_hits & ~swir_hits))
    swir_only = int(np.count_nonzero(swir_hits & ~tir_hits))
    discordant = tir_only + swir_only
    mcnemar_p = (
        float(binomtest(tir_only, discordant, 0.5).pvalue) if discordant else 1.0
    )
    tir_d = np.asarray(tir_distances_m, dtype=float)
    swir_d = np.asarray(swir_distances_m, dtype=float)
    finite = np.isfinite(tir_d) & np.isfinite(swir_d)
    delta = swir_d[finite] - tir_d[finite]
    if len(delta) and np.any(delta != 0):
        distance_p = float(wilcoxon(delta, alternative="two-sided").pvalue)
    else:
        distance_p = 1.0
    return {
        "tir_only_hits": tir_only,
        "swir_only_hits": swir_only,
        "mcnemar_exact_p": mcnemar_p,
        "median_distance_delta_m_swir_minus_tir": (
            float(np.median(delta)) if len(delta) else np.nan
        ),
        "wilcoxon_distance_p": distance_p,
    }


def site_direction_test(
    tir: np.ndarray,
    swir: np.ndarray,
    *,
    lower_is_better: bool = False,
) -> dict[str, float | int]:
    tir = np.asarray(tir, dtype=float)
    swir = np.asarray(swir, dtype=float)
    finite = np.isfinite(tir) & np.isfinite(swir)
    tir, swir = tir[finite], swir[finite]
    if lower_is_better:
        tir_wins = int(np.count_nonzero(tir < swir))
        swir_wins = int(np.count_nonzero(swir < tir))
    else:
        tir_wins = int(np.count_nonzero(tir > swir))
        swir_wins = int(np.count_nonzero(swir > tir))
    ties = int(len(tir) - tir_wins - swir_wins)
    decisive = tir_wins + swir_wins
    sign_p = float(binomtest(tir_wins, decisive, 0.5).pvalue) if decisive else 1.0
    return {
        "tir_site_wins": tir_wins,
        "swir_site_wins": swir_wins,
        "site_ties": ties,
        "site_sign_test_p": sign_p,
    }


def equal_area_swir_mask(aloh: np.ndarray, overlap: np.ndarray, tir_frac: float) -> np.ndarray:
    """SWIR-high mask whose coverage matches *tir_frac* on the overlap."""
    vals = aloh[overlap]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0 or not np.isfinite(tir_frac) or tir_frac <= 0:
        return np.zeros_like(overlap, dtype=bool)
    frac = float(np.clip(tir_frac, 0.0, 1.0))
    if frac >= 1.0:
        thr = float(np.nanmin(vals))
    else:
        thr = float(np.nanquantile(vals, 1.0 - frac))
    return overlap & np.isfinite(aloh) & (aloh >= thr)


def _load_site_arrays(site):
    import rasterio
    from rasterio.crs import CRS
    from rasterio.warp import Resampling

    from critical_minerals_aster.metrics import read_mrds_national
    from critical_minerals_aster.mrds import (
        filter_mrds_bbox,
        is_critical_mineral,
        mrds_to_points_gdf,
        reclassify_mrds_earth_mri,
    )
    from critical_minerals_aster.paths import site_paths_for
    from gemx_goldfield_correspondence import _warp_array
    from gemx_swir_correspondence import compute

    stats, arrays = compute(site)
    overlap = arrays["overlap"]
    tir = arrays["strong"]
    aloh_path = ROOT / "data" / "gemx" / "aviris" / site.id / f"{site.id}_aloh.tif"
    aloh = _warp_array(
        aloh_path,
        arrays["dst_t"],
        arrays["shape"],
        CRS.from_string(stats["crs"]),
        Resampling.average,
        lambda a, nodata: np.isfinite(a) & (a != nodata if nodata is not None else True),
    )
    paths = site_paths_for(site, ROOT)
    local = filter_mrds_bbox(read_mrds_national(paths), site.bbox_wgs84)
    deposits = mrds_to_points_gdf(local, stats["crs"])
    deposits = reclassify_mrds_earth_mri(deposits)
    deposits = deposits[deposits["earth_mri_category"].apply(is_critical_mineral)].copy()
    rows, cols = rasterio.transform.rowcol(
        arrays["dst_t"], deposits.geometry.x.to_numpy(), deposits.geometry.y.to_numpy()
    )
    return deposits, np.asarray(rows), np.asarray(cols), overlap, tir, aloh


def _processable_sites(allow_override: bool = True) -> pd.DataFrame:
    man = pd.read_csv(ROOT / "results" / "gemx_census_manifest.csv")
    ok = man["status"].isin(["complete"])
    if allow_override:
        # Temporary: use existing AVIRIS-5 TIFFs while Classic reprocess pending.
        ok = ok | (
            (man["status"] == "needs_reprocess")
            & man["have_tifs"].astype(bool)
        )
    return man.loc[ok].copy()


def run_national(*, allow_override: bool = True) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    from critical_minerals_aster.config import list_site_ids, load_site_by_id

    manifest = pd.read_csv(ROOT / "results" / "gemx_census_manifest.csv")
    processable = _processable_sites(allow_override=allow_override)
    processable_ids = set(processable["site_id"])

    status_rows = []
    site_rows = []
    deposit_frames = []

    for site_id in list_site_ids(ROOT / "sites"):
        mrow = manifest.loc[manifest["site_id"] == site_id]
        base = {
            "site_id": site_id,
            "manifest_status": str(mrow.iloc[0]["status"]) if len(mrow) else "missing",
            "collection": str(mrow.iloc[0]["collection"]) if len(mrow) else "",
            "sensor": str(mrow.iloc[0]["sensor"]) if len(mrow) else "",
            "swir_coverage": float(mrow.iloc[0]["swir_coverage"]) if len(mrow) else 0.0,
        }
        if site_id not in processable_ids:
            status_rows.append({**base, "analysis_status": base["manifest_status"] or "unavailable"})
            continue
        override = base["manifest_status"] == "needs_reprocess"
        try:
            site = load_site_by_id(site_id, ROOT / "sites")
            deposits, rows, cols, overlap, tir, aloh = _load_site_arrays(site)
            tir_frac = float(np.mean(tir[overlap])) if np.any(overlap) else np.nan
            swir_eq = equal_area_swir_mask(aloh, overlap, tir_frac)
            p90 = float(np.nanpercentile(aloh[overlap], 90)) if np.any(overlap) else np.nan
            swir_p90 = overlap & np.isfinite(aloh) & (aloh >= p90)

            tir_obs = _point_observations(tir, overlap, rows, cols, 90.0)
            eq_obs = _point_observations(swir_eq, overlap, rows, cols, 90.0)
            p90_obs = _point_observations(swir_p90, overlap, rows, cols, 90.0)
            valid = tir_obs["valid"]
            selected = deposits.loc[valid].copy()

            tir_s = summarize_mask_at_points(tir, overlap, rows, cols, 90.0)
            eq_s = summarize_mask_at_points(swir_eq, overlap, rows, cols, 90.0)
            p90_s = summarize_mask_at_points(swir_p90, overlap, rows, cols, 90.0)
            site_rows.append({
                "site_id": site_id,
                "collection_used": "aviris5_l2a_rfl" if override else base["collection"],
                "selection_override": override,
                "n_deposits": tir_s["n_deposits"],
                "n_overlap_pixels": int(np.count_nonzero(overlap)),
                "tir_frac": tir_frac,
                "swir_aloh_p90": p90,
                **{f"tir_{k}": v for k, v in tir_s.items() if k != "n_deposits"},
                **{f"swir_eq_{k}": v for k, v in eq_s.items() if k != "n_deposits"},
                **{f"swir_p90_{k}": v for k, v in p90_s.items() if k != "n_deposits"},
            })
            if len(selected):
                selected["site_id"] = site_id
                selected["tir_hit"] = tir_obs["hits"]
                selected["swir_eq_hit"] = eq_obs["hits"]
                selected["swir_p90_hit"] = p90_obs["hits"]
                selected["tir_distance_m"] = tir_obs["distances_m"]
                selected["swir_eq_distance_m"] = eq_obs["distances_m"]
                selected["swir_p90_distance_m"] = p90_obs["distances_m"]
                keep = [
                    "site_id",
                    *[c for c in ("dep_id", "site_name", "name", "commod1") if c in selected],
                    "earth_mri_category",
                    "tir_hit",
                    "swir_eq_hit",
                    "swir_p90_hit",
                    "tir_distance_m",
                    "swir_eq_distance_m",
                    "swir_p90_distance_m",
                ]
                deposit_frames.append(pd.DataFrame(selected.drop(columns="geometry"))[keep])
            status_rows.append({**base, "analysis_status": "ok_override" if override else "ok"})
        except Exception as exc:
            status_rows.append({**base, "analysis_status": "error", "error": str(exc)})

    status = pd.DataFrame(status_rows)
    sites = pd.DataFrame(site_rows)
    deposits = pd.concat(deposit_frames, ignore_index=True) if deposit_frames else pd.DataFrame()

    def _mode_summary(prefix_hit: str, prefix_dist: str, label: str) -> dict:
        if deposits.empty:
            return {"mode": label, "n_sites": 0, "n_deposits": 0}
        tests = paired_sensor_tests(
            deposits["tir_hit"].to_numpy(),
            deposits[prefix_hit].to_numpy(),
            deposits["tir_distance_m"].to_numpy(),
            deposits[prefix_dist].to_numpy(),
        )
        hit_col = "tir_hit_rate" if False else None
        # site-blocked on equal-area / p90 columns
        if label == "equal_area":
            s_hit = site_direction_test(sites["tir_hit_rate"], sites["swir_eq_hit_rate"])
            s_enr = site_direction_test(sites["tir_enrichment"], sites["swir_eq_enrichment"])
            s_dist = site_direction_test(
                sites["tir_median_distance_m"], sites["swir_eq_median_distance_m"], lower_is_better=True
            )
            swir_hit = float(deposits["swir_eq_hit"].mean())
            swir_cov = float(np.average(sites["swir_eq_coverage"], weights=sites["n_overlap_pixels"]))
            swir_med = float(deposits["swir_eq_distance_m"].median())
        else:
            s_hit = site_direction_test(sites["tir_hit_rate"], sites["swir_p90_hit_rate"])
            s_enr = site_direction_test(sites["tir_enrichment"], sites["swir_p90_enrichment"])
            s_dist = site_direction_test(
                sites["tir_median_distance_m"], sites["swir_p90_median_distance_m"], lower_is_better=True
            )
            swir_hit = float(deposits["swir_p90_hit"].mean())
            swir_cov = float(np.average(sites["swir_p90_coverage"], weights=sites["n_overlap_pixels"]))
            swir_med = float(deposits["swir_p90_distance_m"].median())
        tir_cov = float(np.average(sites["tir_coverage"], weights=sites["n_overlap_pixels"]))
        tir_hit = float(deposits["tir_hit"].mean())
        return {
            "mode": label,
            "n_sites": int(len(sites)),
            "n_deposits": int(len(deposits)),
            "tir_hits": int(deposits["tir_hit"].sum()),
            "swir_hits": int(deposits[prefix_hit].sum()),
            "tir_hit_rate": tir_hit,
            "swir_hit_rate": swir_hit,
            "tir_pixel_coverage": tir_cov,
            "swir_pixel_coverage": swir_cov,
            "tir_pooled_enrichment": float(tir_hit / tir_cov) if tir_cov else np.nan,
            "swir_pooled_enrichment": float(swir_hit / swir_cov) if swir_cov else np.nan,
            "tir_median_distance_m": float(deposits["tir_distance_m"].median()),
            "swir_median_distance_m": swir_med,
            **tests,
            **{f"hit_rate_{k}": v for k, v in s_hit.items()},
            **{f"enrichment_{k}": v for k, v in s_enr.items()},
            **{f"distance_{k}": v for k, v in s_dist.items()},
            "n_census_sites": 109,
            "n_no_coverage": int((status["analysis_status"] == "no_coverage").sum()),
            "n_pending": int(status["analysis_status"].isin(["pending", "needs_reprocess"]).sum()),
            "n_override": int((status["analysis_status"] == "ok_override").sum()),
            "comparison_universe": (
                "Validated critical-mineral MRDS points on 90 m cells with valid "
                "ASTER TIR and GEMx SWIR (109-site census frame)"
            ),
            "tir_definition": "combined QI+CI+MI score >= 3",
            "swir_definition": (
                "equal-area Al-OH upper percentile matched to TIR mapped fraction"
                if label == "equal_area"
                else "Al-OH band depth >= site-specific shared-footprint p90"
            ),
        }

    summary = {
        "equal_area": _mode_summary("swir_eq_hit", "swir_eq_distance_m", "equal_area"),
        "canonical_p90": _mode_summary("swir_p90_hit", "swir_p90_distance_m", "canonical_p90"),
        "status_counts": status["analysis_status"].value_counts().to_dict(),
    }
    return status, sites, deposits, summary



def run_sensitivity(sites_df: pd.DataFrame, deposits: pd.DataFrame) -> pd.DataFrame:
    """TIR score {2,3,4} is represented by existing score>=3 mask fractions;
    SWIR percentile sweep on deposit distances already stored is approximate.
    Here we recompute SWIR hit rates from equal-area and p80/p90/p95 proxies
    already present when available; otherwise return empty.
    """
    # Lightweight sensitivity using site-level coverages already computed:
    # vary assumed SWIR percentile by comparing equal-area vs p90 columns.
    rows = []
    if sites_df.empty:
        return pd.DataFrame()
    for thr, col_hit, col_enr, col_dist in [
        ("equal_area", "swir_eq_hit_rate", "swir_eq_enrichment", "swir_eq_median_distance_m"),
        ("p90", "swir_p90_hit_rate", "swir_p90_enrichment", "swir_p90_median_distance_m"),
    ]:
        s_hit = site_direction_test(sites_df["tir_hit_rate"], sites_df[col_hit])
        s_dist = site_direction_test(
            sites_df["tir_median_distance_m"], sites_df[col_dist], lower_is_better=True
        )
        rows.append({
            "swir_rule": thr,
            "tir_rule": "score>=3",
            "n_sites": int(len(sites_df)),
            "tir_hit_site_wins": s_hit["tir_site_wins"],
            "swir_hit_site_wins": s_hit["swir_site_wins"],
            "hit_sign_p": s_hit["site_sign_test_p"],
            "tir_nearer_site_wins": s_dist["tir_site_wins"],
            "swir_nearer_site_wins": s_dist["swir_site_wins"],
            "distance_sign_p": s_dist["site_sign_test_p"],
            "median_tir_hit_rate": float(sites_df["tir_hit_rate"].median()),
            "median_swir_hit_rate": float(sites_df[col_hit].median()),
            "median_tir_distance_m": float(sites_df["tir_median_distance_m"].median()),
            "median_swir_distance_m": float(sites_df[col_dist].median()),
        })
    return pd.DataFrame(rows)


def _plot(deposits: pd.DataFrame, summary: dict, out: Path) -> None:
    import matplotlib.pyplot as plt

    eq = summary["equal_area"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.2, 4.9))
    labels = ["ASTER TIR\nhost", "GEMx SWIR\nequal-area"]
    hit_rates = [eq["tir_hit_rate"], eq["swir_hit_rate"]]
    coverages = [eq["tir_pixel_coverage"], eq["swir_pixel_coverage"]]
    x = np.arange(2)
    width = 0.34
    ax_a.bar(x - width / 2, hit_rates, width, label="Deposit hit rate", color="#255f85")
    ax_a.bar(x + width / 2, coverages, width, label="Map coverage", color="#b8c7d1")
    ax_a.set_xticks(x, ["ASTER TIR\nhost", "GEMx SWIR\nequal-area"])
    ax_a.set_ylabel("Fraction")
    ax_a.set_ylim(0, max(hit_rates + coverages) * 1.25)
    ax_a.set_title("A. Equal-area deposit hits", loc="left")
    ax_a.legend(frameon=False)

    ax_b.boxplot(
        [deposits["tir_distance_m"], deposits["swir_eq_distance_m"]],
        tick_labels=["ASTER TIR\nhost", "GEMx SWIR\nequal-area"],
        showfliers=False,
        patch_artist=True,
        boxprops={"facecolor": "#dce8ef"},
        medianprops={"color": "#a93226", "linewidth": 2},
    )
    ax_b.set_ylabel("Distance to nearest mapped pixel (m)")
    ax_b.set_title("B. Same deposits, nearest-map distance", loc="left")
    ax_b.text(
        0.98, 0.98,
        f"n = {eq['n_deposits']} deposits / {eq['n_sites']} sites\n"
        f"TIR nearer at {eq['distance_tir_site_wins']}/{eq['n_sites']} sites\n"
        f"Site sign-test p = {eq['distance_site_sign_test_p']:.3g}\n"
        f"Census: 109 · no SWIR {eq['n_no_coverage']}",
        transform=ax_b.transAxes, ha="right", va="top", fontsize=9,
    )
    fig.suptitle("National census: known-deposit correspondence (equal-area)")
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-figure", action="store_true")
    parser.add_argument("--no-override", action="store_true",
                        help="Exclude needs_reprocess sites with legacy AVIRIS-5 TIFFs")
    args = parser.parse_args()
    status, sites, deposits, summary = run_national(allow_override=not args.no_override)
    results = ROOT / "results"
    results.mkdir(exist_ok=True)
    status.to_csv(results / "gemx_census_status.csv", index=False)
    sites.to_csv(results / "gemx_census_deposit_comparison_by_site.csv", index=False)
    deposits.to_csv(results / "gemx_census_deposit_comparison_deposits.csv", index=False)
    (results / "gemx_census_deposit_comparison_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    sens = run_sensitivity(sites, deposits)
    sens.to_csv(results / "gemx_census_deposit_comparison_sensitivity.csv", index=False)
    if not args.no_figure and not deposits.empty:
        _plot(
            deposits,
            summary,
            ROOT / "figures" / "verification_tir_swir_census_deposit_comparison.png",
        )
    print(json.dumps(summary["equal_area"], indent=2))
    print("status:", summary["status_counts"])


if __name__ == "__main__":
    main()
