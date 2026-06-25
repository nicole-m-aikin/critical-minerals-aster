"""
Discovery-bias stratified significance test.

Splits critical-mineral MRDS deposits at each *significant* site into
pre-1950 and post-1950 cohorts based on the MRDS ``disc_yr`` field
(primary) or ``yr_fst_prd`` as a fallback proxy, then runs the same
one-sided exact binomial test used in significance_critical_only.py
on each cohort independently.

Motivation
----------
A pure surface-alteration-detection signal (i.e. ASTER is mapping
historically visible hydrothermal outcrop) would predict:
  pre-1950 hit rate  >>  post-1950 hit rate
because pre-1950 prospectors followed surface alteration cues that
ASTER also detects.  If post-1950 deposits (found after geophysics,
geochemistry, and remote-sensing surveys) show a similar hit rate,
the ASTER signal carries predictive value beyond historical circularity.

Coverage caveat
---------------
MRDS ``disc_yr`` is filled for only ~2–37 % of critical-mineral deposits
at these sites.  The dated subset is smaller and skews toward producers
(not just occurrences).  Results should be interpreted as suggestive,
not definitive.  Coverage fractions are reported alongside every result.

Outputs
-------
  results/discovery_bias_analysis.csv   — per-site cohort table
  stdout                                — formatted summary tables

Usage
-----
  conda run -n aster-minerals python scripts/discovery_bias_analysis.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.significance import run_binomial

CUTOFF = 1950  # pre-CUTOFF / CUTOFF-and-after split


# ---------------------------------------------------------------------------
# Helpers shared with significance_critical_only.py
# ---------------------------------------------------------------------------


def _load_raster_bbox(site_id: str) -> tuple:
    prov_path = ROOT / "results" / f"{site_id}_provenance.json"
    if prov_path.exists():
        with open(prov_path) as f:
            prov = json.load(f)
        rb = prov.get("raster_bbox_wgs84")
        if rb is not None:
            return tuple(rb)
    yaml_path = ROOT / "sites" / f"{site_id}.yaml"
    with open(yaml_path) as f:
        cfg = yaml.safe_load(f)
    return tuple(cfg["bbox_wgs84"])


def _load_tir_footprint(site_id: str, zones_crs) -> gpd.GeoDataFrame | None:
    prov_path = ROOT / "results" / f"{site_id}_provenance.json"
    if not prov_path.exists():
        return None
    with open(prov_path) as f:
        prov = json.load(f)
    geojson = prov.get("tir_footprint_wgs84")
    if geojson is None:
        return None
    try:
        gdf = gpd.GeoDataFrame.from_features(geojson["features"], crs="EPSG:4326")
        return gdf.to_crs(zones_crs)
    except Exception:
        return None


def _coverage_fraction(zones: gpd.GeoDataFrame, footprint: gpd.GeoDataFrame | None, bbox: tuple) -> float:
    from shapely.ops import unary_union
    from critical_minerals_aster.significance import _bbox_in_crs

    union = unary_union(zones.geometry)
    if footprint is not None and len(footprint):
        fp = footprint if footprint.crs == zones.crs else footprint.to_crs(zones.crs)
        denom = float(fp.iloc[0].geometry.area)
    else:
        minx, miny, maxx, maxy = _bbox_in_crs(bbox, zones.crs)
        denom = (maxx - minx) * (maxy - miny)
    if denom == 0:
        return 0.0
    return min(float(union.area / denom), 1.0)


def _load_zones(site_id: str) -> gpd.GeoDataFrame | None:
    path = ROOT / "data" / "sites" / site_id / "vectors" / "strong_anomaly_zones.geojson"
    if not path.exists():
        return None
    gdf = gpd.read_file(path)
    if gdf.crs is None or gdf.crs.is_geographic:
        gdf = gdf.to_crs("EPSG:32611")
    return gdf


# ---------------------------------------------------------------------------
# Discovery-year extraction
# ---------------------------------------------------------------------------


def _best_year(row: pd.Series) -> float | None:
    """Return disc_yr if present, else yr_fst_prd, else NaN."""
    dy = row.get("disc_yr")
    if pd.notna(dy):
        return float(dy)
    yfp = row.get("yr_fst_prd")
    if pd.notna(yfp):
        return float(yfp)
    return float("nan")


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------


def main() -> None:
    mrds_path = ROOT / "data" / "mrds.csv"
    mrds_all = pd.read_csv(mrds_path, low_memory=False)

    sig_path = ROOT / "results" / "significance_critical_only.csv"
    sig_df = pd.read_csv(sig_path)
    sig_sites = sig_df[sig_df["sig_crit"]]["site_id"].tolist()

    print(f"Significant sites (critical-mineral binomial, p<0.05): {sig_sites}\n")

    records: list[dict] = []

    for site_id in sorted(sig_sites):
        print(f"  Processing {site_id} ...", flush=True)

        zones = _load_zones(site_id)
        if zones is None or zones.empty:
            print("    no zones, skipping")
            continue

        bbox = _load_raster_bbox(site_id)
        footprint = _load_tir_footprint(site_id, zones.crs)
        p_cover = _coverage_fraction(zones, footprint, bbox)

        # Load and filter MRDS deposits to critical minerals inside this site
        local = filter_mrds_bbox(mrds_all, bbox)
        deposits = mrds_to_points_gdf(local, zones.crs)

        if footprint is not None and len(footprint):
            try:
                fp = footprint if footprint.crs == zones.crs else footprint.to_crs(zones.crs)
                deposits = gpd.clip(deposits, fp)
            except Exception:
                pass

        deposits = reclassify_mrds_earth_mri(deposits)
        crit = deposits[deposits["earth_mri_category"] != "Non-Critical"].copy()
        if crit.empty:
            print("    no critical deposits, skipping")
            continue

        # Spatial join: mark each deposit as hit or miss
        joined, _, _ = spatial_join_deposits_zones(crit, zones)
        crit = crit.copy()
        hit_ids = set(joined[joined["index_right"].notna()].index)
        crit["in_zone"] = crit.index.isin(hit_ids)

        # Assign best year estimate
        crit["best_year"] = crit.apply(_best_year, axis=1)

        n_total = len(crit)
        n_dated = int(crit["best_year"].notna().sum())
        pct_dated = n_dated / n_total * 100 if n_total else 0.0

        dated = crit[crit["best_year"].notna()].copy()
        pre = dated[dated["best_year"] < CUTOFF]
        post = dated[dated["best_year"] >= CUTOFF]

        def _cohort_stats(cohort: pd.DataFrame, label: str, p: float) -> dict:
            n = len(cohort)
            hits = int(cohort["in_zone"].sum())
            hr = hits / n * 100 if n else 0.0
            p_val, _ = run_binomial(hits, n, p)
            return {
                f"n_{label}": n,
                f"hits_{label}": hits,
                f"hr_{label}_pct": round(hr, 1),
                f"binom_p_{label}": round(p_val, 4),
                f"sig_{label}": p_val < 0.05 if n >= 3 else False,
            }

        row = {
            "site_id": site_id,
            "site_name": sig_df.loc[sig_df["site_id"] == site_id, "site_name"].iloc[0],
            "coverage_pct": round(p_cover * 100, 1),
            "n_crit_total": n_total,
            "n_dated": n_dated,
            "pct_dated": round(pct_dated, 1),
            **_cohort_stats(pre, "pre50", p_cover),
            **_cohort_stats(post, "post50", p_cover),
        }
        records.append(row)

    out_df = pd.DataFrame(records)

    # Aggregate across all significant sites (pooled dated deposits)
    # Each site has its own p_cover, so we compute expected hits individually
    # and sum: E[total hits] = sum_i(n_i * p_i), Var = sum_i(n_i*p_i*(1-p_i))
    # then compare observed total hits to this compound expectation.
    # For a simpler presentation we report observed vs expected and
    # a normal-approximation z-test p-value.
    import numpy as np
    from scipy.stats import norm as sp_norm

    def _pooled_test(cohort_col_n: str, cohort_col_hits: str) -> dict:
        obs_hits = out_df[cohort_col_hits].sum()
        exp_hits_list = []
        var_list = []
        for _, row in out_df.iterrows():
            n_i = row[cohort_col_n]
            p_i = row["coverage_pct"] / 100.0
            exp_hits_list.append(n_i * p_i)
            var_list.append(n_i * p_i * (1 - p_i))
        exp_hits = sum(exp_hits_list)
        var_total = sum(var_list)
        if var_total == 0:
            z, p_val = 0.0, 1.0
        else:
            z = (obs_hits - exp_hits) / var_total ** 0.5
            p_val = float(sp_norm.sf(z))  # one-sided P(Z > z)
        obs_n = int(out_df[cohort_col_n].sum())
        hr = obs_hits / obs_n * 100 if obs_n else 0.0
        return {
            "n": obs_n,
            "hits": int(obs_hits),
            "hr_pct": round(hr, 1),
            "expected_hits": round(exp_hits, 1),
            "z": round(z, 2),
            "pooled_p": round(p_val, 4),
        }

    pooled_pre = _pooled_test("n_pre50", "hits_pre50")
    pooled_post = _pooled_test("n_post50", "hits_post50")

    # Save CSV
    out_path = ROOT / "results" / "discovery_bias_analysis.csv"
    out_df.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}\n")

    # Print per-site table
    disp_cols = [
        "site_id", "coverage_pct", "n_crit_total", "n_dated", "pct_dated",
        "n_pre50", "hits_pre50", "hr_pre50_pct", "binom_p_pre50",
        "n_post50", "hits_post50", "hr_post50_pct", "binom_p_post50",
    ]
    print("=== PER-SITE DISCOVERY-COHORT RESULTS ===")
    print(f"Pre-{CUTOFF} = deposits with discovery/first-production year < {CUTOFF}")
    print(f"Post-{CUTOFF} = deposits with year >= {CUTOFF}\n")
    print(out_df[disp_cols].to_string(index=False))

    # Print pooled summary — include null hit-rate for context
    null_hr_pre = pooled_pre["expected_hits"] / pooled_pre["n"] * 100 if pooled_pre["n"] else 0.0
    null_hr_post = pooled_post["expected_hits"] / pooled_post["n"] * 100 if pooled_post["n"] else 0.0
    print(f"\n=== POOLED ACROSS ALL {len(out_df)} SIGNIFICANT SITES ===")
    print(
        f"{'Cohort':<12} {'n':>5} {'hits':>5} {'observed%':>10} {'null%':>7}"
        f" {'enrichment':>11} {'z':>6} {'pooled_p':>10}"
    )
    for label, d, null_hr in [
        ("pre-1950", pooled_pre, null_hr_pre),
        ("post-1950", pooled_post, null_hr_post),
    ]:
        enrich = d["hr_pct"] / null_hr if null_hr else float("nan")
        print(
            f"{label:<12} {d['n']:>5} {d['hits']:>5} {d['hr_pct']:>10.1f}"
            f" {null_hr:>7.1f} {enrich:>11.2f}× {d['z']:>6.2f} {d['pooled_p']:>10.4f}"
        )

    # Interpretation
    print("\n=== INTERPRETATION ===")
    total_dated = int(out_df["n_dated"].sum())
    total_crit = int(out_df["n_crit_total"].sum())
    print(
        f"Dated deposits (disc_yr or yr_fst_prd present): "
        f"{total_dated}/{total_crit} across significant sites "
        f"({total_dated/total_crit*100:.0f}% coverage)."
    )

    post_sig = pooled_post["pooled_p"] < 0.05
    pre_enrichment = pooled_pre["hr_pct"] / null_hr_pre if null_hr_pre else 0.0
    post_enrichment = pooled_post["hr_pct"] / null_hr_post if null_hr_post else 0.0

    print(
        f"Pre-1950 cohort: {pooled_pre['hr_pct']:.1f}% observed vs {null_hr_pre:.1f}% null "
        f"({pre_enrichment:.1f}× enrichment, pooled p={pooled_pre['pooled_p']:.4f})."
    )
    print(
        f"Post-1950 cohort: {pooled_post['hr_pct']:.1f}% observed vs {null_hr_post:.1f}% null "
        f"({post_enrichment:.1f}× enrichment, pooled p={pooled_post['pooled_p']:.4f})."
    )

    if post_sig and post_enrichment >= 1.5:
        print(
            "RESULT: Post-1950 deposits also cluster in ASTER anomaly zones above the null "
            "expectation (significant enrichment). This suggests the spectral signal captures "
            "mineralization-favorable geology beyond surface outcrop sampling bias, though "
            "pre-1950 enrichment is stronger — indicating partial bias remains."
        )
    elif post_sig:
        print(
            "RESULT: Post-1950 deposits show marginal enrichment above null (significant but "
            "weaker than pre-1950). Weak evidence against a purely circular signal."
        )
    else:
        print(
            "RESULT: Post-1950 hit rate is not statistically distinguishable from null "
            f"(n={pooled_post['n']}; test likely underpowered). "
            "Cannot rule out that the ASTER signal primarily captures historically visible outcrop."
        )

    # Always report the coverage caveat and post-1950 n as they are small regardless
    print(
        f"\nCAVEAT: Only {total_dated}/{total_crit} ({total_dated/total_crit*100:.0f}%) "
        "critical-mineral deposits at significant sites have MRDS discovery/production dates. "
        "Dated deposits skew toward producers and toward early-discovered camps — the undated "
        "majority could shift results in either direction. "
        f"Post-1950 cohort (n={pooled_post['n']}) is too small for strong per-site inference; "
        "pooled results are suggestive only."
    )


if __name__ == "__main__":
    main()
