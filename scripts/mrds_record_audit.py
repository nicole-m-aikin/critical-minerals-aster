#!/usr/bin/env python
"""
Verification B — independent audit of legacy MRDS showings against USMIN.

CIM (2018) Exploration Best Practice §2.1.1.1 asks that legacy data be
verified by an independent check of surface showings.  Field visits are out
of scope, so the desk analogue here is USGS USMIN: mine and prospect symbols
digitised from historical 7.5'/15' topographic maps, a source with no shared
lineage with MRDS report coordinates.  An MRDS record is "USMIN-verified" at
radius R if a USMIN point (shaft, adit, prospect pit, quarry ...) or polygon
(tailings, dump, open pit, disturbed surface) lies within R metres.

Two analyses per FDR-significant site:

1. Stratified audit sample -- up to 20 critical-mineral MRDS records inside
   the TIR footprint, drawn with ``numpy.random.default_rng(42)`` (fresh per
   site, records sorted by dep_id first) from four strata: score grade
   (A-C vs D-E/blank) x zone hit/miss.  Allocation is equal across strata,
   with shortfalls redistributed; ``weight`` = stratum size / n drawn.
2. Census -- every critical-mineral record: % within 250 / 500 / 1000 m of a
   USMIN feature, and the hit test re-run on the USMIN-corroborated (<=250 m)
   subset with the unchanged site-specific null p0.

Key question: do hits and misses verify at different rates?  A hits-only
positional bias (e.g. hits better located because they are producers on
mapped workings) would bias hit rates.  Tested with Fisher's exact 2x2
(hit/miss x verified/not) pooled across sites, plus a Cochran-Mantel-Haenszel
test stratified by site as a guard against Simpson's-paradox pooling.

Where USMIN returned no features for a site, records are marked
``unverifiable_no_reference`` rather than "not verified".

Outputs
-------
    results/mrds_audit.csv           -- one row per sampled record
    results/mrds_audit_summary.csv   -- one row per site + an ALL_SITES pooled row

Usage
-----
    conda run -n aster-minerals python scripts/mrds_record_audit.py
    conda run -n aster-minerals python scripts/mrds_record_audit.py --sample-size 30
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from shapely.strtree import STRtree
from statsmodels.stats.contingency_tables import StratifiedTable

from critical_minerals_aster.significance import run_binomial, significant_site_sets
from critical_minerals_aster.usmin import (
    SiteInputs,
    concat_gdfs,
    load_site_inputs,
    load_usmin,
    nearest_feature_distance,
    quad_coverage,
)

FDR_SITES, ROBUST_SITES = significant_site_sets(REPO_ROOT)

RADII_M = (250, 500, 1000)
CACHE_DIR = REPO_ROOT / "data" / "usmin"
OUT_RECORDS = REPO_ROOT / "results" / "mrds_audit.csv"
OUT_SUMMARY = REPO_ROOT / "results" / "mrds_audit_summary.csv"
SEED = 42

RECORD_COLS = [
    "dep_id", "site_name", "latitude", "longitude", "commod1", "earth_mri_category",
    "score", "dev_stat", "qc_shared_coords", "qc_score_high", "qc_producer",
]


def grade_group(score) -> str:
    s = str(score).strip().upper() if pd.notna(score) else ""
    return "A-C" if s in ("A", "B", "C") else "D-E"


def allocate(sizes: dict[str, int], total: int) -> dict[str, int]:
    """Equal allocation across strata, redistributing shortfalls round-robin."""
    alloc = {k: 0 for k in sizes}
    remaining = min(total, sum(sizes.values()))
    while remaining > 0:
        open_strata = [k for k in sorted(sizes) if alloc[k] < sizes[k]]
        for k in open_strata:
            if remaining == 0:
                break
            alloc[k] += 1
            remaining -= 1
    return alloc


def stratified_sample(dep: pd.DataFrame, n: int, seed: int = SEED) -> pd.DataFrame:
    dep = dep.sort_values("dep_id", kind="stable")
    strata = {k: g for k, g in dep.groupby("stratum")}
    alloc = allocate({k: len(g) for k, g in strata.items()}, n)
    rng = np.random.default_rng(seed)
    picks = []
    for k in sorted(strata):
        g, m = strata[k], alloc[k]
        if m == 0:
            continue
        idx = rng.choice(len(g), size=m, replace=False)
        s = g.iloc[np.sort(idx)].copy()
        s["stratum_N"] = len(g)
        s["stratum_n"] = m
        s["weight"] = len(g) / m
        picks.append(s)
    return pd.concat(picks) if picks else dep.iloc[:0]


def pct(mask: pd.Series) -> float:
    return round(float(mask.mean()) * 100, 1) if len(mask) else np.nan


def fisher_hit_vs_verified(df: pd.DataFrame, col: str) -> tuple[list[list[int]], float, float]:
    """2x2 [[hit&ver, hit&not], [miss&ver, miss&not]] -> (table, odds ratio, p)."""
    hit = df["inside_zone"].astype(bool)
    ver = df[col].astype(bool)
    t = [[int((hit & ver).sum()), int((hit & ~ver).sum())],
         [int((~hit & ver).sum()), int((~hit & ~ver).sum())]]
    if min(sum(t[0]), sum(t[1])) == 0:
        return t, np.nan, np.nan
    odds, p = fisher_exact(t)
    return t, float(odds), float(p)


def cmh_by_site(df: pd.DataFrame, col: str) -> tuple[float, float]:
    """Mantel-Haenszel pooled OR and CMH p across per-site 2x2 tables."""
    tables = []
    for _, g in df.groupby("site_id"):
        t, _, _ = fisher_hit_vs_verified(g, col)
        arr = np.asarray(t, dtype=float)
        if arr.sum(axis=1).min() > 0 and arr.sum(axis=0).min() > 0:
            tables.append(arr)
    if len(tables) < 2:
        return np.nan, np.nan
    st = StratifiedTable(tables)
    return float(st.oddsratio_pooled), float(st.test_null_odds(correction=True).pvalue)


def site_records(si: SiteInputs, refresh: bool) -> tuple[pd.DataFrame, dict]:
    """All critical records with USMIN distances, plus site coverage metadata."""
    dep = si.deposits.copy()
    coverage = {
        "usmin_status": si.usmin_status,
        "n_usmin_points_fp": 0 if si.usmin_points is None else len(si.usmin_points),
        "n_usmin_polygons_fp": 0 if si.usmin_polygons is None else len(si.usmin_polygons),
    }
    ref = None
    if si.usmin_status != "unavailable":
        # Unclipped bbox layers so records near the footprint edge can match
        # workings just outside it.
        pts, polys = load_usmin(si.raster_bbox, CACHE_DIR, refresh)
        pts["layer"], polys["layer"] = "point", "polygon"
        ref = concat_gdfs([pts, polys]).to_crs(si.crs)
        n_q, n_qf = quad_coverage(si.footprint_wgs84, concat_gdfs([pts, polys]).geometry)
        coverage.update({"n_quads_footprint": n_q, "n_quads_with_usmin": n_qf,
                         "n_usmin_topo_sheets": int(pd.concat([pts.get("topo_name", pd.Series(dtype=str)),
                                                               polys.get("topo_name", pd.Series(dtype=str))])
                                                    .replace("", np.nan).dropna().nunique())})

    has_ref = ref is not None and len(ref) > 0
    coverage["verifiable"] = has_ref
    if has_ref:
        dist = nearest_feature_distance(dep.geometry, ref.geometry)
        tree = STRtree(list(ref.geometry))
        nn = tree.query_nearest(list(dep.geometry), all_matches=False)
        near_type = np.array([""] * len(dep), dtype=object)
        near_layer = np.array([""] * len(dep), dtype=object)
        near_type[nn[0]] = ref["ftr_type"].values[nn[1]]
        near_layer[nn[0]] = ref["layer"].values[nn[1]]
    else:
        dist = np.full(len(dep), np.nan)
        near_type = near_layer = np.array([""] * len(dep), dtype=object)
    dep["dist_usmin_m"] = np.round(dist, 1)
    dep["nearest_usmin_type"] = near_type
    dep["nearest_usmin_layer"] = near_layer
    for r in RADII_M:
        dep[f"verified_{r}m"] = (dep["dist_usmin_m"] <= r) if has_ref else pd.NA
    if has_ref:
        dep["verification_status"] = np.select(
            [dep["dist_usmin_m"] <= 250, dep["dist_usmin_m"] <= 500, dep["dist_usmin_m"] <= 1000],
            ["verified_250m", "verified_500m", "verified_1000m"], default="not_verified",
        )
    else:
        dep["verification_status"] = "unverifiable_no_reference"
    dep["grade_group"] = dep["score"].apply(grade_group)
    dep["hit_status"] = np.where(dep["inside_zone"], "hit", "miss")
    dep["stratum"] = dep["grade_group"] + "|" + dep["hit_status"]
    return dep, coverage


def summarise_site(si: SiteInputs, dep: pd.DataFrame, samp: pd.DataFrame, coverage: dict) -> dict:
    row = {"site_id": si.site_id, "site_name": si.site_name, "robust": si.site_id in ROBUST_SITES,
           **coverage}
    n, h = len(dep), int(dep["inside_zone"].sum())
    p_base, _ = run_binomial(h, n, si.p0)
    row.update({"p0_null": round(si.p0, 4), "n_crit": n, "hits_crit": h,
                "enrichment_all": round(h / n / si.p0, 2) if n and si.p0 else np.nan,
                "p_binomial_all": p_base})

    row["n_sampled"] = len(samp)
    for k in ("A-C|hit", "A-C|miss", "D-E|hit", "D-E|miss"):
        row[f"n_sampled_{k.replace('|', '_')}"] = int((samp["stratum"] == k).sum())
    if not coverage["verifiable"]:
        return row

    sv = samp.copy()
    for r in RADII_M:
        c = f"verified_{r}m"
        sv[c] = sv[c].astype(bool)
        row[f"sample_pct_verified_{r}m"] = pct(sv[c])
        w = sv["weight"]
        row[f"sample_pct_verified_{r}m_weighted"] = (
            round(float((sv[c] * w).sum() / w.sum()) * 100, 1) if len(sv) else np.nan)
    for g in ("A-C", "D-E"):
        sub = sv[sv["grade_group"] == g]
        row[f"sample_pct_verified_250m_grade_{g}"] = pct(sub["verified_250m"])
    for hs in ("hit", "miss"):
        sub = sv[sv["hit_status"] == hs]
        row[f"sample_pct_verified_250m_{hs}"] = pct(sub["verified_250m"])

    dv = dep.copy()
    for r in RADII_M:
        dv[f"verified_{r}m"] = dv[f"verified_{r}m"].astype(bool)
        row[f"all_pct_within_{r}m"] = pct(dv[f"verified_{r}m"])
    for hs in ("hit", "miss"):
        row[f"all_pct_within_250m_{hs}"] = pct(dv.loc[dv["hit_status"] == hs, "verified_250m"])
    for g in ("A-C", "D-E"):
        row[f"all_pct_within_250m_grade_{g}"] = pct(dv.loc[dv["grade_group"] == g, "verified_250m"])
    row["all_median_dist_usmin_m"] = round(float(dv["dist_usmin_m"].median()), 0)
    _, odds, p_f = fisher_hit_vs_verified(dv, "verified_250m")
    row["all_fisher_or_hit_vs_verified_250m"] = odds
    row["all_fisher_p_hit_vs_verified_250m"] = p_f

    corr = dv[dv["verified_250m"]]
    nc, hc = len(corr), int(corr["inside_zone"].sum())
    p_c, exp_c = run_binomial(hc, nc, si.p0)
    row.update({
        "corroborated_n": nc,
        "corroborated_hits": hc,
        "corroborated_hit_rate_pct": round(hc / nc * 100, 2) if nc else np.nan,
        "corroborated_expected_hits": round(exp_c, 1),
        "corroborated_enrichment": round(hc / nc / si.p0, 2) if nc and si.p0 else np.nan,
        "corroborated_p_binomial": p_c,
        "corroborated_sig_05": p_c < 0.05,
    })
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", nargs="+", default=list(FDR_SITES))
    ap.add_argument("--sample-size", type=int, default=20)
    ap.add_argument("--refresh", action="store_true", help="Re-download cached USMIN responses")
    args = ap.parse_args()

    samples, census, summaries = [], [], []
    for site_id in args.sites:
        print(f"{site_id:22s}", end=" ", flush=True)
        si = load_site_inputs(site_id, REPO_ROOT, CACHE_DIR, fetch_osm=False, refresh=args.refresh)
        dep, coverage = site_records(si, args.refresh)
        samp = stratified_sample(dep, args.sample_size)
        summaries.append(summarise_site(si, dep, samp, coverage))
        for frame, bucket in ((samp, samples), (dep, census)):
            f = pd.DataFrame(frame.drop(columns="geometry"))
            f.insert(0, "site_id", site_id)
            bucket.append(f)
        s = summaries[-1]
        print(f"usmin={coverage['usmin_status']}  n_crit={s['n_crit']}  sampled={s['n_sampled']}  "
              f"all<=250m={s.get('all_pct_within_250m', 'NA')}%")

    samp_all = pd.concat(samples, ignore_index=True)
    cens_all = pd.concat(census, ignore_index=True)
    rec_cols = (["site_id"] + [c for c in RECORD_COLS if c in samp_all.columns]
                + ["grade_group", "inside_zone", "hit_status", "stratum", "stratum_N", "stratum_n", "weight",
                   "dist_usmin_m", "nearest_usmin_type", "nearest_usmin_layer"]
                + [f"verified_{r}m" for r in RADII_M] + ["verification_status"])
    samp_all = samp_all[rec_cols].rename(columns={"site_name": "mrds_site_name"})
    OUT_RECORDS.parent.mkdir(parents=True, exist_ok=True)
    samp_all.rename(columns={"verified_250m": "verified_usmin"}).to_csv(OUT_RECORDS, index=False)

    ver_s = samp_all[samp_all["verification_status"] != "unverifiable_no_reference"].copy()
    ver_c = cens_all[cens_all["verification_status"] != "unverifiable_no_reference"].copy()
    for d in (ver_s, ver_c):
        for r in RADII_M:
            d[f"verified_{r}m"] = d[f"verified_{r}m"].astype(bool)
        d["inside_zone"] = d["inside_zone"].astype(bool)

    pooled = {"site_id": "ALL_SITES", "site_name": "Pooled across verifiable sites",
              "n_crit": len(cens_all), "hits_crit": int(cens_all["inside_zone"].sum()),
              "n_sampled": len(samp_all)}
    for r in RADII_M:
        pooled[f"sample_pct_verified_{r}m"] = pct(ver_s[f"verified_{r}m"])
        pooled[f"all_pct_within_{r}m"] = pct(ver_c[f"verified_{r}m"])
    for g in ("A-C", "D-E"):
        pooled[f"sample_pct_verified_250m_grade_{g}"] = pct(ver_s.loc[ver_s["grade_group"] == g, "verified_250m"])
        pooled[f"all_pct_within_250m_grade_{g}"] = pct(ver_c.loc[ver_c["grade_group"] == g, "verified_250m"])
    for hs in ("hit", "miss"):
        pooled[f"sample_pct_verified_250m_{hs}"] = pct(ver_s.loc[ver_s["hit_status"] == hs, "verified_250m"])
        pooled[f"all_pct_within_250m_{hs}"] = pct(ver_c.loc[ver_c["hit_status"] == hs, "verified_250m"])

    print("\nHits vs misses: USMIN verification (Fisher exact, pooled; CMH stratified by site)")
    tests = []
    for label, d in (("sample", ver_s), ("all", ver_c)):
        for r in RADII_M:
            t, odds, p = fisher_hit_vs_verified(d, f"verified_{r}m")
            or_mh, p_cmh = cmh_by_site(d, f"verified_{r}m")
            pooled[f"{label}_fisher_or_hit_vs_verified_{r}m"] = odds
            pooled[f"{label}_fisher_p_hit_vs_verified_{r}m"] = p
            pooled[f"{label}_cmh_or_hit_vs_verified_{r}m"] = or_mh
            pooled[f"{label}_cmh_p_hit_vs_verified_{r}m"] = p_cmh
            tests.append({"set": label, "radius_m": r,
                          "hit_ver": t[0][0], "hit_not": t[0][1], "miss_ver": t[1][0], "miss_not": t[1][1],
                          "pct_ver_hit": round(100 * t[0][0] / max(sum(t[0]), 1), 1),
                          "pct_ver_miss": round(100 * t[1][0] / max(sum(t[1]), 1), 1),
                          "fisher_OR": round(odds, 2), "fisher_p": f"{p:.3g}",
                          "CMH_OR": round(or_mh, 2), "CMH_p": f"{p_cmh:.3g}"})
    print(pd.DataFrame(tests).to_string(index=False))

    summary = pd.DataFrame(summaries + [pooled])
    summary.to_csv(OUT_SUMMARY, index=False)
    print(f"\nSaved {OUT_RECORDS} ({len(samp_all)} records)")
    print(f"Saved {OUT_SUMMARY}")

    view = summary[summary["site_id"] != "ALL_SITES"]
    cols = ["site_id", "usmin_status", "n_quads_footprint", "n_quads_with_usmin", "n_crit", "n_sampled",
            "sample_pct_verified_250m", "sample_pct_verified_1000m",
            "all_pct_within_250m", "all_pct_within_250m_hit", "all_pct_within_250m_miss",
            "enrichment_all", "p_binomial_all",
            "corroborated_n", "corroborated_hits", "corroborated_enrichment", "corroborated_p_binomial"]
    fmt = view[cols].copy()
    for c in ("p_binomial_all", "corroborated_p_binomial"):
        fmt[c] = fmt[c].map(lambda v: f"{v:.2g}")
    print("\n" + fmt.to_string(index=False))

    print("\nPooled verification (sample / all records):")
    for r in RADII_M:
        print(f"  <= {r:4d} m: {pooled[f'sample_pct_verified_{r}m']}% / {pooled[f'all_pct_within_{r}m']}%")
    print(f"  by grade @250 m  A-C: {pooled['sample_pct_verified_250m_grade_A-C']}% / "
          f"{pooled['all_pct_within_250m_grade_A-C']}%   D-E: {pooled['sample_pct_verified_250m_grade_D-E']}% / "
          f"{pooled['all_pct_within_250m_grade_D-E']}%")
    print(f"  by status @250 m hit: {pooled['sample_pct_verified_250m_hit']}% / {pooled['all_pct_within_250m_hit']}%"
          f"   miss: {pooled['sample_pct_verified_250m_miss']}% / {pooled['all_pct_within_250m_miss']}%")

    lost = view[(view["p_binomial_all"] < 0.05) & (view["corroborated_p_binomial"] >= 0.05)]
    print("\nSites significant on all records but NOT on the USMIN-corroborated subset "
          "(n, hits vs expected under p0; expected < 5 = underpowered):")
    for _, r in lost.iterrows():
        print(f"  {r['site_id']}{' (robust)' if r['robust'] else ''}: n={int(r['corroborated_n'])}, "
              f"hits={int(r['corroborated_hits'])} vs {r['corroborated_expected_hits']} expected, "
              f"p={r['corroborated_p_binomial']:.2g}; USMIN quads with features "
              f"{int(r['n_quads_with_usmin'])}/{int(r['n_quads_footprint'])}")
    if lost.empty:
        print("  none")


if __name__ == "__main__":
    main()
