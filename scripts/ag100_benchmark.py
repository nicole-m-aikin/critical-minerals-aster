#!/usr/bin/env python
"""
AG100 benchmark: the strong-anomaly / critical-deposit test on the ASTER GED
2000–2008 mean emissivity (AG100 v003) instead of a single per-scene granule.

Per site (sites/index.yaml):

1. AG100 mean emissivity B10–B14 mosaicked over the site bbox, resampled
   (bilinear) to a 100 m UTM grid; water / NDVI / observation count by
   nearest / bilinear / nearest (critical_minerals_aster.ag100).
2. Valid pixels: inside the bbox, all five emissivities present, not water
   (AG100 LWmap), NDVI ≤ NDVI_VEGETATION_THRESHOLD (0.3), NumObs ≥ 3.
3. Ninomiya quartz / carbonate / mafic indices (spectral.alteration_indices)
   on valid pixels, per-site 70th/90th percentile classes, combined score,
   strong zones = score ≥ 3 (classification.vectorize_strong_zones).
4. Deposits: validated national MRDS → site bbox → critical minerals
   (reclassify_mrds_earth_mri + is_critical_mineral, as in
   scripts/lithology_matched_null.py) → deposits on valid pixels. Hits =
   deposits within a zone polygon (mrds.spatial_join_deposits_zones).
5. Geometric null: p0 = zone pixels / valid pixels, one-sided exact binomial.
6. Bedrock / cover null: SGMC (data/lithology/{site}_sgmc_units.geojson,
   fetched if missing) coarse class ``unconsolidated`` → cover, all other
   mapped classes → bedrock, uncovered → unmapped. Strata are rasterized on
   the analysis grid (overlaps resolve bedrock-first, as in
   lithology.build_strata's alphabetical order); per-stratum p0 from pixel
   counts; exact Poisson-binomial p via lithology.stratified_test.
7. BH-FDR (α = 0.05) across all analysed sites for both nulls.

The comparison panel reads per-scene enrichment / FDR flags from
results/phase4_fdr_corrected_significance.csv as it exists at run time.

Outputs
-------
    results/ag100_benchmark.csv
    figures/verification_ag100_benchmark.png

Usage:
    conda run -n aster-minerals python scripts/ag100_benchmark.py
    conda run -n aster-minerals python scripts/ag100_benchmark.py --sites bisbee holden
    conda run -n aster-minerals python scripts/ag100_benchmark.py --figure-only
"""

from __future__ import annotations

import argparse
import math
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rasterio.crs import CRS
from rasterio.features import geometry_mask, rasterize
from rasterio.transform import from_origin
from rasterio.warp import transform_bounds
from shapely.geometry import box
from statsmodels.stats.multitest import multipletests

from critical_minerals_aster.ag100 import (
    MIN_OBS_COUNT,
    NDVI_VEGETATION_THRESHOLD,
    masks_from_layers,
    mosaic_to_grid,
    read_ag100_mosaic,
)
from critical_minerals_aster.classification import (
    classify_percentiles,
    combined_score,
    vectorize_strong_zones,
)
from critical_minerals_aster.config import list_site_ids, load_site_by_id
from critical_minerals_aster.lithology import add_lith_class, load_sgmc_units, stratified_test
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    is_critical_mineral,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.significance import run_binomial
from critical_minerals_aster.spectral import alteration_indices

REPO_ROOT = Path(__file__).parent.parent
RES_M = 100.0
LOW_PCT, HIGH_PCT = 70.0, 90.0
MIN_SCORE = 3
ALPHA = 0.05
COVER, BEDROCK, UNMAPPED = "cover", "bedrock", "unmapped"
PER_SCENE_CSV = REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv"


def utm_crs(bbox: tuple[float, float, float, float]) -> CRS:
    lon_c, lat_c = (bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2
    zone = int((lon_c + 180) // 6) + 1
    return CRS.from_epsg((32600 if lat_c >= 0 else 32700) + zone)


def site_grid(bbox):
    """100 m UTM grid enclosing the bbox, plus the in-bbox pixel mask."""
    crs = utm_crs(bbox)
    left, bottom, right, top = transform_bounds("EPSG:4326", crs, *bbox, densify_pts=51)
    x0, y1 = math.floor(left / RES_M) * RES_M, math.ceil(top / RES_M) * RES_M
    shape = (math.ceil((y1 - bottom) / RES_M), math.ceil((right - x0) / RES_M))
    transform = from_origin(x0, y1, RES_M, RES_M)
    bbox_poly = gpd.GeoSeries([box(*bbox).segmentize(0.01)], crs="EPSG:4326").to_crs(crs).iloc[0]
    in_bbox = geometry_mask([bbox_poly], out_shape=shape, transform=transform, invert=True)
    return crs, transform, shape, in_bbox


def pixel_index(transform, shape, geoms) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cols, rows = ~transform * (geoms.x.values, geoms.y.values)
    rows, cols = np.floor(rows).astype(int), np.floor(cols).astype(int)
    inside = (rows >= 0) & (rows < shape[0]) & (cols >= 0) & (cols < shape[1])
    return np.where(inside, rows, 0), np.where(inside, cols, 0), inside


def cover_raster(site_id, bbox, crs, transform, shape) -> tuple[np.ndarray, int]:
    """0 unmapped, 1 bedrock, 2 cover on the analysis grid; also the SGMC unit count."""
    units = load_sgmc_units(site_id, bbox, REPO_ROOT)
    out = np.zeros(shape, dtype=np.uint8)
    if not len(units):
        return out, 0
    units = add_lith_class(units.to_crs(crs))
    units = units[units.geometry.notna() & ~units.geometry.is_empty]
    is_cover = units["lith_class"] == "unconsolidated"
    for code, sub in ((2, units[is_cover]), (1, units[~is_cover])):
        if len(sub):
            burned = rasterize(((g, 1) for g in sub.geometry), out_shape=shape, transform=transform,
                               fill=0, dtype=np.uint8)
            out[burned == 1] = code
    return out, len(units)


def analyse_site(site_id: str, mrds: pd.DataFrame) -> dict:
    site = load_site_by_id(site_id, REPO_ROOT / "sites")
    bbox = tuple(site.bbox_wgs84)
    row: dict = {"site_id": site_id, "site_name": site.name}
    crs, transform, shape, in_bbox = site_grid(bbox)
    mosaic = read_ag100_mosaic(bbox, REPO_ROOT, ("emissivity", "water", "ndvi", "obs_count"))
    if mosaic is None:
        row["status"] = "no_ag100"
        return row
    grids = mosaic_to_grid(mosaic, transform, shape, crs)
    masks = masks_from_layers(grids)
    emis = grids["emissivity"]
    ndvi, obs = masks["ndvi"], masks["obs_count"]

    emis_ok = np.isfinite(emis).all(axis=0)
    water = masks["water"] & in_bbox
    veg = (ndvi > NDVI_VEGETATION_THRESHOLD) & in_bbox
    low_obs = ~(obs >= MIN_OBS_COUNT) & in_bbox
    valid = in_bbox & emis_ok & ~masks["water"] & np.isfinite(ndvi) & ~veg & (obs >= MIN_OBS_COUNT)

    n_bbox = int(in_bbox.sum())
    n_valid = int(valid.sum())
    px_km2 = RES_M ** 2 / 1e6
    row.update({
        "status": "ok",
        "utm_crs": crs.to_string(),
        "n_tiles": len(mosaic.tiles),
        "bbox_km2": round(n_bbox * px_km2, 1),
        "valid_km2": round(n_valid * px_km2, 1),
        "valid_frac": round(n_valid / n_bbox, 4),
        "water_frac": round(int(water.sum()) / n_bbox, 4),
        "veg_frac": round(int(veg.sum()) / n_bbox, 4),
        "low_obs_frac": round(int(low_obs.sum()) / n_bbox, 4),
        "emis_missing_frac": round(int((~emis_ok & in_bbox).sum()) / n_bbox, 4),
        "ndvi_median": round(float(np.nanmedian(ndvi[in_bbox])), 3),
        "obs_median": float(np.nanmedian(obs[in_bbox])),
    })
    if n_valid < 100:
        row["status"] = "too_few_valid_pixels"
        return row

    indices = alteration_indices(*emis)
    classes = []
    for name, idx in zip(("quartz", "carbonate", "mafic"), indices):
        idx = np.where(valid, idx, np.nan)
        cls, lo, hi = classify_percentiles(idx, LOW_PCT, HIGH_PCT)
        row[f"{name}_p70"], row[f"{name}_p90"] = round(lo, 5), round(hi, 5)
        classes.append(cls)
    score = combined_score(*classes)
    zone_px = (score >= MIN_SCORE) & valid
    zones = vectorize_strong_zones(score, transform, crs, min_score=MIN_SCORE)
    p0 = int(zone_px.sum()) / n_valid
    row.update({"n_zones": len(zones), "zone_km2": round(int(zone_px.sum()) * px_km2, 2), "p0": round(p0, 4)})

    deposits = mrds_to_points_gdf(filter_mrds_bbox(mrds, bbox), crs)
    deposits = reclassify_mrds_earth_mri(deposits)
    deposits = deposits[deposits["earth_mri_category"].apply(is_critical_mineral)].copy()
    r, c, inside = pixel_index(transform, shape, deposits.geometry)
    keep = inside & valid[r, c]
    n_masked = int((inside & in_bbox[r, c] & ~valid[r, c]).sum())
    deposits, r, c = deposits[keep].copy(), r[keep], c[keep]
    n = len(deposits)
    if n and len(zones):
        _, hits_gdf, _ = spatial_join_deposits_zones(deposits, zones)
        deposits["hit"] = deposits.index.isin(hits_gdf.index.unique())
    else:
        deposits["hit"] = False
    hits = int(deposits["hit"].sum())
    p_geo, expected = run_binomial(hits, n, p0)
    row.update({
        "n_crit_bbox": n + n_masked,
        "n_crit_masked": n_masked,
        "n_crit": n,
        "hits": hits,
        "hits_raster_check": int(zone_px[r, c].sum()) if n else 0,
        "hit_rate_pct": round(100 * hits / n, 2) if n else np.nan,
        "expected": round(expected, 2),
        "enrichment": round(hits / expected, 3) if expected > 0 else np.nan,
        "p_geo": p_geo,
    })

    cover, n_units = cover_raster(site_id, bbox, crs, transform, shape)
    dep_stratum = cover[r, c]
    tab_rows = []
    for code, label in ((0, UNMAPPED), (1, BEDROCK), (2, COVER)):
        s = valid & (cover == code)
        tab_rows.append({
            "stratum": label,
            "area": float(s.sum()),
            "zone_area": float((s & zone_px).sum()),
            "n": int((dep_stratum == code).sum()),
            "hits": int(deposits["hit"].values[dep_stratum == code].sum()),
        })
    tab = pd.DataFrame(tab_rows)
    tab["p0"] = np.where(tab["area"] > 0, tab["zone_area"] / tab["area"].where(tab["area"] > 0, 1.0), 0.0)
    strat = stratified_test(tab)
    t = tab.set_index("stratum")
    row.update({
        "n_sgmc_units": n_units,
        "cover_frac": round(t.loc[COVER, "area"] / n_valid, 4),
        "unmapped_frac": round(t.loc[UNMAPPED, "area"] / n_valid, 4),
        "p0_bedrock": round(t.loc[BEDROCK, "p0"], 4),
        "p0_cover": round(t.loc[COVER, "p0"], 4),
        "n_bedrock": int(t.loc[BEDROCK, "n"]),
        "n_cover": int(t.loc[COVER, "n"]),
        "hits_bedrock": int(t.loc[BEDROCK, "hits"]),
        "hits_cover": int(t.loc[COVER, "hits"]),
        "expected_strat": round(strat["expected"], 2),
        "enrichment_strat": round(strat["enrichment"], 3) if n else np.nan,
        "p_strat": strat["p"],
    })
    return row


def add_fdr(df: pd.DataFrame, p_col: str, q_col: str, sig_col: str) -> None:
    ok = df[p_col].notna()
    df[q_col] = np.nan
    df[sig_col] = False
    if ok.any():
        reject, q, _, _ = multipletests(df.loc[ok, p_col].values, alpha=ALPHA, method="fdr_bh")
        df.loc[ok, q_col] = q
        df.loc[ok, sig_col] = reject


def per_scene_reference() -> tuple[pd.DataFrame | None, str]:
    if not PER_SCENE_CSV.is_file():
        return None, "missing"
    ref = pd.read_csv(PER_SCENE_CSV)
    stamp = datetime.fromtimestamp(PER_SCENE_CSV.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
    cols = [c for c in ("site_id", "n_crit", "hits_crit", "p0_null", "enrichment", "sig_fdr_binomial")
            if c in ref.columns]
    return ref[cols], stamp


def make_figure(df: pd.DataFrame, out_path: Path) -> None:
    ok = df[df["status"] == "ok"].drop(columns=["enrichment_scene", "sig_fdr_scene", "n_crit_scene"],
                                       errors="ignore")
    ref, stamp = per_scene_reference()
    if ref is not None and "enrichment" in ref.columns:
        ok = ok.merge(ref.rename(columns={"enrichment": "enrichment_scene", "sig_fdr_binomial": "sig_scene",
                                          "n_crit": "n_crit_scene"}), on="site_id", how="left")
    else:
        ok["enrichment_scene"], ok["sig_scene"] = np.nan, False
    ok["sig_scene"] = ok["sig_scene"].fillna(False).astype(bool)

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1], height_ratios=[1.35, 1], hspace=0.32, wspace=0.28)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_b = fig.add_subplot(gs[:, 1])

    # A — AG100 vs per-scene enrichment, both under the bedrock/cover stratified null
    # (the per-scene phase4 ``enrichment`` / ``sig_fdr_binomial`` columns are stratified).
    sc = ok[ok["enrichment_strat"].notna() & ok["enrichment_scene"].notna()].copy()
    sc["enrichment"], sc["sig_geo"] = sc["enrichment_strat"], sc["sig_strat"].astype(bool)
    groups = [
        ("Both FDR-sig", sc["sig_geo"] & sc["sig_scene"], "#c0392b", "o"),
        ("AG100 only", sc["sig_geo"] & ~sc["sig_scene"], "#e67e22", "s"),
        ("Per-scene only", ~sc["sig_geo"] & sc["sig_scene"], "#2c7fb8", "D"),
        ("Neither", ~sc["sig_geo"] & ~sc["sig_scene"], "0.65", "o"),
    ]
    lim = float(np.nanmax([sc["enrichment"].max(), sc["enrichment_scene"].max(), 2.0])) * 1.08 if len(sc) else 3
    ax_a.plot([0, lim], [0, lim], ls="--", color="0.5", lw=1, zorder=0, label="1:1")
    ax_a.axhline(1, color="0.85", lw=0.8, zorder=0)
    ax_a.axvline(1, color="0.85", lw=0.8, zorder=0)
    for label, m, colour, marker in groups:
        sub = sc[m]
        if len(sub):
            ax_a.scatter(sub["enrichment_scene"], sub["enrichment"], s=18 + 2.2 * np.sqrt(sub["n_crit"]),
                         c=colour, marker=marker, edgecolors="k", linewidths=0.4, alpha=0.85,
                         label=f"{label} (n={len(sub)})", zorder=3)
    for _, r in sc[sc["sig_geo"] | sc["sig_scene"]].iterrows():
        ax_a.annotate(r["site_id"], (r["enrichment_scene"], r["enrichment"]), xytext=(3, 3),
                      textcoords="offset points", fontsize=6.5, color="0.25")
    rho = sc[["enrichment", "enrichment_scene"]].corr(method="spearman").iloc[0, 1] if len(sc) > 2 else np.nan
    ax_a.set_xlim(0, lim)
    ax_a.set_ylim(0, lim)
    ax_a.set_xlabel(f"Per-scene AST_05 enrichment, stratified (phase4 file, {stamp})")
    ax_a.set_ylabel("AG100 mean-emissivity enrichment (stratified null)")
    ax_a.set_title(f"A. AG100 vs per-scene enrichment, critical-mineral deposits\n"
                   f"Spearman ρ = {rho:.2f} ({len(sc)} sites with both)", loc="left", fontsize=11)
    ax_a.legend(fontsize=8, loc="upper right", title="BH-FDR q < 0.05 (stratified)", title_fontsize=8,
                framealpha=0.9)
    ax_a.grid(alpha=0.25)

    # B — AG100 enrichment with FDR flags, for every site flagged anywhere
    flag = ok[ok["sig_geo"] | ok["sig_strat"] | ok["sig_scene"]].sort_values("enrichment")
    y = np.arange(len(flag))
    ax_b.axvline(1, color="0.5", lw=1, ls="--")
    for col, sig, colour, marker, dy, label in (
        ("enrichment", "sig_geo", "#c0392b", "o", 0.17, "Geometric null"),
        ("enrichment_strat", "sig_strat", "#6a3d9a", "D", -0.17, "Bedrock/cover stratified null"),
    ):
        filled = flag[sig].values
        ax_b.scatter(flag[col].values[filled], y[filled] + dy, marker=marker, s=40, c=colour,
                     edgecolors=colour, label=f"{label}: FDR-sig")
        ax_b.scatter(flag[col].values[~filled], y[~filled] + dy, marker=marker, s=40, facecolors="white",
                     edgecolors=colour, label=f"{label}: not sig")
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([f"{s}  (n={int(n)})" for s, n in zip(flag["site_id"], flag["n_crit"])], fontsize=7.5)
    for tick, s in zip(ax_b.get_yticklabels(), flag["sig_scene"]):
        if s:
            tick.set_fontweight("bold")
    ax_b.set_ylim(-0.7, len(flag) - 0.3)
    ax_b.set_xlabel("AG100 enrichment (observed / expected hits)")
    ax_b.set_title("B. AG100 enrichment, sites FDR-significant on any test\n(bold label = per-scene FDR-sig)",
                   loc="left", fontsize=11)
    ax_b.legend(fontsize=8, loc="lower right")
    ax_b.grid(axis="x", alpha=0.25)

    # C — mask fractions
    cols = [("water_frac", "Water (LWmap)"), ("veg_frac", f"NDVI > {NDVI_VEGETATION_THRESHOLD}"),
            ("low_obs_frac", f"NumObs < {MIN_OBS_COUNT}"), ("valid_frac", "Valid (analysed)")]
    data = [100 * ok[c].dropna().values for c, _ in cols]
    bp = ax_c.boxplot(data, vert=False, widths=0.55, patch_artist=True, showfliers=True,
                      flierprops={"markersize": 3})
    for patch, colour in zip(bp["boxes"], ("#74a9cf", "#78c679", "#bdbdbd", "#fdae6b")):
        patch.set_facecolor(colour)
    ax_c.set_yticks(range(1, len(cols) + 1))
    ax_c.set_yticklabels([f"{lab}\nmedian {np.median(d):.1f}%" for (_, lab), d in zip(cols, data)], fontsize=8.5)
    ax_c.set_xlabel(f"% of site bbox ({len(ok)} sites)")
    ax_c.set_title("C. AG100 mask fractions per site", loc="left", fontsize=11)
    ax_c.grid(axis="x", alpha=0.25)

    n_ok = len(ok)
    fig.suptitle(
        f"AG100 v003 mean-emissivity benchmark — FDR-significant: geometric {int(ok['sig_geo'].sum())}/{n_ok}, "
        f"bedrock/cover stratified {int(ok['sig_strat'].sum())}/{n_ok}, per-scene {int(ok['sig_scene'].sum())}",
        fontsize=13, fontweight="bold",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--sites", nargs="*", help="subset of site ids (default: all in sites/index.yaml)")
    ap.add_argument("--figure-only", action="store_true", help="redraw the figure from results/ag100_benchmark.csv")
    args = ap.parse_args()

    out_csv = REPO_ROOT / "results" / "ag100_benchmark.csv"
    fig_path = REPO_ROOT / "figures" / "verification_ag100_benchmark.png"
    if args.figure_only:
        make_figure(pd.read_csv(out_csv), fig_path)
        print(f"Saved {fig_path}")
        return

    site_ids = args.sites or list_site_ids(REPO_ROOT / "sites")
    first = load_site_by_id(site_ids[0], REPO_ROOT / "sites")
    mrds = read_mrds_national(site_paths_for(first, REPO_ROOT))

    rows = []
    for site_id in site_ids:
        t0 = time.time()
        print(f"{site_id:26s}", end=" ", flush=True)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FutureWarning)
                row = analyse_site(site_id, mrds)
        except Exception as exc:
            print(f"FAILED — {exc!r}")
            rows.append({"site_id": site_id, "status": f"error: {exc!r}"})
            continue
        rows.append(row)
        if row["status"] != "ok":
            print(row["status"])
            continue
        print(
            f"valid={row['valid_frac']:.2f} water={row['water_frac']:.3f} veg={row['veg_frac']:.2f} "
            f"n={row['n_crit']:4d} hits={row['hits']:3d} p0={row['p0']:.3f} E={row['enrichment']} "
            f"p={row['p_geo']:.3g} E_strat={row['enrichment_strat']} p_strat={row['p_strat']:.3g} "
            f"({time.time() - t0:.1f}s)"
        )

    df = pd.DataFrame(rows)
    ok = df["status"] == "ok"
    for p_col, q_col, sig_col in (("p_geo", "q_geo", "sig_geo"), ("p_strat", "q_strat", "sig_strat")):
        sub = df[ok].copy()
        add_fdr(sub, p_col, q_col, sig_col)
        df[q_col] = sub[q_col]
        df[sig_col] = sub[sig_col]
        df[sig_col] = df[sig_col].fillna(False).astype(bool)

    ref, stamp = per_scene_reference()
    if ref is not None:
        r = ref.set_index("site_id")
        for src, dst in (("enrichment", "enrichment_scene"), ("sig_fdr_binomial", "sig_fdr_scene"),
                         ("n_crit", "n_crit_scene")):
            if src in r.columns:
                df[dst] = df["site_id"].map(r[src])
    df["ndvi_threshold"] = NDVI_VEGETATION_THRESHOLD
    df["min_obs"] = MIN_OBS_COUNT
    lead = ["site_id", "site_name", "status", "n_crit", "hits", "p0", "expected", "enrichment", "p_geo",
            "q_geo", "sig_geo", "expected_strat", "enrichment_strat", "p_strat", "q_strat", "sig_strat"]
    lead = [c for c in lead if c in df.columns]
    df = df[lead + [c for c in df.columns if c not in lead]]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\nSaved {out_csv}")

    good = df[ok]
    print(f"\nSites analysed: {int(ok.sum())}/{len(df)}  (per-scene reference file: {stamp})")
    for c in ("water_frac", "veg_frac", "low_obs_frac", "valid_frac"):
        print(f"  median {c}: {100 * good[c].median():.2f}%   (max {100 * good[c].max():.1f}%)")
    print(f"FDR-significant geometric: {int(good['sig_geo'].sum())} — {sorted(good.loc[good['sig_geo'], 'site_id'])}")
    print(f"FDR-significant stratified: {int(good['sig_strat'].sum())} — "
          f"{sorted(good.loc[good['sig_strat'], 'site_id'])}")
    mism = good[good["hits"] != good["hits_raster_check"]]
    print(f"Polygon vs raster hit mismatches: {len(mism)}")
    make_figure(df, fig_path)
    print(f"Saved {fig_path}")


if __name__ == "__main__":
    main()
