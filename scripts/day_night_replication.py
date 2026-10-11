#!/usr/bin/env python3
"""Day vs night AST_05 replication of the strong-anomaly result.

Emissivity is a surface property, so a real lithological / alteration signal
should survive a change of acquisition (day → night: different surface
temperature, different TES conditions). For every site with both a day
mosaic (``aster/{id}_mosaic_*``, AST_05) and a night mosaic
(``aster_night/{id}_night_mosaic_*``, built by
``rebuild_emissivity_mosaics.py --night``):

1. Night QI / CI / MI are resampled (nearest) onto the day grid; only pixels
   valid in both are used.
2. Spearman correlation per index (pixel level).
3. Both acquisitions are classified identically on the common pixels
   (site percentiles, combined score >= strong_score_min) and the strong-zone
   masks compared (Jaccard; share of day zone pixels also strong at night).
4. Bedrock/cover-stratified null on the common pixels for each acquisition:
   per-deposit p_i = strong-pixel share of the deposit's stratum, exact
   Poisson-binomial tail (same null as site_specific_null_significance.py,
   evaluated on rasters so day and night share one denominator).

Night results are a replication check only; they never enter the primary
score or the FDR table.

Outputs:
    results/day_night_replication.csv
    figures/verification_day_night_replication.png
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from critical_minerals_aster.classification import classify_percentiles, combined_score  # noqa: E402
from critical_minerals_aster.config import ClassificationParams  # noqa: E402
from critical_minerals_aster.lithology import (  # noqa: E402
    STRATUM_CODES,
    assign_cover_stratum,
    rasterize_cover_strata,
    stratified_null,
)
from site_specific_null_significance import load_site_geometry, stratified_site  # noqa: E402

warnings.filterwarnings("ignore")

INDICES = ("silica", "carbonate", "mafic")
OUT_CSV = REPO / "results" / "day_night_replication.csv"
OUT_FIG = REPO / "figures" / "verification_day_night_replication.png"


def _read(path: Path):
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float64)
        arr[arr == 0] = np.nan
        return arr, src.transform, src.crs


def _night_on_day_grid(night_path: Path, transform, crs, shape) -> np.ndarray:
    out = np.full(shape, np.nan, dtype=np.float64)
    with rasterio.open(night_path) as src:
        arr = src.read(1).astype(np.float64)
        arr[arr == 0] = np.nan
        reproject(arr, out, src_transform=src.transform, src_crs=src.crs,
                  dst_transform=transform, dst_crs=crs, resampling=Resampling.nearest,
                  src_nodata=np.nan, dst_nodata=np.nan)
    return out


def _strong(layers: dict[str, np.ndarray], common: np.ndarray, cp: ClassificationParams) -> np.ndarray:
    cls = []
    for n in INDICES:
        a = np.where(common, layers[n], np.nan)
        c, _, _ = classify_percentiles(a, cp.low_pct, cp.high_pct)
        cls.append(c)
    return (combined_score(*cls) >= cp.strong_score_min) & common


def _stratified(strong: np.ndarray, common: np.ndarray, strata_grid: np.ndarray,
                rows: np.ndarray, cols: np.ndarray, labels: pd.Series) -> dict:
    p = {}
    for name, code in STRATUM_CODES.items():
        m = common & (strata_grid == code)
        p[name] = float(strong[m].mean()) if m.any() else 0.0
    hits = int(strong[rows, cols].sum())
    p_i = labels.map(p).to_numpy(dtype=float)
    return stratified_null(p_i, hits) | {"zone_frac": float(strong[common].mean())}


def site_row(site_id: str) -> dict | None:
    day_dir = REPO / "data" / "sites" / site_id / "aster"
    night_dir = REPO / "data" / "sites" / site_id / "aster_night"
    sel_day = day_dir / f"{site_id}_mosaic_granules.json"
    sel_night = night_dir / f"{site_id}_night_mosaic_granules.json"
    if not (sel_day.is_file() and sel_night.is_file()):
        return None
    if json.loads(sel_day.read_text()).get("product") != "AST_05":
        return None
    night_meta = json.loads(sel_night.read_text())

    site, zones, _, fp_geom, raster_bbox, deposits = load_site_geometry(site_id)
    cp = site.classification or ClassificationParams()
    day, night = {}, {}
    for n in INDICES:
        day[n], transform, crs = _read(day_dir / f"{site_id}_mosaic_ratio_{n}.tif")
        night[n] = _night_on_day_grid(night_dir / f"{site_id}_night_mosaic_ratio_{n}.tif",
                                      transform, crs, day[n].shape)
    common = np.ones(day["silica"].shape, dtype=bool)
    for n in INDICES:
        common &= np.isfinite(day[n]) & np.isfinite(night[n])
    day_valid = np.isfinite(day["silica"])
    row = {
        "site_id": site_id,
        "night_tier": night_meta.get("tier"),
        "night_months": ";".join(str(g.get("month")) for g in night_meta.get("granules", [])),
        "common_frac_of_day_valid": round(float(common.sum() / max(day_valid.sum(), 1)), 4),
    }
    if common.sum() < 1000:
        return row | {"note": "too few common pixels"}

    rng = np.random.default_rng(0)
    idx = np.flatnonzero(common)
    take = rng.choice(idx, size=min(200_000, idx.size), replace=False)
    for n in INDICES:
        row[f"spearman_{n}"] = round(float(spearmanr(day[n].ravel()[take], night[n].ravel()[take])[0]), 4)

    s_day, s_night = _strong(day, common, cp), _strong(night, common, cp)
    union = (s_day | s_night).sum()
    row["jaccard_strong"] = round(float((s_day & s_night).sum() / union), 4) if union else np.nan
    row["day_zone_also_night"] = round(float((s_day & s_night).sum() / max(s_day.sum(), 1)), 4)
    # Jaccard expected if day and night zones were independent at these sizes.
    a, b = s_day[common].mean(), s_night[common].mean()
    row["jaccard_independent"] = round(float(a * b / (a + b - a * b)), 4) if (a + b) else np.nan

    strata, _, _, _ = stratified_site(site_id, zones, fp_geom, raster_bbox, deposits)
    strata_grid = rasterize_cover_strata(strata, transform, common.shape)
    if len(deposits):
        dep = deposits.to_crs(crs)
        rr, cc = rasterio.transform.rowcol(transform, dep.geometry.x.values, dep.geometry.y.values)
        rr, cc = np.asarray(rr), np.asarray(cc)
        inside = (rr >= 0) & (rr < common.shape[0]) & (cc >= 0) & (cc < common.shape[1])
        keep = np.zeros(len(dep), dtype=bool)
        keep[inside] = common[rr[inside], cc[inside]]
        dep, rr, cc = dep[keep], rr[keep], cc[keep]
        labels = assign_cover_stratum(dep, {k: v for k, v in strata.items()})
    else:
        rr = cc = np.array([], dtype=int)
        labels = pd.Series(dtype=object)
    row["n_crit_common"] = int(len(labels))
    for tag, strong in (("day", s_day), ("night", s_night)):
        res = _stratified(strong, common, strata_grid, rr, cc, labels)
        row[f"{tag}_zone_frac"] = round(res["zone_frac"], 4)
        row[f"{tag}_hits"] = res["hits"]
        row[f"{tag}_expected"] = round(res["expected"], 3)
        row[f"{tag}_enrichment"] = round(res["enrichment"], 3) if np.isfinite(res["enrichment"]) else np.nan
        row[f"{tag}_p"] = res["p"]
    row["both_hits"] = int((s_day[rr, cc] & s_night[rr, cc]).sum()) if len(rr) else 0
    return row


def figure(df: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    ok = df.dropna(subset=["day_enrichment", "night_enrichment"])
    ok = ok[ok["n_crit_common"] >= 5]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    sig = (ok["day_p"] < 0.05) | (ok["night_p"] < 0.05)
    ax.scatter(ok.loc[~sig, "day_enrichment"], ok.loc[~sig, "night_enrichment"], s=18, c="0.6",
               label="neither p < 0.05")
    ax.scatter(ok.loc[sig, "day_enrichment"], ok.loc[sig, "night_enrichment"], s=28, c="C3",
               label="day or night p < 0.05")
    both = (ok["day_p"] < 0.05) & (ok["night_p"] < 0.05)
    ax.scatter(ok.loc[both, "day_enrichment"], ok.loc[both, "night_enrichment"], s=60,
               facecolors="none", edgecolors="k", label="both p < 0.05 (listed)")
    listing = "\n".join(
        f"{r.site_id}: {r.day_enrichment:.2f} / {r.night_enrichment:.2f}"
        for r in ok[both].sort_values("day_p").itertuples()
    )
    if listing:
        ax.text(0.98, 0.02, "Both p < 0.05 (day / night):\n" + listing, transform=ax.transAxes,
                fontsize=7, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, edgecolor="0.7"))
    lim = max(3.0, float(np.nanmax(ok[["day_enrichment", "night_enrichment"]].values)) * 1.05)
    ax.plot([0, lim], [0, lim], "k--", lw=0.8)
    ax.axhline(1, c="0.8", lw=0.8)
    ax.axvline(1, c="0.8", lw=0.8)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Day enrichment (hits / stratified expected)")
    ax.set_ylabel("Night enrichment")
    ax.set_title(f"Deposit enrichment, day vs night\n(AST_05, {len(ok)} sites with ≥5 deposits)")
    ax.legend(fontsize=8, loc="upper left")

    ax = axes[1]
    data = [df[f"spearman_{n}"].dropna() for n in INDICES]
    ax.boxplot(data, labels=["QI (quartz)", "CI (carbonate)", "MI (mafic)"])
    ax.axhline(0, c="0.8", lw=0.8)
    ax.set_ylabel("Pixel Spearman ρ, day vs night")
    ax.set_title("Index agreement across acquisitions")

    ax = axes[2]
    d = df.dropna(subset=["jaccard_strong", "jaccard_independent"])
    ax.scatter(d["jaccard_independent"], d["jaccard_strong"], s=18)
    m = max(0.3, float(d[["jaccard_strong", "jaccard_independent"]].max().max()) * 1.05)
    ax.plot([0, m], [0, m], "k--", lw=0.8, label="no better than independent")
    ax.set_xlim(0, m)
    ax.set_ylim(0, m)
    ax.set_xlabel("Jaccard expected if independent")
    ax.set_ylabel("Observed Jaccard (strong zones)")
    ax.set_title("Strong-zone overlap, day vs night")
    ax.legend(fontsize=8)

    fig.tight_layout()
    OUT_FIG.parent.mkdir(exist_ok=True)
    fig.savefig(OUT_FIG, dpi=150)
    plt.close(fig)


def main() -> int:
    import yaml

    if "--figure-only" in sys.argv[1:]:
        figure(pd.read_csv(OUT_CSV))
        print(f"Wrote {OUT_FIG.relative_to(REPO)}")
        return 0
    site_ids = yaml.safe_load((REPO / "sites" / "index.yaml").read_text())["sites"]
    rows = []
    for sid in site_ids:
        try:
            r = site_row(sid)
        except Exception as exc:
            print(f"{sid}: FAILED ({exc!r})", file=sys.stderr)
            continue
        if r is not None:
            rows.append(r)
            print(f"{sid}: jaccard={r.get('jaccard_strong')} "
                  f"day_enr={r.get('day_enrichment')} night_enr={r.get('night_enrichment')}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV.relative_to(REPO)} ({len(df)} sites)")
    if len(df):
        figure(df)
        print(f"Wrote {OUT_FIG.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
