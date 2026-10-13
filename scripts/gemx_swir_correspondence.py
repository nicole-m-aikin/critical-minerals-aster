"""Overlap test: ASTER strong zones vs AVIRIS-5 doublets, any site.

Same pre-specified rules as the Goldfield test, on the shared 90 m grid
where both TIR and SWIR are valid.  Not Tetracorder.

Verdict heuristic (Al-OH primary) — see ``gemx_site_select.correspondence_verdict``:
overlapping if Jaccard ≥ 0.15 or lift ≥ 1.3; complementary if Jaccard < 0.15
and lift < 1.1; else insufficient (mixed / inconclusive); plus no_coverage / error.

Usage:
    conda run -n aster-minerals python scripts/gemx_swir_correspondence.py --site marysvale
    conda run -n aster-minerals python scripts/gemx_swir_correspondence.py --primary-fdr
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from fig_goldfield_coverage_compare import (  # noqa: E402
    _mask_to_geom,
    _valid_mask,
    _valid_on_grid,
)
from gemx_goldfield_correspondence import (  # noqa: E402
    STRONG_MIN,
    PX_M,
    _index_block,
    _mwu,
    _warp_array,
    _warp_class,
    median_enrichment,
)
from gemx_site_select import (  # noqa: E402
    add_site_select_args,
    aviris5_coverage_by_site,
    correspondence_verdict,
    resolve_pins,
    sites_from_namespace,
)
from gemx_swir_pilot import (  # noqa: E402
    _epsg_label,
    _hillshade_overlay,
    _load_hillshade,
    _load_site_basemap,
    _site_id_from_argv,
    _town_xy,
    swir_tifs_exist,
)


def _site_grid(site, crs_epsg: int):
    from rasterio.crs import CRS as RasterioCRS
    from rasterio.transform import from_bounds
    from rasterio.warp import transform_bounds

    crs = f"EPSG:{crs_epsg}"
    west, south, east, north = transform_bounds(
        "EPSG:4326", crs, *site.bbox_wgs84, densify_pts=21,
    )
    cols = max(1, int(round((east - west) / PX_M)))
    rows = max(1, int(round((north - south) / PX_M)))
    dst_t = from_bounds(west, south, east, north, cols, rows)
    return dst_t, (rows, cols), RasterioCRS.from_epsg(crs_epsg)


def _crs_epsg_from_aster(aster_dir: Path, site_id: str) -> int:
    import rasterio

    path = aster_dir / f"{site_id}_mosaic_TIR_B10.tif"
    with rasterio.open(path) as src:
        epsg = src.crs.to_epsg() if src.crs else None
    if not epsg:
        raise RuntimeError(f"no EPSG on {path}")
    return int(epsg)


def _aster_dir(site) -> Path:
    from critical_minerals_aster.paths import site_paths_for

    return site_paths_for(site, ROOT).aster_dir


def ensure_mosaic_cls(site) -> Path:
    """Write ``*_mosaic_cls_*.tif`` from existing ratio mosaics if missing.

    Bodie (and any later site) may have ratio GeoTIFFs without the
    classification rasters correspondence needs.  Percentiles are the
    site YAML 70/90 cuts, same as the pipeline.  Does not download ASTER.
    """
    from critical_minerals_aster.classification import classify_percentiles
    from critical_minerals_aster.spectral import load_ratio_mosaic

    aster_dir = _aster_dir(site)
    mosaic_id = f"{site.id}_mosaic"
    needed = [
        aster_dir / f"{mosaic_id}_cls_silica.tif",
        aster_dir / f"{mosaic_id}_cls_carbonate.tif",
        aster_dir / f"{mosaic_id}_cls_mafic.tif",
        aster_dir / f"{mosaic_id}_cls_valid.tif",
    ]
    if all(p.is_file() for p in needed):
        return aster_dir
    silica, carbonate, mafic, b10, transform, crs = load_ratio_mosaic(
        aster_dir, mosaic_id,
    )
    cp = site.classification
    if cp is None:
        raise RuntimeError(f"{site.id} has no classification params")
    silica_cls, _, _ = classify_percentiles(silica, cp.low_pct, cp.high_pct)
    carbonate_cls, _, _ = classify_percentiles(carbonate, cp.low_pct, cp.high_pct)
    mafic_cls, _, _ = classify_percentiles(mafic, cp.low_pct, cp.high_pct)
    valid = np.isfinite(b10).astype(np.uint8)
    import rasterio

    profile = {
        "driver": "GTiff",
        "dtype": "uint8",
        "count": 1,
        "height": silica_cls.shape[0],
        "width": silica_cls.shape[1],
        "transform": transform,
        "crs": crs,
        "nodata": 0,
        "compress": "lzw",
    }
    for name, arr in (
        ("silica", silica_cls),
        ("carbonate", carbonate_cls),
        ("mafic", mafic_cls),
    ):
        with rasterio.open(aster_dir / f"{mosaic_id}_cls_{name}.tif", "w", **profile) as dst:
            dst.write(arr, 1)
    with rasterio.open(aster_dir / f"{mosaic_id}_cls_valid.tif", "w", **profile) as dst:
        dst.write(valid, 1)
    return aster_dir


def compute(site) -> tuple[dict, dict]:
    from rasterio.warp import Resampling
    from scipy.ndimage import distance_transform_edt

    pins = resolve_pins(site, sites_dir=ROOT / "sites")
    aster_dir = ensure_mosaic_cls(site)
    swir_dir = ROOT / "data" / "gemx" / "aviris" / site.id
    crs_epsg = _crs_epsg_from_aster(aster_dir, site.id)
    dst_t, shape, dst_crs = _site_grid(site, crs_epsg)
    crs = f"EPSG:{crs_epsg}"

    def tir_valid(a, nodata):
        return _valid_mask(a, nodata) & (a > 0)

    def swir_valid(a, nodata):
        return _valid_mask(a, nodata)

    tir_ok = _valid_on_grid(
        aster_dir / f"{site.id}_mosaic_TIR_B10.tif",
        tir_valid, dst_t, shape, dst_crs,
    )
    swir_ok = _valid_on_grid(
        swir_dir / f"{site.id}_aloh.tif",
        swir_valid, dst_t, shape, dst_crs,
    )
    cls_ok = _valid_on_grid(
        aster_dir / f"{site.id}_mosaic_cls_valid.tif",
        lambda a, nd: a > 0, dst_t, shape, dst_crs,
    )
    overlap = tir_ok & swir_ok & cls_ok

    silica = _warp_class(aster_dir / f"{site.id}_mosaic_cls_silica.tif", dst_t, shape, dst_crs)
    carb = _warp_class(aster_dir / f"{site.id}_mosaic_cls_carbonate.tif", dst_t, shape, dst_crs)
    mafic = _warp_class(aster_dir / f"{site.id}_mosaic_cls_mafic.tif", dst_t, shape, dst_crs)
    combined = np.where(cls_ok, silica.astype(np.int16) + carb + mafic, 0)
    strong = overlap & (combined >= STRONG_MIN)
    silica_strong = overlap & (silica == 2)

    depths = {}
    for key in ("aloh", "kaolinite", "alunite"):
        depths[key] = _warp_array(
            swir_dir / f"{site.id}_{key}.tif",
            dst_t, shape, dst_crs, Resampling.average, swir_valid,
        )

    dist_m = distance_transform_edt(~strong) * PX_M
    dist_m = np.where(overlap, dist_m, np.nan)

    indices = {}
    high_counts = {}
    for key, arr in depths.items():
        indices[key] = _index_block(arr, overlap, strong, key)
        p90 = indices[key]["overlap_p90"]
        high = overlap & np.isfinite(arr) & (arr >= p90)
        high_counts[key] = int(np.count_nonzero(high))
        d_high = dist_m[high]
        d_all = dist_m[overlap]
        indices[key]["median_dist_m_high"] = float(np.nanmedian(d_high)) if d_high.size else None
        indices[key]["median_dist_m_overlap"] = float(np.nanmedian(d_all)) if d_all.size else None
        indices[key]["frac_high_within_500m"] = (
            float(np.mean(d_high <= 500.0)) if d_high.size else None
        )
        indices[key]["frac_overlap_within_500m"] = (
            float(np.mean(d_all <= 500.0)) if d_all.size else None
        )
        indices[key]["n_high"] = high_counts[key]

    silica_idx = {
        key: _index_block(arr, overlap, silica_strong, key)
        for key, arr in depths.items()
    }

    alun = depths["alunite"]
    alun_p90 = float(np.nanpercentile(alun[overlap], 90))
    high_alun = overlap & np.isfinite(alun) & (alun >= alun_p90)
    score_in = combined[high_alun].astype(np.float64)
    score_out = combined[overlap & ~high_alun].astype(np.float64)

    from pyproj import Transformer

    to_utm = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    pin_stats = {}
    for label, lonlat in (("town", pins["town"]), ("target", pins["target"])):
        x, y = to_utm.transform(*lonlat)
        # Affine inverse returns (col, row), not (row, col).
        c, r = ~dst_t * (x, y)
        r, c = int(r), int(c)
        inside = 0 <= r < shape[0] and 0 <= c < shape[1]
        pin_stats[label] = {
            "label": pins["town_label"] if label == "town" else pins["target_label"],
            "in_overlap": bool(inside and overlap[r, c]),
            "in_strong_zone": bool(inside and strong[r, c]),
            "in_silica_strong": bool(inside and silica_strong[r, c]),
            "alunite": float(alun[r, c]) if inside and np.isfinite(alun[r, c]) else None,
            "alunite_high": bool(inside and high_alun[r, c]),
        }

    stats = {
        "site_id": site.id,
        "crs": crs,
        "hypothesis": "complementary (not coincident) on the overlap",
        "universe": "90 m cells with valid TIR and valid SWIR",
        "strong_score_min": STRONG_MIN,
        "n_overlap": int(np.count_nonzero(overlap)),
        "n_strong": int(np.count_nonzero(strong)),
        "n_aloh_high": high_counts.get("aloh"),
        "n_silica_strong": int(np.count_nonzero(silica_strong)),
        "frac_strong_on_overlap": float(np.mean(strong[overlap])) if np.any(overlap) else None,
        "frac_silica_strong_on_overlap": (
            float(np.mean(silica_strong[overlap])) if np.any(overlap) else None
        ),
        "mean_combined_high_alunite": float(np.mean(score_in)) if score_in.size else None,
        "mean_combined_rest": float(np.mean(score_out)) if score_out.size else None,
        "median_combined_high_alunite": float(np.median(score_in)) if score_in.size else None,
        "median_combined_rest": float(np.median(score_out)) if score_out.size else None,
        "combined_enrichment_in_high_alunite": median_enrichment(score_in, score_out),
        "reverse_mwu": _mwu(score_in, score_out),
        "indices_in_strong": indices,
        "indices_in_silica_strong": silica_idx,
        "pins": pin_stats,
        "note": (
            "Lift = P(high SWIR | TIR group) / P(high SWIR on overlap). "
            "1 = independent.  Not Tetracorder."
        ),
    }
    arrays = {
        "overlap": overlap,
        "combined": combined,
        "strong": strong,
        "alunite": alun,
        "alunite_high": high_alun,
        "dst_t": dst_t,
        "shape": shape,
        "dst_crs": dst_crs,
        "crs": crs,
        "tir_ok": tir_ok,
        "swir_ok": swir_ok,
    }
    return stats, arrays


def save_figure(site, stats: dict, arrays: dict) -> Path:
    import geopandas as gpd
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    from matplotlib.colors import ListedColormap, PowerNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib_scalebar.scalebar import ScaleBar
    from rasterio.plot import plotting_extent

    pins = resolve_pins(site, sites_dir=ROOT / "sites")
    crs = arrays["crs"]
    rgb, bbox_extent = _load_site_basemap(site.id)
    hillshade, hs_extent = _load_hillshade(site, crs)
    town_x, town_y = _town_xy(crs, *pins["town"])
    tgt_x, tgt_y = _town_xy(crs, *pins["target"])
    extent = list(plotting_extent(
        np.empty(arrays["shape"]), transform=arrays["dst_t"],
    ))
    overlap = arrays["overlap"]
    strong = arrays["strong"]
    high = arrays["alunite_high"]
    alun = np.where(overlap, arrays["alunite"], np.nan)
    combo = np.zeros(arrays["shape"], dtype=np.uint8)
    combo[strong & ~high] = 1
    combo[high & ~strong] = 2
    combo[strong & high] = 3
    combo = np.ma.masked_where(~overlap | (combo == 0), combo)
    score = np.where(overlap, arrays["combined"].astype(np.float32), np.nan)
    layers = {
        "tir_fp": _mask_to_geom(arrays["tir_ok"], arrays["dst_t"], arrays["dst_crs"]),
        "swir_fp": _mask_to_geom(arrays["swir_ok"], arrays["dst_t"], arrays["dst_crs"]),
        "overlap": _mask_to_geom(overlap, arrays["dst_t"], arrays["dst_crs"]),
    }
    lift = stats["indices_in_strong"]["alunite"]["lift_p90"]
    jac = stats["indices_in_strong"]["alunite"]["jaccard_p90"]
    fig, axes = plt.subplots(1, 3, figsize=(14.6, 6.6), constrained_layout=True)
    titles = (
        "a  ASTER combined score\n0–6 on the overlap; strong is ≥ 3",
        "b  AVIRIS-5 alunite-like\nband depth, overlap only",
        f"c  Overlap test\nlift {lift:.2f} · Jaccard {jac:.3f}",
    )
    epsg = _epsg_label(crs)

    def _base(ax):
        ax.set_facecolor("#d0d0d0")
        ax.imshow(rgb, extent=bbox_extent, origin="upper", zorder=0,
                  interpolation="bilinear")
        if hillshade is not None:
            ax.imshow(
                _hillshade_overlay(hillshade),
                extent=hs_extent, origin="upper", zorder=1,
                interpolation="bilinear",
            )
        ax.set_xlim(bbox_extent[0], bbox_extent[1])
        ax.set_ylim(bbox_extent[2], bbox_extent[3])
        ax.set_aspect("equal")
        ax.plot(
            town_x, town_y, marker="*", color="black", markersize=11,
            markeredgecolor="white", markeredgewidth=0.5, zorder=6,
        )
        town_off = pins.get("town_offset", (-8, 8))
        ax.annotate(
            pins["town_label"], (town_x, town_y),
            textcoords="offset points", xytext=town_off,
            fontsize=8, color="#111111", ha="right",
        )
        ax.plot(
            tgt_x, tgt_y, marker="^", color="black", markersize=7,
            markeredgecolor="white", markeredgewidth=0.4, zorder=6,
        )
        tgt_off = pins.get("target_offset", (8, 6))
        ax.annotate(
            pins["target_label"], (tgt_x, tgt_y),
            textcoords="offset points", xytext=tgt_off,
            fontsize=7, color="#111111", ha="left",
        )
        ax.set_xlabel(f"Easting (km, {epsg})")
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
        )
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
        )
        ax.tick_params(labelsize=7)

    def _plot_geom(ax, geom, **kwargs):
        if geom is None or geom.is_empty:
            return
        gpd.GeoDataFrame(geometry=[geom], crs=crs).plot(ax=ax, **kwargs)

    _base(axes[0])
    score_cmap = ListedColormap([
        (1.00, 1.00, 1.00, 0.00),
        (0.96, 0.85, 0.75, 0.20),
        (0.90, 0.70, 0.55, 0.30),
        (0.75, 0.22, 0.17, 0.80),
        (0.57, 0.00, 0.00, 0.88),
        (0.40, 0.00, 0.00, 0.92),
        (0.29, 0.00, 0.00, 0.95),
    ])
    axes[0].imshow(
        np.ma.masked_invalid(score), extent=extent, origin="upper",
        cmap=score_cmap, vmin=-0.5, vmax=6.5, interpolation="nearest", zorder=2,
    )
    _plot_geom(axes[0], layers["tir_fp"], facecolor="none", edgecolor="#1a252f",
               linewidth=1.6, linestyle="--", alpha=0.85, zorder=4)

    _base(axes[1])
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad(alpha=0.0)
    im = axes[1].imshow(
        np.ma.masked_invalid(alun), extent=extent, origin="upper",
        cmap=cmap, norm=PowerNorm(gamma=0.55, vmin=0.0, vmax=0.12),
        interpolation="nearest", zorder=2,
    )
    _plot_geom(axes[1], layers["swir_fp"], facecolor="none", edgecolor="#0b5345",
               linewidth=1.6, linestyle="--", alpha=0.85, zorder=4)

    _base(axes[2])
    class_cmap = ListedColormap(["#c0392b", "#1a7a6d", "#6c3483"])
    axes[2].imshow(
        combo, extent=extent, origin="upper", interpolation="nearest",
        cmap=class_cmap, vmin=1, vmax=3, zorder=2, alpha=0.75,
    )
    _plot_geom(axes[2], layers["overlap"], facecolor="none", edgecolor="#1a252f",
               linewidth=1.2, linestyle=":", alpha=0.7, zorder=4)

    for i, ax in enumerate(axes):
        ax.set_title(titles[i], fontsize=10)
    axes[0].set_ylabel(f"Northing (km, {epsg})")
    axes[0].add_artist(ScaleBar(
        1, "m", length_fraction=0.15, location="lower left",
        box_alpha=0.7, font_properties={"size": 8},
    ))
    axes[0].annotate(
        "N\n↑", xy=(0.04, 0.96), xycoords="axes fraction",
        ha="center", va="top", fontsize=11, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                  alpha=0.8, edgecolor="gray"),
    )
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04,
                 label="Alunite-like band depth")
    fig.legend(
        handles=[
            Patch(facecolor="#c0392b", edgecolor="none", alpha=0.5,
                  label="TIR strong only"),
            Patch(facecolor="#1a7a6d", edgecolor="none", alpha=0.5,
                  label="SWIR high only (alunite p90)"),
            Patch(facecolor="#6c3483", edgecolor="none", alpha=0.6,
                  label="Both"),
            Line2D([0], [0], color="#1a252f", linestyle="--", linewidth=1.6,
                   label="TIR / overlap boundary"),
        ],
        loc="lower center", ncol=4, frameon=False, fontsize=8,
        bbox_to_anchor=(0.5, -0.04),
    )
    fig.suptitle(
        f"{site.name} · TIR strong zones vs SWIR alunite-like · "
        "overlap only · not mineral ID",
        fontsize=11,
    )
    out = ROOT / "figures" / f"gemx_{site.id}_correspondence.png"
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def summary_row_from_stats(
    stats: dict,
    *,
    status: str = "ok",
    swir_coverage: float | None = None,
) -> dict:
    """Flatten correspondence stats into a host-vs-alteration CSV row."""
    idx = stats.get("indices_in_strong") or {}
    aloh = idx.get("aloh") or {}
    kao = idx.get("kaolinite") or {}
    alun = idx.get("alunite") or {}
    j_aloh = aloh.get("jaccard_p90")
    lift_aloh = aloh.get("lift_p90")
    n_strong = stats.get("n_strong")
    n_overlap = stats.get("n_overlap")
    n_swir_high = aloh.get("n_high", stats.get("n_aloh_high"))
    row = {
        "site_id": stats.get("site_id"),
        "status": status,
        "swir_coverage": swir_coverage,
        "jaccard_aloh": j_aloh,
        "lift_aloh": lift_aloh,
        "jaccard_kaolinite": kao.get("jaccard_p90"),
        "jaccard_alunite": alun.get("jaccard_p90"),
        "n_tir_strong": n_strong,
        "n_swir_high": n_swir_high,
        "n_overlap": n_overlap,
    }
    row["verdict"] = correspondence_verdict(
        status=status,
        swir_coverage=swir_coverage,
        jaccard_aloh=j_aloh,
        lift_aloh=lift_aloh,
        n_overlap=n_overlap,
        n_tir_strong=n_strong,
    )
    return row


def run_site_correspondence(
    site_id: str,
    *,
    write_figure: bool = True,
    write_json: bool = True,
) -> dict:
    """Correspondence for one site; returns a status/verdict CSV row."""
    from critical_minerals_aster.config import load_site_by_id

    cov_map = aviris5_coverage_by_site(ROOT)
    cov = cov_map.get(site_id, 0.0)
    base = {
        "site_id": site_id,
        "status": "ok",
        "swir_coverage": cov,
        "jaccard_aloh": None,
        "lift_aloh": None,
        "jaccard_kaolinite": None,
        "jaccard_alunite": None,
        "n_tir_strong": None,
        "n_swir_high": None,
        "verdict": "no_coverage",
    }
    if cov <= 0 and not swir_tifs_exist(site_id):
        base["status"] = "no_coverage"
        base["verdict"] = "no_coverage"
        print(f"[{site_id}] no AVIRIS-5 coverage", flush=True)
        return base
    if not swir_tifs_exist(site_id):
        base["status"] = "insufficient"
        base["verdict"] = "insufficient"
        print(f"[{site_id}] SWIR GeoTIFFs missing", flush=True)
        return base

    site = load_site_by_id(site_id, ROOT / "sites")
    try:
        aster_dir = _aster_dir(site)
        b10 = aster_dir / f"{site_id}_mosaic_TIR_B10.tif"
        if not b10.is_file():
            base["status"] = "insufficient"
            base["verdict"] = "insufficient"
            print(f"[{site_id}] TIR mosaic missing ({b10})", flush=True)
            return base
        stats, arrays = compute(site)
        row = summary_row_from_stats(stats, status="ok", swir_coverage=cov)
        if write_json:
            stats_path = ROOT / "results" / f"{site_id}_correspondence.json"
            stats_path.write_text(json.dumps(stats, indent=2))
            print(f"[{site_id}] Wrote {stats_path.relative_to(ROOT)}", flush=True)
        if write_figure:
            try:
                fig_path = save_figure(site, stats, arrays)
                print(f"[{site_id}] Wrote {fig_path.relative_to(ROOT)}", flush=True)
            except FileNotFoundError as exc:
                print(f"[{site_id}] figure skipped: {exc}", flush=True)
        return row
    except Exception as exc:
        base["status"] = "error"
        base["verdict"] = "error"
        base["error"] = str(exc)
        print(f"[{site_id}] error: {exc}", flush=True)
        return base


def main() -> None:
    from critical_minerals_aster.config import load_site_by_id

    parser = argparse.ArgumentParser(
        description="ASTER strong zones vs AVIRIS-5 SWIR band-depth correspondence",
    )
    add_site_select_args(parser)
    parser.add_argument("--skip-figure", action="store_true")
    args = parser.parse_args()
    if not any([args.site, args.sites, args.primary_fdr, args.all_covered]):
        # Legacy default
        site_ids = [_site_id_from_argv("marysvale")]
    else:
        site_ids = sites_from_namespace(args, ROOT)

    rows = []
    for sid in site_ids:
        rows.append(
            run_site_correspondence(
                sid,
                write_figure=not args.skip_figure,
            )
        )
    if len(rows) == 1 and rows[0].get("status") == "ok":
        # Keep single-site console dump of full JSON when present
        site = load_site_by_id(rows[0]["site_id"], ROOT / "sites")
        stats_path = ROOT / "results" / f"{site.id}_correspondence.json"
        if stats_path.is_file():
            print(stats_path.read_text())
    elif len(rows) > 1:
        import csv

        out = ROOT / "results" / "gemx_host_vs_alteration.csv"
        fields = [
            "site_id", "status", "swir_coverage", "jaccard_aloh", "lift_aloh",
            "jaccard_kaolinite", "jaccard_alunite", "n_tir_strong", "n_swir_high",
            "verdict",
        ]
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"Wrote {out.relative_to(ROOT)} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
