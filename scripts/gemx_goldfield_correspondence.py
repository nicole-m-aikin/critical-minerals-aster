"""Goldfield overlap test: do ASTER strong zones sit on AVIRIS doublets?

Pre-specified before looking at the numbers
------------------------------------------
Universe: 90 m cells where *both* TIR and SWIR have valid pixels (the 59%
overlap).  No comparison outside that mask.

TIR "strong": pipeline combined score >= 3 (same rule as the zone polygons).
Also report strong silica (class == 2) alone.

SWIR "high": band depth >= p90 *of the overlap*, computed separately for
Al-OH, kaolinite-like, and alunite-like.  Absolute depth > 0.10 is a
sensitivity check, not the primary cut.

Hypothesis under test (user): the two are complementary, not coincident.
That predicts lift ~ 1 or < 1, low Jaccard, and high-SWIR pixels no closer
to TIR zones than the rest of the overlap.  A silica-core / clay-halo
would show low pixel overlap but *shorter* distance-to-zone than random.

This is not Tetracorder and not a drill target.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from gemx_goldfield_pilot import (  # noqa: E402
    _hillshade_overlay,
    _load_hillshade,
    _load_site_basemap,
    _town_xy,
)
from fig_goldfield_coverage_compare import (  # noqa: E402
    CRS,
    DISTRICT,
    TOWN,
    _mask_to_geom,
    _valid_mask,
    _valid_on_grid,
)

ASTER_DIR = ROOT / "data" / "sites" / "goldfield" / "aster"
SWIR_DIR = ROOT / "data" / "gemx" / "aviris" / "goldfield"
STRONG_MIN = 3
PX_M = 90.0


def median_enrichment(inside: np.ndarray, outside: np.ndarray) -> float:
    """median(inside) / median(outside); NaN if either side is empty or 0."""
    if inside.size == 0 or outside.size == 0:
        return float("nan")
    den = float(np.median(outside))
    if den == 0.0:
        return float("nan")
    return float(np.median(inside) / den)


def binary_lift(event: np.ndarray, group: np.ndarray) -> float:
    """P(event|group) / P(event) on the same universe.  1 = independent."""
    if event.size == 0 or not np.any(group):
        return float("nan")
    base = float(np.mean(event))
    if base == 0.0:
        return float("nan")
    return float(np.mean(event[group]) / base)


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    both = np.count_nonzero(a & b)
    union = np.count_nonzero(a | b)
    if union == 0:
        return float("nan")
    return both / union


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """+1 if every x > every y; 0 if indistinguishable.  Uses MW U."""
    from scipy.stats import mannwhitneyu

    if x.size == 0 or y.size == 0:
        return float("nan")
    u = mannwhitneyu(x, y, alternative="two-sided").statistic
    return float((2.0 * u) / (x.size * y.size) - 1.0)


def _site_grid(site):
    from rasterio.crs import CRS as RasterioCRS
    from rasterio.transform import from_bounds
    from rasterio.warp import transform_bounds

    west, south, east, north = transform_bounds(
        "EPSG:4326", CRS, *site.bbox_wgs84, densify_pts=21,
    )
    cols = max(1, int(round((east - west) / PX_M)))
    rows = max(1, int(round((north - south) / PX_M)))
    dst_t = from_bounds(west, south, east, north, cols, rows)
    return dst_t, (rows, cols), RasterioCRS.from_epsg(32611)


def _warp_array(path: Path, dst_t, shape, dst_crs, resampling, valid_fn=None):
    import rasterio
    from rasterio.warp import reproject

    rows, cols = shape
    with rasterio.open(path) as src:
        src_a = src.read(1).astype(np.float32)
        if valid_fn is not None:
            src_a = np.where(valid_fn(src_a, src.nodata), src_a, np.nan)
        dst = np.full((rows, cols), np.nan, dtype=np.float32)
        reproject(
            src_a, dst,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=dst_t, dst_crs=dst_crs,
            resampling=resampling,
            src_nodata=np.nan,
            dst_nodata=np.nan,
        )
    return dst


def _warp_class(path: Path, dst_t, shape, dst_crs):
    """Nearest-neighbour class warp.  0 is a real class; do not treat as nodata."""
    import rasterio
    from rasterio.warp import reproject, Resampling

    rows, cols = shape
    with rasterio.open(path) as src:
        src_a = src.read(1).astype(np.float32)
        dst = np.zeros((rows, cols), dtype=np.float32)
        reproject(
            src_a, dst,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=dst_t, dst_crs=dst_crs,
            resampling=Resampling.nearest,
        )
    return dst.astype(np.uint8)


def _mwu(x, y) -> dict:
    from scipy.stats import mannwhitneyu

    if x.size == 0 or y.size == 0:
        return {"p_two_sided": None, "cliffs_delta": None, "n_in": int(x.size), "n_out": int(y.size)}
    res = mannwhitneyu(x, y, alternative="two-sided")
    return {
        "p_two_sided": float(res.pvalue),
        "cliffs_delta": cliffs_delta(x, y),
        "n_in": int(x.size),
        "n_out": int(y.size),
    }


def _index_block(depth: np.ndarray, overlap: np.ndarray, group: np.ndarray, name: str) -> dict:
    d = depth[overlap]
    g = group[overlap]
    finite = np.isfinite(d)
    inside = d[g & finite]
    outside = d[~g & finite]
    p90 = float(np.nanpercentile(d, 90))
    high = np.isfinite(d) & (d >= p90)
    deep = np.isfinite(d) & (d > 0.10)
    mwu = _mwu(inside, outside)
    return {
        "name": name,
        "n_finite_in": int(inside.size),
        "n_finite_out": int(outside.size),
        "p50_in": float(np.median(inside)) if inside.size else None,
        "p50_out": float(np.median(outside)) if outside.size else None,
        "p90_in": float(np.percentile(inside, 90)) if inside.size else None,
        "p90_out": float(np.percentile(outside, 90)) if outside.size else None,
        "median_enrichment": median_enrichment(inside, outside),
        "overlap_p90": p90,
        "frac_gt_0_10_in": float(np.mean(inside > 0.10)) if inside.size else None,
        "frac_gt_0_10_out": float(np.mean(outside > 0.10)) if outside.size else None,
        "lift_p90": binary_lift(high, g),
        "lift_gt_0_10": binary_lift(deep, g),
        "jaccard_p90": jaccard(g, high),
        "frac_high_inside_group": float(np.mean(g[high])) if np.any(high) else None,
        "frac_group_on_overlap": float(np.mean(g)),
        **mwu,
    }


def compute(site) -> tuple[dict, dict]:
    from rasterio.warp import Resampling
    from scipy.ndimage import distance_transform_edt

    dst_t, shape, dst_crs = _site_grid(site)

    def tir_valid(a, nodata):
        return _valid_mask(a, nodata) & (a > 0)

    def swir_valid(a, nodata):
        return _valid_mask(a, nodata)

    tir_ok = _valid_on_grid(
        ASTER_DIR / "goldfield_mosaic_TIR_B10.tif",
        tir_valid, dst_t, shape, dst_crs,
    )
    swir_ok = _valid_on_grid(
        SWIR_DIR / "goldfield_aloh.tif",
        swir_valid, dst_t, shape, dst_crs,
    )
    cls_ok = _valid_on_grid(
        ASTER_DIR / "goldfield_mosaic_cls_valid.tif",
        lambda a, nd: a > 0, dst_t, shape, dst_crs,
    )
    overlap = tir_ok & swir_ok & cls_ok

    silica = _warp_class(ASTER_DIR / "goldfield_mosaic_cls_silica.tif", dst_t, shape, dst_crs)
    carb = _warp_class(ASTER_DIR / "goldfield_mosaic_cls_carbonate.tif", dst_t, shape, dst_crs)
    mafic = _warp_class(ASTER_DIR / "goldfield_mosaic_cls_mafic.tif", dst_t, shape, dst_crs)
    combined = np.where(cls_ok, silica.astype(np.int16) + carb + mafic, 0)
    strong = overlap & (combined >= STRONG_MIN)
    silica_strong = overlap & (silica == 2)

    depths = {}
    for key in ("aloh", "kaolinite", "alunite"):
        depths[key] = _warp_array(
            SWIR_DIR / f"goldfield_{key}.tif",
            dst_t, shape, dst_crs, Resampling.average, swir_valid,
        )

    dist_m = distance_transform_edt(~strong) * PX_M
    dist_m = np.where(overlap, dist_m, np.nan)

    indices = {}
    for key, arr in depths.items():
        indices[key] = _index_block(arr, overlap, strong, key)
        p90 = indices[key]["overlap_p90"]
        high = overlap & np.isfinite(arr) & (arr >= p90)
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

    silica_idx = {
        key: _index_block(arr, overlap, silica_strong, key)
        for key, arr in depths.items()
    }

    # Reverse: combined score inside high-alunite vs rest of overlap.
    alun = depths["alunite"]
    alun_p90 = float(np.nanpercentile(alun[overlap], 90))
    high_alun = overlap & np.isfinite(alun) & (alun >= alun_p90)
    score_in = combined[high_alun].astype(np.float64)
    score_out = combined[overlap & ~high_alun].astype(np.float64)

    from pyproj import Transformer

    to_utm = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    pins = {}
    for label, lonlat in (("town", TOWN), ("district", DISTRICT)):
        x, y = to_utm.transform(*lonlat)
        # Affine inverse returns (col, row), not (row, col).
        c, r = ~dst_t * (x, y)
        r, c = int(r), int(c)
        inside = 0 <= r < shape[0] and 0 <= c < shape[1]
        pins[label] = {
            "in_overlap": bool(inside and overlap[r, c]),
            "in_strong_zone": bool(inside and strong[r, c]),
            "in_silica_strong": bool(inside and silica_strong[r, c]),
            "alunite": float(alun[r, c]) if inside and np.isfinite(alun[r, c]) else None,
            "alunite_high": bool(inside and high_alun[r, c]),
        }

    stats = {
        "site_id": site.id,
        "hypothesis": "complementary (not coincident) on the overlap",
        "universe": "90 m cells with valid TIR and valid SWIR",
        "strong_score_min": STRONG_MIN,
        "n_overlap": int(np.count_nonzero(overlap)),
        "n_strong": int(np.count_nonzero(strong)),
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
        "pins": pins,
        "note": (
            "Lift = P(high SWIR | TIR group) / P(high SWIR on overlap). "
            "1 = independent.  Cliff's delta > 0 means SWIR is deeper inside "
            "the TIR group.  Not Tetracorder."
        ),
    }
    arrays = {
        "overlap": overlap,
        "combined": combined,
        "strong": strong,
        "silica_strong": silica_strong,
        "alunite": alun,
        "alunite_high": high_alun,
        "kaolinite": depths["kaolinite"],
        "dst_t": dst_t,
        "shape": shape,
        "dst_crs": dst_crs,
        "tir_ok": tir_ok,
        "swir_ok": swir_ok,
    }
    return stats, arrays


def save_figure(site, stats: dict, arrays: dict) -> Path:
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    from matplotlib.colors import ListedColormap, PowerNorm
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib_scalebar.scalebar import ScaleBar
    from rasterio.plot import plotting_extent

    rgb, bbox_extent = _load_site_basemap()
    hillshade, hs_extent = _load_hillshade(site, CRS)
    town_x, town_y = _town_xy(CRS, *TOWN)
    dist_x, dist_y = _town_xy(CRS, *DISTRICT)
    extent = list(plotting_extent(
        np.empty(arrays["shape"]), transform=arrays["dst_t"],
    ))
    overlap = arrays["overlap"]
    strong = arrays["strong"]
    high = arrays["alunite_high"]
    alun = np.where(overlap, arrays["alunite"], np.nan)
    # 0 empty / 1 TIR only / 2 SWIR only / 3 both — drawn as a raster so
    # single 90 m zone cells are not deleted by polygon simplify.
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
        ax.annotate(
            "Goldfield", (town_x, town_y),
            textcoords="offset points", xytext=(-8, 8),
            fontsize=8, color="#111111", ha="right",
        )
        ax.plot(
            dist_x, dist_y, marker="^", color="black", markersize=7,
            markeredgecolor="white", markeredgewidth=0.4, zorder=6,
        )
        ax.set_xlabel("Easting (km, EPSG:32611)")
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
        )
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
        )
        ax.tick_params(labelsize=7)

    import geopandas as gpd

    def _plot_geom(ax, geom, **kwargs):
        if geom is None or geom.is_empty:
            return
        gpd.GeoDataFrame(geometry=[geom], crs=CRS).plot(ax=ax, **kwargs)

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
    axes[0].set_ylabel("Northing (km, EPSG:32611)")
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
        "Goldfield–Cuprite · TIR strong zones vs SWIR alunite-like · "
        "overlap only · not mineral ID",
        fontsize=11,
    )
    out = ROOT / "figures" / "gemx_goldfield_correspondence.png"
    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def main() -> None:
    from critical_minerals_aster.config import load_site_by_id

    site = load_site_by_id("goldfield", ROOT / "sites")
    stats, arrays = compute(site)
    stats_path = ROOT / "results" / "goldfield_correspondence.json"
    stats_path.write_text(json.dumps(stats, indent=2))
    fig_path = save_figure(site, stats, arrays)
    print(json.dumps(stats, indent=2))
    print(f"Wrote {stats_path.relative_to(ROOT)}")
    print(f"Wrote {fig_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
