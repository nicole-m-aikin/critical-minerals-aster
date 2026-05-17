"""End-to-end per-site processing (classification, vectors, metrics, provenance)."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from critical_minerals_aster.classification import (
    classify_percentiles,
    combined_score,
    vectorize_strong_zones,
)
from critical_minerals_aster.config import BBox, SiteConfig, search_bbox
from critical_minerals_aster.metrics import compute_site_summary, write_site_summary
from critical_minerals_aster.paths import SitePaths, site_paths_for
from critical_minerals_aster.spectral import (
    alteration_ratios,
    clip_bands_to_bbox,
    extract_granule_id,
    load_tir_bands_10_14,
    raster_bbox_wgs84,
    select_granule,
)
from critical_minerals_aster.structure import (
    annotate_deposits_with_structure,
    load_structure_layers,
    structure_buffer_union,
)

_GRANULE_ID_RE = re.compile(r"(AST_L1T_\d+_\d+)")


def _bbox_annotation(site: "SiteConfig") -> str:
    """Return a compact geographic extent string for figure annotations."""
    west, south, east, north = site.bbox_wgs84
    ew_dir = "W" if west < 0 else "E"
    return (
        f"{south:.2f}°–{north:.2f}°N · "
        f"{abs(west):.2f}°–{abs(east):.2f}°{ew_dir}"
    )


def resolve_granule_id(site: SiteConfig, paths: SitePaths) -> str:
    if site.granule_id:
        return site.granule_id
    aster_dir = paths.aster_dir
    if not aster_dir.is_dir():
        raise FileNotFoundError(f"ASTER directory not found: {aster_dir}")
    # Check for a pre-built mosaic before scanning for individual granule files.
    mosaic_b10 = aster_dir / f"{site.id}_mosaic_TIR_B10.tif"
    if mosaic_b10.is_file():
        return f"{site.id}_mosaic"
    for name in sorted(aster_dir.iterdir()):
        match = _GRANULE_ID_RE.search(name.name)
        if match and "TIR_B10" in name.name:
            return match.group(1)
    raise ValueError(f"No ASTER granule TIR files found under {aster_dir}")


def run_classification(
    site: SiteConfig, paths: SitePaths, granule_id: str
) -> tuple[
    gpd.GeoDataFrame,
    np.ndarray,  # silica
    np.ndarray,  # carbonate
    np.ndarray,  # mafic
    np.ndarray,  # silica_cls
    np.ndarray,  # carbonate_cls
    np.ndarray,  # mafic_cls
    np.ndarray,  # combined
    BBox,        # raster_bbox — WGS84 extent of the analysed (clipped) raster
    Any,         # transform — affine transform of the clipped raster
    tuple[int, int],  # shape — (rows, cols) of the clipped raster
    Any,         # crs — coordinate reference system of the raster
]:
    """Classify, vectorize, return zones, class maps, raster extent, and raster metadata.

    The raster_bbox is the WGS84 bounding box of the ASTER data *actually
    analysed* (i.e. after bbox clipping).  It is the intersection of the ASTER
    granule footprint and site.bbox_wgs84, so MRDS deposit queries should use
    it rather than the raw site bbox to avoid counting deposits that fall
    outside the TIR coverage.

    Elements 10–12 (transform, shape, crs) describe the clipped raster pixel
    grid so callers can reproject auxiliary data (e.g. a DEM hillshade) to
    exactly the same extent.
    """
    _, _, b12, b13, b14, _, transform, crs = load_tir_bands_10_14(
        paths.aster_dir, granule_id
    )
    # Clip to site bbox so percentile thresholds and zone counts are
    # site-specific rather than whole-scene artifacts.  Shared-granule sites
    # (e.g. goldfield/silver_peak on the same ASTER swath) would otherwise
    # produce identical zone polygons from the full 60-90 km scene.
    (b12, b13, b14), transform = clip_bands_to_bbox(
        [b12, b13, b14], transform, crs, site.bbox_wgs84
    )
    # Record the actual raster extent AFTER clipping so downstream MRDS
    # queries are constrained to the true TIR coverage area, not the full
    # (possibly larger) site.bbox_wgs84.
    raster_bbox: BBox = raster_bbox_wgs84(transform, b12.shape, crs)

    silica, carbonate, mafic = alteration_ratios(b12, b13, b14)

    cp = site.classification
    assert cp is not None
    silica_cls, _, _ = classify_percentiles(silica, cp.low_pct, cp.high_pct)
    carbonate_cls, _, _ = classify_percentiles(carbonate, cp.low_pct, cp.high_pct)
    mafic_cls, _, _ = classify_percentiles(mafic, cp.low_pct, cp.high_pct)
    combined = combined_score(silica_cls, carbonate_cls, mafic_cls)

    zones = vectorize_strong_zones(
        combined, transform, crs, min_score=cp.strong_score_min
    )
    return (
        zones,
        silica,
        carbonate,
        mafic,
        silica_cls,
        carbonate_cls,
        mafic_cls,
        combined,
        raster_bbox,
        transform,
        b12.shape,
        crs,
    )


def compute_global_limits(
    site_ids: list[str],
    repo_root: Path,
    *,
    low_pct: float = 2,
    high_pct: float = 98,
    subsample: int = 10,
) -> dict[str, tuple[float, float]]:
    """Compute cross-site percentile limits for each band ratio.

    Loads TIR bands for every site in *site_ids*, computes the three
    alteration ratios, and collects every *subsample*-th finite pixel from
    all sites.  Returns 2nd–98th percentile limits suitable for passing as
    the *global_limits* argument to :func:`save_band_ratio_figure`.

    Parameters
    ----------
    site_ids:
        List of site IDs to include (must have ASTER data on disk).
    repo_root:
        Repository root path.
    low_pct / high_pct:
        Percentile bounds (default 2 / 98).
    subsample:
        Take every N-th pixel to keep memory usage bounded.
    """
    from critical_minerals_aster.config import load_site_by_id

    sites_dir = repo_root / "sites"
    all_silica: list[np.ndarray] = []
    all_carbonate: list[np.ndarray] = []
    all_mafic: list[np.ndarray] = []

    for site_id in site_ids:
        try:
            site = load_site_by_id(site_id, sites_dir)
            paths = site_paths_for(site, repo_root)
            granule_id = resolve_granule_id(site, paths)
            _, _, b12, b13, b14, _, transform, crs = load_tir_bands_10_14(
                paths.aster_dir, granule_id
            )
            (b12, b13, b14), _ = clip_bands_to_bbox(
                [b12, b13, b14], transform, crs, site.bbox_wgs84
            )
            silica, carbonate, mafic = alteration_ratios(b12, b13, b14)
            flat = silica.ravel()
            all_silica.append(flat[np.isfinite(flat)][::subsample])
            flat = carbonate.ravel()
            all_carbonate.append(flat[np.isfinite(flat)][::subsample])
            flat = mafic.ravel()
            all_mafic.append(flat[np.isfinite(flat)][::subsample])
        except Exception as exc:
            print(f"  [global_limits] skipping {site_id}: {exc}", file=sys.stderr)

    def _limits(arrays: list[np.ndarray]) -> tuple[float, float]:
        combined = np.concatenate(arrays) if arrays else np.array([])
        if combined.size == 0:
            return (0.0, 1.0)
        return (
            float(np.percentile(combined, low_pct)),
            float(np.percentile(combined, high_pct)),
        )

    return {
        "silica": _limits(all_silica),
        "carbonate": _limits(all_carbonate),
        "mafic": _limits(all_mafic),
    }


def save_band_ratio_figure(
    site: SiteConfig,
    paths: SitePaths,
    silica: np.ndarray,
    carbonate: np.ndarray,
    mafic: np.ndarray,
    hillshade: np.ndarray | None = None,
    global_limits: dict[str, tuple[float, float]] | None = None,
) -> None:
    """Save Figure 01 — TIR band ratio panels.

    Parameters
    ----------
    global_limits:
        Optional dict with keys ``"silica"``, ``"carbonate"``, ``"mafic"``
        mapping to ``(vmin, vmax)`` tuples.  When supplied the same colorbar
        range is used for every site, making cross-site comparisons valid.
        When *None* (default) per-site 2nd–98th percentile limits are used.
    """
    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    ratios = [
        (silica, "Silica/quartz (B13/B14)", "magma", "silica"),
        (carbonate, "Carbonate/dolomite (B13/B12)", "YlOrBr", "carbonate"),
        (mafic, "Mafic minerals (B12/B13)", "PuBu", "mafic"),
    ]
    for ax, (ratio, title, cmap, key) in zip(axes, ratios):
        if global_limits is not None and key in global_limits:
            vmin, vmax = global_limits[key]
        else:
            vmin, vmax = _percentile_limits(ratio, 2, 98)
        im = ax.imshow(
            ratio,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        if hillshade is not None:
            _hs_cmap = plt.cm.gray.copy()
            _hs_cmap.set_bad(alpha=0.0)
            ax.imshow(hillshade, cmap=_hs_cmap, alpha=0.25, vmin=0, vmax=1)
        ax.set_title(title)
        ax.axis("off")
        plt.colorbar(im, ax=ax, shrink=0.8)
    plt.suptitle(f"ASTER TIR Band Ratios — {site.name}", fontsize=11, y=1.02)
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    fig.text(
        0.5, 0.01, _bbox_annotation(site),
        ha="center", fontsize=7.5, color="#555555",
    )
    plt.savefig(
        paths.figures_dir / "01_tir_band_ratios.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig)


def _percentile_limits(
    ratio: np.ndarray,
    low_pct: float,
    high_pct: float,
) -> tuple[float | None, float | None]:
    finite = ratio[np.isfinite(ratio)]
    if finite.size == 0:
        return None, None
    return (
        float(np.percentile(finite, low_pct)),
        float(np.percentile(finite, high_pct)),
    )


def _normalize_ratio_channel(
    ratio: np.ndarray,
    low_pct: float,
    high_pct: float,
    scale: float = 1.0,
) -> np.ndarray:
    p_low, p_high = _percentile_limits(ratio, low_pct, high_pct)
    if p_low is None or p_high is None or p_high == p_low:
        return np.zeros_like(ratio, dtype=float)
    return np.clip((ratio - p_low) / (p_high - p_low), 0, 1) * scale


def save_composite_figure(
    site: SiteConfig,
    paths: SitePaths,
    silica: np.ndarray,
    carbonate: np.ndarray,
    mafic: np.ndarray,
) -> None:
    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    rgb = np.dstack(
        [
            _normalize_ratio_channel(silica, 20, 80, scale=0.6),
            _normalize_ratio_channel(carbonate, 2, 98),
            _normalize_ratio_channel(mafic, 20, 80, scale=0.8),
        ]
    )
    nan_mask = np.isnan(silica) | np.isnan(carbonate) | np.isnan(mafic)
    rgb[nan_mask] = 0

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(rgb)
    ax.set_title(
        f"False-color composite — {site.name}\n"
        "Red=silica · Green=carbonate · Blue=mafic",
        fontsize=12,
    )
    ax.axis("off")
    fig.text(
        0.5, 0.01, _bbox_annotation(site),
        ha="center", fontsize=7.5, color="#555555",
    )
    plt.tight_layout()
    plt.savefig(
        paths.figures_dir / "00_composite_rgb.png",
        dpi=150,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_classification_figure(
    site: SiteConfig,
    paths: SitePaths,
    silica_cls: np.ndarray,
    carbonate_cls: np.ndarray,
    mafic_cls: np.ndarray,
    combined: np.ndarray,
    hillshade: np.ndarray | None = None,
) -> None:
    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    arrays = [silica_cls, carbonate_cls, mafic_cls, combined]
    titles = ["Silica classes", "Carbonate classes", "Mafic classes", "Combined score"]
    for ax, arr, title in zip(axes, arrays, titles):
        im = ax.imshow(arr, cmap="YlOrRd" if title != "Combined score" else "viridis")
        if hillshade is not None:
            _hs_cmap = plt.cm.gray.copy()
            _hs_cmap.set_bad(alpha=0.0)
            ax.imshow(hillshade, cmap=_hs_cmap, alpha=0.25, vmin=0, vmax=1)
        ax.set_title(title)
        ax.axis("off")
        plt.colorbar(im, ax=ax, shrink=0.8)
    plt.suptitle(f"Alteration classification — {site.name}", fontsize=13)
    plt.tight_layout()
    fig.text(
        0.5, 0.01, _bbox_annotation(site),
        ha="center", fontsize=7.5, color="#555555",
    )
    plt.savefig(paths.figures_dir / "02_classification.png", dpi=150, bbox_inches="tight")
    plt.close()


def save_deposit_overlay_figure(
    site: SiteConfig,
    paths: SitePaths,
    zones: gpd.GeoDataFrame,
    deposits: gpd.GeoDataFrame,
    repo_root: Path,
    hillshade: np.ndarray | None = None,
    raster_transform: "rasterio.Affine | None" = None,
    raster_shape: "tuple[int, int] | None" = None,
    structs: "gpd.GeoDataFrame | None" = None,
    n_on_structure: "int | None" = None,
    n_total_deposits: "int | None" = None,
    hs_transform: "rasterio.Affine | None" = None,
    hs_shape: "tuple[int, int] | None" = None,
) -> None:
    """Figure 03 — spatial map of strong anomaly zones with MRDS deposit overlay."""
    import matplotlib.lines as mlines
    import matplotlib.ticker

    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 10))

    if len(zones) > 0:
        small = zones[zones["area_km2"] < 10]
        large = zones[zones["area_km2"] >= 10]
        if len(small):
            small.plot(ax=ax, color="#c0392b", alpha=0.55, linewidth=0)
        if len(large):
            large.plot(ax=ax, color="#4a0000", alpha=0.90, linewidth=0)

    # Draw structural corridor buffer fill before deposit points.
    has_structure = False
    _buffer_m: float = 500.0
    if site.structure_layers:
        if structs is None:
            target_crs = zones.crs if len(zones) else deposits.crs
            structs = load_structure_layers(site, repo_root, target_crs)
        _buffer_m = site.structure_layers[0].buffer_m
        if structs is not None and not structs.empty:
            union_geom = structure_buffer_union(structs, _buffer_m)
            if union_geom is not None:
                gpd.GeoSeries([union_geom], crs=structs.crs).plot(
                    ax=ax, color="#e67e22", alpha=0.25, linewidth=0, zorder=1
                )
                has_structure = True

    outside = deposits[~deposits["inside_zone"]] if "inside_zone" in deposits.columns else deposits
    inside = deposits[deposits["inside_zone"]] if "inside_zone" in deposits.columns else gpd.GeoDataFrame()

    # Dense-deposit handling: reduce marker size/alpha for large sets; cap display at 500.
    _n_outside = len(outside)
    _outside_capped_note = ""
    if _n_outside > 500:
        outside = outside.iloc[:500]
        _outside_capped_note = f" (showing 500 of {_n_outside})"
    _outside_ms = 8 if _n_outside > 100 else 20
    _outside_alpha = 0.5 if _n_outside > 100 else 0.8

    if len(outside):
        outside.plot(ax=ax, color="steelblue", markersize=_outside_ms, alpha=_outside_alpha, zorder=3)
    if len(inside):
        inside.plot(ax=ax, color="gold", markersize=70, marker="*", zorder=4,
                    edgecolor="black", linewidth=0.5)

    # Compute view bounds from the hillshade grid (covers full site.bbox_wgs84)
    # so axes limits always match the configured site area.  Fall back to the
    # ASTER raster extent, then zone bounds if neither is available.
    _limit_transform = hs_transform if hs_transform is not None else raster_transform
    _limit_shape = hs_shape if hs_shape is not None else raster_shape
    if _limit_transform is not None and _limit_shape is not None:
        _r = _limit_transform
        _rc, _cc = _limit_shape
        _rx0, _rx1 = _r.c, _r.c + _r.a * _cc
        _ry0, _ry1 = _r.f + _r.e * _rc, _r.f
        _mx = (_rx1 - _rx0) * 0.01
        _my = (_ry1 - _ry0) * 0.01
        _xlim: tuple[float, float] = (_rx0 - _mx, _rx1 + _mx)
        _ylim: tuple[float, float] = (_ry0 - _my, _ry1 + _my)
    elif len(zones) > 0:
        _zb = zones.total_bounds  # xmin, ymin, xmax, ymax
        _zm = max(_zb[2] - _zb[0], _zb[3] - _zb[1]) * 0.02
        _xlim = (_zb[0] - _zm, _zb[2] + _zm)
        _ylim = (_zb[1] - _zm, _zb[3] + _zm)
    else:
        _xlim = ax.get_xlim()
        _ylim = ax.get_ylim()

    # Light neutral background — visible where ASTER doesn't cover the bbox
    # (swath rotation gaps) and at the axes margin.  #f0f0f0 means even the
    # darkest hillshade shadow (alpha-blended over this) stays above medium gray.
    ax.set_facecolor("#f0f0f0")

    # Hillshade now covers site.bbox_wgs84 (full configured site area).
    # Use hs_transform/hs_shape for extent; fall back to raster grid if absent.
    _hs_t = hs_transform if hs_transform is not None else raster_transform
    _hs_s = hs_shape if hs_shape is not None else raster_shape
    if hillshade is not None and _hs_t is not None and _hs_s is not None:
        rows, cols = _hs_s
        t = _hs_t
        _hs_extent = (t.c, t.c + t.a * cols, t.f + t.e * rows, t.f)
        _hs_cmap = plt.cm.gray.copy()
        _hs_cmap.set_bad(alpha=0.0)  # NaN nodata pixels → show background color
        ax.imshow(
            hillshade, cmap=_hs_cmap, alpha=0.35, vmin=0, vmax=1,
            extent=_hs_extent, origin="upper", zorder=0,
        )

    # Restore limits — hillshade imshow may have expanded the view to its own extent.
    ax.set_xlim(_xlim)
    ax.set_ylim(_ylim)

    # Format UTM tick labels as km integers for readability.
    km_fmt = matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
    ax.xaxis.set_major_formatter(km_fmt)
    ax.yaxis.set_major_formatter(km_fmt)

    _crs_label = zones.crs if len(zones) else (deposits.crs if len(deposits) else None)
    _epsg = _crs_label.to_epsg() if _crs_label is not None else None
    _epsg_suffix = f" (EPSG:{_epsg})" if _epsg else ""
    ax.set_xlabel(f"Easting (km{_epsg_suffix})")
    ax.set_ylabel("Northing (km)")

    legend_elements = [
        mlines.Line2D([], [], color="#c0392b", linewidth=6, alpha=0.5, label="Strong anomaly zone (< 10 km²)"),
        mlines.Line2D([], [], color="#4a0000", linewidth=6, alpha=0.9, label="Strong anomaly zone (≥ 10 km²)"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor="steelblue",
                      markersize=9, label=f"MRDS deposit (outside zone, n={_n_outside}{_outside_capped_note})"),
        mlines.Line2D([], [], marker="*", color="w", markerfacecolor="gold",
                      markeredgecolor="black", markersize=13,
                      label=f"MRDS deposit (inside zone, n={len(inside)})"),
    ]
    if has_structure:
        legend_elements.append(
            mlines.Line2D(
                [], [], color="#e67e22", linewidth=6, alpha=0.4,
                label=f"Structural corridor (±{_buffer_m:.0f} m)",
            )
        )
    ax.legend(handles=legend_elements, loc="upper right", framealpha=0.95, fontsize=9)

    # Structure fraction annotation — placed above the scale bar to avoid overlap.
    # y=0.09 clears the matplotlib-scalebar "lower left" widget (~0.02–0.07 height).
    if has_structure and n_on_structure is not None and n_total_deposits is not None:
        _pct = (n_on_structure / n_total_deposits * 100) if n_total_deposits > 0 else 0.0
        ax.text(
            0.02, 0.09,
            f"{n_on_structure}/{n_total_deposits} deposits ({_pct:.0f}%) within {_buffer_m:.0f} m of structure",
            transform=ax.transAxes,
            fontsize=8.5,
            va="bottom",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="none"),
            zorder=10,
        )

    ax.set_title(
        f"Strong alteration zones vs MRDS deposits\n{site.name}\n{_bbox_annotation(site)}",
        fontsize=12,
    )

    # Scale bar using matplotlib-scalebar.
    try:
        from matplotlib_scalebar.scalebar import ScaleBar
        ax.add_artist(ScaleBar(1, "m", length_fraction=0.2, location="lower left",
                               box_alpha=0.7, font_properties={"size": 9}))
    except Exception:
        # Fallback: manual 10-km line scale bar in the lower-left corner.
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        bar_len_m = 10_000
        bar_x0 = xlim[0] + (xlim[1] - xlim[0]) * 0.05
        bar_y = ylim[0] + (ylim[1] - ylim[0]) * 0.04
        ax.plot([bar_x0, bar_x0 + bar_len_m], [bar_y, bar_y],
                color="black", linewidth=3, solid_capstyle="butt", zorder=10)
        ax.text(bar_x0 + bar_len_m / 2, bar_y + (ylim[1] - ylim[0]) * 0.015,
                "10 km", ha="center", va="bottom", fontsize=8, zorder=10)

    # North arrow — simple text annotation in the upper-left corner.
    ax.annotate(
        "N\n↑",
        xy=(0.04, 0.96),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=13,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="gray"),
    )

    plt.tight_layout()
    plt.savefig(paths.figures_dir / "03_deposit_overlay.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_structure_proximity_figure(
    site: SiteConfig,
    paths: SitePaths,
    annotated: gpd.GeoDataFrame,
) -> None:
    """Figure 05b — strip chart of deposit-to-structure distances by commodity group.

    One horizontal strip per commodity group.  Each deposit is a point at its
    ``nearest_structure_m`` value (log x-axis if the range spans > 10×).
    Points are coloured gold (inside anomaly zone) or steelblue (outside).
    A vertical dashed line marks the ``buffer_m`` threshold from the first
    structure layer.

    Only written when the site has structure data and ``nearest_structure_m``
    is present in *annotated*.
    """
    if "nearest_structure_m" not in annotated.columns:
        return
    if "commodity_group" not in annotated.columns:
        return
    if annotated["nearest_structure_m"].isna().all():
        return

    buffer_m = site.structure_layers[0].buffer_m if site.structure_layers else 500.0

    groups = (
        annotated.groupby("commodity_group")["nearest_structure_m"]
        .median()
        .sort_values()
        .index.tolist()
    )
    # Drop groups with no structure data at all.
    groups = [g for g in groups if not annotated.loc[annotated["commodity_group"] == g, "nearest_structure_m"].isna().all()]
    if not groups:
        return

    n_groups = len(groups)
    fig, ax = plt.subplots(figsize=(10, max(3, n_groups * 0.65 + 1.5)))

    # Use log scale when range spans more than 10× — typical for fault-distance data.
    valid = annotated["nearest_structure_m"].dropna()
    use_log = (valid.max() / max(valid.min(), 1)) > 10

    # Jitter positions on y so overlapping points are readable.
    rng = np.random.default_rng(42)
    y_positions = {grp: i for i, grp in enumerate(groups)}

    for grp in groups:
        sub = annotated[annotated["commodity_group"] == grp].copy()
        sub = sub.dropna(subset=["nearest_structure_m"])
        if sub.empty:
            continue
        y_base = y_positions[grp]
        jitter = rng.uniform(-0.25, 0.25, size=len(sub))

        inside_mask = sub.get("inside_zone", pd.Series(False, index=sub.index)).fillna(False)
        ax.scatter(
            sub.loc[~inside_mask, "nearest_structure_m"],
            y_base + jitter[~inside_mask.values],
            color="steelblue", alpha=0.7, s=18, linewidths=0, zorder=3,
        )
        ax.scatter(
            sub.loc[inside_mask, "nearest_structure_m"],
            y_base + jitter[inside_mask.values],
            color="gold", alpha=0.9, s=30, marker="*",
            edgecolors="black", linewidths=0.4, zorder=4,
        )

    ax.axvline(buffer_m, color="#e67e22", linestyle="--", linewidth=1.5,
               alpha=0.8, label=f"Buffer threshold ({buffer_m:.0f} m)", zorder=2)

    ax.set_yticks(range(n_groups))
    ax.set_yticklabels(groups, fontsize=9)
    ax.set_ylim(-0.6, n_groups - 0.4)

    if use_log:
        ax.set_xscale("log")
        ax.set_xlabel("Distance to nearest structure (m, log scale)", fontsize=10)
    else:
        ax.set_xlabel("Distance to nearest structure (m)", fontsize=10)

    import matplotlib.lines as mlines
    legend_elements = [
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor="steelblue",
                      markersize=7, label="Outside anomaly zone"),
        mlines.Line2D([], [], marker="*", color="w", markerfacecolor="gold",
                      markeredgecolor="black", markersize=10, label="Inside anomaly zone"),
        mlines.Line2D([], [], color="#e67e22", linestyle="--", linewidth=1.5,
                      label=f"Buffer threshold ({buffer_m:.0f} m)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8.5, framealpha=0.9)
    ax.grid(axis="x", alpha=0.3, zorder=0)

    ax.set_title(
        f"Structural proximity by commodity group\n{site.name}", fontsize=12
    )
    plt.tight_layout()

    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    out = paths.figures_dir / "05_structure_proximity.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_commodity_correlation_figure(
    site: SiteConfig,
    paths: SitePaths,
    deposits: gpd.GeoDataFrame,
) -> None:
    """Figure 04 — horizontal stacked bar chart of hit rate by commodity group."""
    if "inside_zone" not in deposits.columns or "commodity_group" not in deposits.columns:
        return
    if len(deposits) == 0:
        return

    paths.figures_dir.mkdir(parents=True, exist_ok=True)

    ct = pd.crosstab(deposits["commodity_group"], deposits["inside_zone"])
    ct.columns = [str(c) for c in ct.columns]
    outside_col = "False" if "False" in ct.columns else ct.columns[0]
    inside_col = "True" if "True" in ct.columns else (ct.columns[1] if len(ct.columns) > 1 else None)

    ct = ct.rename(columns={outside_col: "Outside zone", inside_col: "Inside zone"} if inside_col else {outside_col: "Outside zone"})
    if "Inside zone" not in ct.columns:
        ct["Inside zone"] = 0
    if "Outside zone" not in ct.columns:
        ct["Outside zone"] = 0

    ct["Total"] = ct["Outside zone"] + ct["Inside zone"]
    ct["% inside"] = (ct["Inside zone"] / ct["Total"] * 100).round(1)
    ct = ct[ct["Total"] > 0].sort_values("% inside", ascending=True)

    if len(ct) == 0:
        return

    fig, ax = plt.subplots(figsize=(10, max(4, len(ct) * 0.55)))
    fig.patch.set_facecolor("#f5f0e8")
    ax.set_facecolor("#f5f0e8")

    y = range(len(ct))
    width = 0.6
    ax.barh(y, ct["Outside zone"], width, color="#7f8c8d", alpha=0.8, label="Outside zone")
    ax.barh(y, ct["Inside zone"], width, left=ct["Outside zone"],
            color="#e74c3c", alpha=0.9, label="Inside zone")

    ax.set_yticks(list(y))
    ax.set_yticklabels(ct.index.tolist(), fontsize=10)
    ax.set_xlabel("Number of MRDS deposits", fontsize=11)
    ax.set_title(f"Commodity correlation with anomaly zones\n{site.name}", fontsize=13)
    ax.legend(loc="lower right", framealpha=0.9)

    for i, (_, row) in enumerate(ct.iterrows()):
        ax.text(
            row["Total"] + ct["Total"].max() * 0.01,
            i,
            f"{row['% inside']:.0f}% inside",
            va="center",
            fontsize=9,
            color="#333333",
        )

    ax.set_xlim(0, ct["Total"].max() * 1.25)
    plt.tight_layout()
    plt.savefig(
        paths.figures_dir / "04_commodity_correlation.png",
        dpi=150,
        bbox_inches="tight",
        facecolor="#f5f0e8",
    )
    plt.close(fig)


def write_provenance(
    paths: SitePaths,
    granule_id: str,
    extra: dict | None = None,
) -> None:
    paths.results_dir.mkdir(parents=True, exist_ok=True)
    prov: dict[str, Any] = {
        "site_id": paths.site.id,
        "granule_id": granule_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "packages": {},
    }
    try:
        import rasterio

        prov["packages"]["rasterio"] = rasterio.__version__
    except Exception:
        pass
    try:
        prov["packages"]["geopandas"] = gpd.__version__
    except Exception:
        pass
    try:
        prov["git_commit"] = (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=paths.repo_root,
                stderr=subprocess.DEVNULL,
            )
            .decode()
            .strip()
        )
    except Exception:
        prov["git_commit"] = None
    if extra:
        prov.update(extra)
    paths.site_provenance_json.write_text(json.dumps(prov, indent=2))


def download_aster(
    site: SiteConfig,
    paths: SitePaths,
    interactive_login: bool = True,
) -> str:
    """Search EarthData and download best granule; returns granule_id."""
    import earthaccess

    if interactive_login:
        earthaccess.login(strategy="netrc")

    bbox = search_bbox(site)
    results = earthaccess.search_data(
        short_name="AST_L1T",
        bounding_box=bbox,
        temporal=(site.temporal_start, site.temporal_end),
        count=10,
    )
    target = select_granule(results, site.bbox_wgs84, granule_id_override=site.granule_id)
    granule_id = extract_granule_id(target)
    paths.aster_dir.mkdir(parents=True, exist_ok=True)
    earthaccess.download(target, str(paths.aster_dir))
    return granule_id


def _feathered_mosaic(corrected_paths: list[Path], out_path: Path) -> None:
    """Merge granules with distance-weighted (feathered) blending at overlaps.

    For each granule:
    1. Reproject to a common grid (union extent at the reference resolution).
    2. Compute a per-pixel distance-to-nearest-nodata weight via EDT.
    3. Normalize each granule's weights to [0, 1] so mosaic scale is preserved.
    4. Blend: output = sum(w_i * v_i) / sum(w_i) at every pixel.

    All corrected granule files are expected to use 0 as the nodata sentinel.
    """
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.transform import from_bounds
    from rasterio.warp import reproject as warp_reproject
    from rasterio.warp import transform_bounds
    from scipy.ndimage import distance_transform_edt

    datasets = [rasterio.open(p) for p in corrected_paths]
    try:
        ref_ds = datasets[0]
        dst_crs = ref_ds.crs
        res_x = abs(ref_ds.transform.a)   # pixel width in CRS units
        res_y = abs(ref_ds.transform.e)   # pixel height in CRS units

        # Build union extent in the reference CRS.
        all_bounds = [
            transform_bounds(ds.crs, dst_crs, *ds.bounds) for ds in datasets
        ]
        left   = min(b[0] for b in all_bounds)
        bottom = min(b[1] for b in all_bounds)
        right  = max(b[2] for b in all_bounds)
        top    = max(b[3] for b in all_bounds)

        out_width  = max(1, int(round((right - left) / res_x)))
        out_height = max(1, int(round((top - bottom) / res_y)))
        out_transform = from_bounds(left, bottom, right, top, out_width, out_height)

        acc_val = np.zeros((out_height, out_width), dtype=np.float64)
        acc_w   = np.zeros((out_height, out_width), dtype=np.float64)

        for ds in datasets:
            src_arr = ds.read(1).astype(np.float32)
            # Build a float32 valid mask (1.0 = valid pixel).
            valid_src = ((src_arr > 0) & np.isfinite(src_arr)).astype(np.float32)
            # Zero out invalid pixels so bilinear resampling stays clean.
            src_arr[valid_src < 0.5] = 0.0

            # Reproject data to the common grid.
            dst_val = np.zeros((out_height, out_width), dtype=np.float32)
            warp_reproject(
                source=src_arr,
                destination=dst_val,
                src_transform=ds.transform,
                src_crs=ds.crs,
                dst_transform=out_transform,
                dst_crs=dst_crs,
                resampling=Resampling.bilinear,
                src_nodata=0.0,
                dst_nodata=0.0,
            )

            # Reproject valid mask (nearest-neighbor to keep it binary).
            dst_valid_f = np.zeros((out_height, out_width), dtype=np.float32)
            warp_reproject(
                source=valid_src,
                destination=dst_valid_f,
                src_transform=ds.transform,
                src_crs=ds.crs,
                dst_transform=out_transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest,
                src_nodata=0.0,
                dst_nodata=0.0,
            )
            dst_valid = dst_valid_f > 0.5

            # EDT: distance (pixels) from each valid pixel to nearest nodata.
            edt = distance_transform_edt(dst_valid)
            max_dist = float(edt.max())
            weight = (edt / max_dist) if max_dist > 0 else dst_valid.astype(np.float64)

            valid_and_finite = dst_valid & (np.abs(dst_val) > 0)
            acc_val[valid_and_finite] += weight[valid_and_finite] * dst_val[valid_and_finite]
            acc_w[valid_and_finite]   += weight[valid_and_finite]

        # Weighted average; pixels never covered by any granule stay 0 (nodata).
        out_arr = np.zeros((out_height, out_width), dtype=np.float32)
        nonzero_w = acc_w > 0
        out_arr[nonzero_w] = (acc_val[nonzero_w] / acc_w[nonzero_w]).astype(np.float32)

        profile = datasets[0].profile.copy()
        profile.update(
            height=out_height,
            width=out_width,
            transform=out_transform,
            crs=dst_crs,
            dtype=rasterio.float32,
            count=1,
            nodata=0.0,
        )
        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(out_arr, 1)
    finally:
        for ds in datasets:
            ds.close()


def _normalise_to_reference(
    src: np.ndarray,
    ref: np.ndarray,
    min_pixels: int = 100,
) -> np.ndarray:
    """Linearly rescale *src* so its mean/std match *ref* over valid pixels.

    "Valid" means finite and > 0 (ASTER uses 0 as nodata).  If either array
    has fewer than *min_pixels* valid samples the src is returned unchanged so
    the caller can fall back gracefully.  The correction is multiplicative-
    additive (mean + std match), which removes both additive path-radiance
    offsets and multiplicative gain differences between acquisitions.
    """
    src_valid = np.isfinite(src) & (src > 0)
    ref_valid = np.isfinite(ref) & (ref > 0)
    sv = src[src_valid]
    rv = ref[ref_valid]
    if sv.size < min_pixels or rv.size < min_pixels:
        return src
    s_mu, s_sd = float(sv.mean()), float(sv.std())
    r_mu, r_sd = float(rv.mean()), float(rv.std())
    if s_sd < 1e-6:
        return src
    out = src.astype(np.float64)
    out[src_valid] = (src[src_valid] - s_mu) / s_sd * r_sd + r_mu
    # Clip to [0, ∞) — negative radiance is non-physical.
    out = np.clip(out, 0, None)
    return out.astype(src.dtype)


def download_and_mosaic_aster(
    site: SiteConfig,
    paths: SitePaths,
    interactive_login: bool = True,
) -> str:
    """Download all bbox-covering ASTER granules, merge per-band, return mosaic granule_id.

    Always searches the full archive for every granule covering >5% of the site
    bbox (and has ≥3 TIR bands), regardless of whether ``site.granule_id`` is
    pinned.  Falls back to single-granule download only when no granule meets
    the coverage threshold.

    Before merging, each secondary granule is linearly normalised to the
    reference (first) granule's mean/std per band, eliminating the seam
    artefacts that arise from different acquisition dates, solar angles, and
    atmospheric conditions.
    """
    import shutil
    import tempfile

    import earthaccess
    import rasterio

    from critical_minerals_aster.spectral import score_granule

    if interactive_login:
        earthaccess.login(strategy="netrc")

    bbox = search_bbox(site)
    results = earthaccess.search_data(
        short_name="AST_L1T",
        bounding_box=bbox,
        temporal=(site.temporal_start, site.temporal_end),
        count=20,
    )

    covering = []
    for g in results:
        try:
            coverage, _, band_count = score_granule(g, site.bbox_wgs84)
            # Skip bundles larger than the per-site cap (default 20 MB keeps
            # TIR-only extracts; sites where only full VNIR+SWIR+TIR bundles
            # exist set a higher cap in their YAML).
            # 30% coverage threshold avoids edge-clipping granules that would
            # leave large nodata gaps in the merged output.
            if coverage > 0.30 and band_count >= 3 and g.size() < site.max_bundle_mb:
                covering.append(g)
        except Exception:
            continue

    if not covering:
        # No granule meets the threshold — fall back to best single granule.
        return download_aster(site, paths, interactive_login=False)

    paths.aster_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        for granule in covering:
            try:
                earthaccess.download(granule, tmpdir)
            except Exception as exc:
                print(f"  Warning: could not download granule: {exc}", file=sys.stderr)

        mosaic_id = f"{site.id}_mosaic"
        for band_num in (10, 11, 12, 13, 14):
            pattern = f"*_TIR_B{band_num}.tif"
            band_files = sorted(Path(tmpdir).glob(pattern))
            if not band_files:
                continue
            if len(band_files) == 1:
                shutil.copy(
                    band_files[0],
                    paths.aster_dir / f"{mosaic_id}_TIR_B{band_num}.tif",
                )
                continue

            # Read the reference (first) granule and normalise each secondary
            # granule to match its mean/std before merging.  Writing corrected
            # arrays to temporary GeoTIFFs lets _feathered_mosaic reproject and
            # distance-weight-blend them without modifying the originals.
            with rasterio.open(band_files[0]) as ds0:
                ref_arr = ds0.read(1).astype(np.float64)
                ref_profile = ds0.profile.copy()
            ref_arr[ref_arr == 0] = np.nan

            corrected_paths: list[Path] = [band_files[0]]
            for idx, bf in enumerate(band_files[1:], start=1):
                with rasterio.open(bf) as ds:
                    src_arr = ds.read(1).astype(np.float64)
                    src_profile = ds.profile.copy()
                src_arr[src_arr == 0] = np.nan
                normed = _normalise_to_reference(src_arr, ref_arr)
                # Write corrected array back to a temp GeoTIFF (restoring 0 nodata).
                tmp_path = Path(tmpdir) / f"normed_{idx}_B{band_num}.tif"
                write_arr = normed.copy()
                write_arr[~np.isfinite(write_arr)] = 0
                src_profile.update(dtype=rasterio.float32)
                with rasterio.open(tmp_path, "w", **src_profile) as dst:
                    dst.write(write_arr.astype(np.float32), 1)
                corrected_paths.append(tmp_path)

            out_path = paths.aster_dir / f"{mosaic_id}_TIR_B{band_num}.tif"
            _feathered_mosaic(corrected_paths, out_path)

    return mosaic_id


def auto_fetch_structure(
    site: SiteConfig,
    repo_root: Path,
    target_crs,
    *,
    buffer_m: float = 500.0,
    timeout: int = 30,
) -> "gpd.GeoDataFrame | None":
    """Auto-download fault data for a site that has no configured structure layers.

    Query order:
    1. USGS Quaternary Faults API (earthquake.usgs.gov) — fast, US-wide.
    2. USGS SGMC FeatureServer (ArcGIS REST) — all geological ages, all 48 states.

    Result is cached to ``data/structures/{site_id}_faults_auto.geojson`` so
    subsequent runs skip the network call.  Returns a GeoDataFrame in
    *target_crs* suitable for direct use as a structure layer, or ``None`` when
    both sources return no features.

    The function intentionally does **not** mutate the site YAML — callers use
    the returned GDF directly.  Run ``scripts/download_usgs_faults.py`` or
    ``scripts/download_sgmc_structures.py`` manually (or via ``run_site`` with
    ``--download-structures``) to persist the result permanently.
    """
    import json
    import urllib.parse
    import urllib.request

    out_dir = repo_root / "data" / "structures"
    out_path = out_dir / f"{site.id}_faults_auto.geojson"

    # ---- Serve from cache when available --------------------------------
    if out_path.is_file() and out_path.stat().st_size > 100:
        try:
            gdf = gpd.read_file(out_path)
            if not gdf.empty:
                if gdf.crs is None:
                    gdf = gdf.set_crs("EPSG:4326")
                return gdf.to_crs(target_crs)
        except Exception:
            pass  # fall through to re-fetch

    west, south, east, north = site.bbox_wgs84

    def _http_get(url: str) -> dict | None:
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "critical-minerals-aster/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            return None if "error" in data else data
        except Exception:
            return None

    # ---- 1. USGS Quaternary Faults API ----------------------------------
    _QFAULTS_URL = (
        "https://earthquake.usgs.gov/arcgis/rest/services/haz/Qfaults/MapServer/21/query"
    )
    params: dict = {
        "where": "1=1",
        "geometry": f"{west},{south},{east},{north}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": "4326",
        "outSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "FAULT_NAME,AGE,SLIP_RATE",
        "f": "geojson",
    }
    data = _http_get(_QFAULTS_URL + "?" + urllib.parse.urlencode(params))
    features = (data or {}).get("features", [])
    source_label = "USGS Quaternary Faults"

    # ---- 2. SGMC fallback -----------------------------------------------
    if not features:
        _SGMC_URL = (
            "https://services.arcgis.com/v01gqwM5QqNysAAi/arcgis/rest/services"
            "/SB_5888bf4fe4b05ccb964bab9d_USGS_SGMC_feature/FeatureServer/1/query"
        )
        _FAULT_RULE_IDS = (
            "11,12,13,21,22,23,24,29,30,31,33,34,35,36,"
            "42,43,44,45,46,47,48,49,50,51,52,53,54,62,63,64,65,66"
        )
        sgmc_params: dict = {
            "where": f"RuleID IN ({_FAULT_RULE_IDS})",
            "geometry": f"{west},{south},{east},{north}",
            "geometryType": "esriGeometryEnvelope",
            "inSR": "4326",
            "outSR": "4326",
            "spatialRel": "esriSpatialRelIntersects",
            "outFields": "STATE,DESCRIPTION,RuleID",
            "resultRecordCount": 2000,
            "f": "geojson",
        }
        sgmc_data = _http_get(_SGMC_URL + "?" + urllib.parse.urlencode(sgmc_params))
        features = (sgmc_data or {}).get("features", [])
        source_label = "USGS SGMC (all ages)"

    if not features:
        print(
            f"  [auto_fetch_structure] No fault features found for {site.id}; "
            "structure annotation skipped.",
            file=sys.stderr,
        )
        return None

    # Keep only line geometries.
    line_features = [
        f for f in features
        if (f.get("geometry") or {}).get("type", "") in ("LineString", "MultiLineString")
    ]
    if not line_features:
        return None

    # Cache to disk.
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps({"type": "FeatureCollection", "features": line_features}, indent=2)
    )
    print(
        f"  [auto_fetch_structure] {site.id}: {len(line_features)} features "
        f"from {source_label} → {out_path.relative_to(repo_root)}",
        file=sys.stderr,
    )

    gdf = gpd.read_file(out_path)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    return gdf.to_crs(target_crs)


def run_site(
    site: SiteConfig,
    repo_root: Path,
    *,
    download: bool = False,
    skip_figures: bool = False,
    global_limits: dict[str, tuple[float, float]] | None = None,
    skip_existing: bool = False,
) -> pd.DataFrame:
    """Run classification, write vectors/summary/provenance; optional EarthData download.

    Parameters
    ----------
    global_limits:
        Cross-site colorbar limits for Figure 01.  Pass the output of
        :func:`compute_global_limits` to make band-ratio colorbars
        comparable across sites.  When *None* per-site percentiles are used.
    """
    paths = site_paths_for(site, repo_root)

    if skip_existing:
        overlay_exists = (paths.figures_dir / "03_deposit_overlay.png").exists()
        prov_exists = paths.site_provenance_json.exists()
        if overlay_exists and prov_exists:
            print(
                f"Skipping {site.id} (outputs exist, use --force to regenerate)",
                file=sys.stderr,
            )
            return pd.DataFrame()

    if download:
        granule_id = download_aster(site, paths)
    else:
        granule_id = resolve_granule_id(site, paths)

    (
        zones,
        silica,
        carbonate,
        mafic,
        silica_cls,
        carbonate_cls,
        mafic_cls,
        combined,
        raster_bbox,  # WGS84 extent of the analysed (clipped) ASTER data
        raster_transform,
        raster_shape,
        raster_crs,
    ) = run_classification(site, paths, granule_id)
    paths.vectors_dir.mkdir(parents=True, exist_ok=True)
    zones.to_file(paths.strong_zones_geojson, driver="GeoJSON")

    # Compute hillshade for structural geology context in figs 01 and 02.
    # Gracefully degrades to None if DEM download/computation fails.
    hillshade: np.ndarray | None = None
    hs_transform = raster_transform  # fallback: use ASTER grid if terrain fails
    hs_shape = raster_shape
    if not skip_figures:
        try:
            from critical_minerals_aster.terrain import compute_hillshade_for_site

            _hs_result = compute_hillshade_for_site(
                site, paths, raster_transform, raster_shape, raster_crs
            )
            if _hs_result is not None:
                hillshade, hs_transform, hs_shape = _hs_result
        except Exception as exc:
            print(f"  [terrain] Hillshade skipped for {site.id}: {exc}", file=sys.stderr)

    if not skip_figures:
        save_composite_figure(site, paths, silica, carbonate, mafic)
        save_band_ratio_figure(
            site, paths, silica, carbonate, mafic,
            hillshade=hillshade, global_limits=global_limits,
        )
        save_classification_figure(
            site, paths, silica_cls, carbonate_cls, mafic_cls, combined, hillshade=hillshade
        )

    # Always compute MRDS join so figures 03+04 have deposit data.
    # Uses raster_bbox so only deposits within the TIR coverage area are shown.
    # Degrades gracefully if mrds.csv is missing.
    _deposits_gdf: gpd.GeoDataFrame | None = None
    try:
        from critical_minerals_aster.metrics import filter_mrds_bbox, read_mrds_national, simplify_commodity
        from critical_minerals_aster.mrds import mrds_to_points_gdf, spatial_join_deposits_zones

        mrds = read_mrds_national(paths)
        local = filter_mrds_bbox(mrds, raster_bbox)  # constrained to actual TIR coverage
        _deposits_gdf = mrds_to_points_gdf(local, zones.crs)
        joined, hits, _ = spatial_join_deposits_zones(_deposits_gdf, zones)
        hit_ids = joined[joined["index_right"].notna()].index.unique()
        _deposits_gdf["inside_zone"] = _deposits_gdf.index.isin(hit_ids)
        _deposits_gdf["commodity_group"] = _deposits_gdf["commod1"].apply(simplify_commodity)
    except FileNotFoundError:
        pass  # mrds.csv not downloaded yet — skip deposit figures

    # Load structure layers once — reused for annotation metrics AND figure 03.
    # When no layers are configured, attempt a lazy network fetch (USGS Quaternary
    # Faults → SGMC fallback) and cache the result to data/structures/.
    _structs_gdf: gpd.GeoDataFrame | None = None
    _target_crs = zones.crs if len(zones) else (_deposits_gdf.crs if _deposits_gdf is not None else None)
    if _deposits_gdf is not None and _target_crs is not None:
        if site.structure_layers:
            _structs_gdf = load_structure_layers(site, repo_root, _target_crs)
        else:
            _structs_gdf = auto_fetch_structure(site, repo_root, _target_crs)

    # Compute structure metrics before summary so they appear in the summary CSV.
    n_on_structure: int | None = None
    mean_nearest_m: float | None = None
    _annotated_gdf: gpd.GeoDataFrame | None = None
    provenance_extra: dict[str, Any] = {
        "n_zones": len(zones),
        "raster_bbox_wgs84": list(raster_bbox),
    }
    # Use whichever structure source is available — configured layers take
    # priority; auto-fetched GDF (from lazy SGMC fetch) is used as fallback.
    _buffer_m_for_annot = (
        site.structure_layers[0].buffer_m if site.structure_layers else 500.0
    )
    if _deposits_gdf is not None and _structs_gdf is not None and not _structs_gdf.empty:
        if site.structure_layers:
            _annotated_gdf = annotate_deposits_with_structure(
                _deposits_gdf, site, paths, structs=_structs_gdf
            )
        else:
            # Auto-fetched structures: annotate manually without a StructureLayer config.
            from critical_minerals_aster.structure import (
                nearest_structure_distance_m,
                points_on_structure,
            )
            _annotated_gdf = _deposits_gdf.copy()
            _annotated_gdf["nearest_structure_m"] = nearest_structure_distance_m(
                _annotated_gdf, _structs_gdf
            )
            _annotated_gdf["on_structure"] = points_on_structure(
                _annotated_gdf, _structs_gdf, _buffer_m_for_annot
            )
        n_on_structure = int(_annotated_gdf["on_structure"].sum())
        mean_nearest_m = float(_annotated_gdf["nearest_structure_m"].mean())
        provenance_extra["n_deposits_on_structure"] = n_on_structure
        provenance_extra["mean_nearest_structure_m"] = mean_nearest_m

    # Use raster_bbox (actual TIR coverage) instead of site.bbox_wgs84 so
    # MRDS deposits that lie outside the ASTER scene footprint are excluded.
    summary = compute_site_summary(
        site, paths, zones, granule_id, mrds_bbox=raster_bbox,
        n_on_structure=n_on_structure, mean_nearest_m=mean_nearest_m,
        annotated_deposits=_annotated_gdf,
    )

    if not skip_figures and _deposits_gdf is not None:
        save_deposit_overlay_figure(
            site, paths, zones, _deposits_gdf, repo_root,
            hillshade=hillshade,
            raster_transform=raster_transform,
            raster_shape=raster_shape,
            structs=_structs_gdf,
            n_on_structure=n_on_structure,
            n_total_deposits=len(_deposits_gdf),
            hs_transform=hs_transform,
            hs_shape=hs_shape,
        )
        save_commodity_correlation_figure(site, paths, _deposits_gdf)
        # Structure proximity strip chart — only when structure annotation exists.
        if _annotated_gdf is not None and "commodity_group" in _deposits_gdf.columns:
            # Merge commodity_group onto the annotated frame.
            _prox_df = _annotated_gdf.copy()
            _prox_df["commodity_group"] = _deposits_gdf["commodity_group"]
            _prox_df["inside_zone"] = _deposits_gdf.get("inside_zone", False)
            save_structure_proximity_figure(site, paths, _prox_df)

    write_site_summary(summary, paths.site_summary_csv)
    write_provenance(paths, granule_id, provenance_extra)
    return summary


def run_batch(
    site_ids: list[str],
    repo_root: Path,
    *,
    download: bool = False,
    skip_figures: bool = False,
    skip_existing: bool = False,
) -> list[pd.DataFrame]:
    from critical_minerals_aster.config import load_site_by_id

    sites_dir = repo_root / "sites"
    outputs: list[pd.DataFrame] = []
    for site_id in site_ids:
        site = load_site_by_id(site_id, sites_dir)
        try:
            outputs.append(
                run_site(
                    site,
                    repo_root,
                    download=download,
                    skip_figures=skip_figures,
                    skip_existing=skip_existing,
                )
            )
        except FileNotFoundError as exc:
            print(f"Skipping {site_id}: {exc}", file=sys.stderr)
    return outputs


def _run_site_worker(
    site_id: str,
    repo_root: Path,
    download: bool,
    skip_figures: bool,
    skip_existing: bool,
) -> pd.DataFrame:
    """Top-level function for ProcessPoolExecutor (must be picklable)."""
    from critical_minerals_aster.config import load_site_by_id

    site = load_site_by_id(site_id, repo_root / "sites")
    return run_site(
        site,
        repo_root,
        download=download,
        skip_figures=skip_figures,
        skip_existing=skip_existing,
    )


def run_batch_parallel(
    site_ids: list[str],
    repo_root: Path,
    *,
    workers: int = 2,
    download: bool = False,
    skip_figures: bool = False,
    skip_existing: bool = False,
) -> list[pd.DataFrame]:
    """Run sites in parallel using ProcessPoolExecutor."""
    import concurrent.futures

    outputs: list[pd.DataFrame] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                _run_site_worker,
                site_id,
                repo_root,
                download,
                skip_figures,
                skip_existing,
            ): site_id
            for site_id in site_ids
        }
        for future in concurrent.futures.as_completed(futures):
            site_id = futures[future]
            try:
                outputs.append(future.result())
            except Exception as exc:
                print(f"Error processing {site_id}: {exc}", file=sys.stderr)
    return outputs
