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
    band_ratio,
    clip_bands_to_bbox,
    extract_granule_id,
    load_ratio_mosaic,
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
    "gpd.GeoDataFrame | None",  # tir_footprint — valid-pixel polygon in raster CRS
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

    tir_footprint is the polygon of valid (non-NaN) pixels derived from B10.
    It differs from raster_bbox in that it captures the actual scene boundary
    (diagonal ASTER swath edge) rather than the bounding rectangle.  Use it
    for all deposit and structure clipping instead of the bbox.
    """
    from critical_minerals_aster.spectral import compute_valid_data_footprint

    _cls_silica_path = paths.aster_dir / f"{granule_id}_cls_silica.tif"
    _ratio_silica_path = paths.aster_dir / f"{granule_id}_ratio_silica.tif"

    if _cls_silica_path.exists():
        # Per-granule classification mosaic path: percentile thresholds were applied
        # per-granule before merging, so the classified arrays are loaded directly
        # and classify_percentiles() is skipped.
        from critical_minerals_aster.spectral import load_classification_mosaic
        silica_cls, carbonate_cls, mafic_cls, valid_mask, transform, crs = load_classification_mosaic(
            paths.aster_dir, granule_id
        )
        (valid_mask, silica_cls, carbonate_cls, mafic_cls), transform = clip_bands_to_bbox(
            [valid_mask, silica_cls, carbonate_cls, mafic_cls], transform, crs, site.bbox_wgs84
        )
        # Ratio mosaic is still loaded for visualization figures (composite, band ratios).
        if _ratio_silica_path.exists():
            silica, carbonate, mafic, _, ratio_transform, ratio_crs = load_ratio_mosaic(
                paths.aster_dir, granule_id
            )
            (silica, carbonate, mafic), _ = clip_bands_to_bbox(
                [silica, carbonate, mafic], ratio_transform, ratio_crs, site.bbox_wgs84
            )
        else:
            silica = np.full(silica_cls.shape, np.nan, dtype=float)
            carbonate = np.full(silica_cls.shape, np.nan, dtype=float)
            mafic = np.full(silica_cls.shape, np.nan, dtype=float)

        raster_bbox: BBox = raster_bbox_wgs84(transform, silica_cls.shape, crs)
        tir_footprint = compute_valid_data_footprint(valid_mask, transform, crs)

        cp = site.classification
        assert cp is not None
        combined = combined_score(silica_cls, carbonate_cls, mafic_cls)
        zones = vectorize_strong_zones(combined, transform, crs, min_score=cp.strong_score_min)
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
            silica_cls.shape,
            crs,
            tir_footprint,
        )

    if _ratio_silica_path.exists():
        # Ratio mosaic path: download_and_mosaic_aster writes per-ratio GeoTIFFs so that
        # normalization operates on the actual signal (ratios) rather than individual
        # bands whose independent correction factors would distort the ratios.
        silica, carbonate, mafic, b10, transform, crs = load_ratio_mosaic(
            paths.aster_dir, granule_id
        )
        (b10, silica, carbonate, mafic), transform = clip_bands_to_bbox(
            [b10, silica, carbonate, mafic], transform, crs, site.bbox_wgs84
        )
    else:
        b10, _, b12, b13, b14, _, transform, crs = load_tir_bands_10_14(
            paths.aster_dir, granule_id
        )
        # Clip to site bbox so percentile thresholds and zone counts are
        # site-specific rather than whole-scene artifacts.  Shared-granule sites
        # (e.g. goldfield/silver_peak on the same ASTER swath) would otherwise
        # produce identical zone polygons from the full 60-90 km scene.
        (b10, b12, b13, b14), transform = clip_bands_to_bbox(
            [b10, b12, b13, b14], transform, crs, site.bbox_wgs84
        )
        silica, carbonate, mafic = alteration_ratios(b12, b13, b14)

    raster_bbox = raster_bbox_wgs84(transform, silica.shape, crs)
    tir_footprint = compute_valid_data_footprint(b10, transform, crs)

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
        silica.shape,
        crs,
        tir_footprint,
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


FIG03_VERSION = 10


def load_tir_footprint_from_provenance(
    site_id: str,
    repo_root: Path,
    zones_crs: Any,
) -> gpd.GeoDataFrame | None:
    """Load valid-pixel footprint polygon from site provenance JSON."""
    prov_path = repo_root / "results" / f"{site_id}_provenance.json"
    if not prov_path.exists():
        return None
    try:
        prov = json.loads(prov_path.read_text())
        geojson = prov.get("tir_footprint_wgs84")
        if geojson is None:
            return None
        gdf = gpd.GeoDataFrame.from_features(geojson["features"], crs="EPSG:4326")
        return gdf.to_crs(zones_crs)
    except Exception:
        return None


def fig03_outputs_current(repo_root: Path, paths: SitePaths) -> bool:
    """True when dual-panel fig 03 and matching provenance version exist."""
    overlay = paths.figures_dir / "03_deposit_overlay.png"
    prov_path = paths.site_provenance_json
    if not overlay.exists() or not prov_path.exists():
        return False
    try:
        prov = json.loads(prov_path.read_text())
        return prov.get("fig03_version") == FIG03_VERSION
    except Exception:
        return False


def _fig03_view_limits(
    hs_transform: Any,
    hs_shape: tuple[int, int] | None,
    zones: gpd.GeoDataFrame,
    ax_fallback: Any,
) -> tuple[tuple[float, float], tuple[float, float]]:
    if hs_transform is not None and hs_shape is not None:
        _r = hs_transform
        _rc, _cc = hs_shape
        _rx0, _rx1 = _r.c, _r.c + _r.a * _cc
        _ry0, _ry1 = _r.f + _r.e * _rc, _r.f
        _mx = (_rx1 - _rx0) * 0.01
        _my = (_ry1 - _ry0) * 0.01
        return (_rx0 - _mx, _rx1 + _mx), (_ry0 - _my, _ry1 + _my)
    if len(zones) > 0:
        _zb = zones.total_bounds
        _zm = max(_zb[2] - _zb[0], _zb[3] - _zb[1]) * 0.02
        return (_zb[0] - _zm, _zb[2] + _zm), (_zb[1] - _zm, _zb[3] + _zm)
    return ax_fallback.get_xlim(), ax_fallback.get_ylim()


def _fig03_hs_extent(hs_transform: Any, hs_shape: tuple[int, int]) -> tuple[float, float, float, float]:
    rows, cols = hs_shape
    t = hs_transform
    return (t.c, t.c + t.a * cols, t.f + t.e * rows, t.f)


def _fig03_draw_zones(ax: Any, zones: gpd.GeoDataFrame) -> None:
    """Draw strong anomaly zones with fill and crisp outlines."""
    if len(zones) == 0:
        return
    small = zones[zones["area_km2"] < 10]
    large = zones[zones["area_km2"] >= 10]
    if len(small):
        small.plot(ax=ax, color="#922b21", alpha=0.65, linewidth=0, zorder=2)
        small.plot(ax=ax, facecolor="none", edgecolor="#c0392b", linewidth=0.5, alpha=0.80, zorder=2)
    if len(large):
        large.plot(ax=ax, color="#4a0000", alpha=0.90, linewidth=0, zorder=2)
        large.plot(ax=ax, facecolor="none", edgecolor="#4a0000", linewidth=1.2, alpha=0.95, zorder=2)


def _fig03_draw_structure(
    ax: Any,
    site: SiteConfig,
    repo_root: Path,
    zones: gpd.GeoDataFrame,
    deposits: gpd.GeoDataFrame,
    tir_footprint: gpd.GeoDataFrame | None,
    structs: gpd.GeoDataFrame | None,
) -> tuple[bool, float]:
    has_structure = False
    buffer_m = site.structure_layers[0].buffer_m if site.structure_layers else 500.0
    if structs is None and site.structure_layers:
        target_crs = zones.crs if len(zones) else deposits.crs
        structs = load_structure_layers(site, repo_root, target_crs)
    if structs is None or structs.empty:
        return has_structure, buffer_m
    union_geom = structure_buffer_union(structs, buffer_m)
    if union_geom is not None:
        if tir_footprint is not None and len(tir_footprint):
            try:
                fp_for_clip = tir_footprint
                if tir_footprint.crs != structs.crs:
                    fp_for_clip = tir_footprint.to_crs(structs.crs)
                union_geom = union_geom.intersection(fp_for_clip.union_all())
            except Exception:
                pass
        if union_geom is not None and not union_geom.is_empty:
            # Filled buffer corridor — already clipped to TIR footprint above.
            gpd.GeoSeries([union_geom], crs=structs.crs).plot(
                ax=ax, color="#e67e22", alpha=0.30, linewidth=0, zorder=1,
            )
            has_structure = True
    # Draw individual fault/structure lines clipped to the TIR footprint so they
    # do not extend into no-data areas outside the scene boundary.
    if has_structure and tir_footprint is not None and len(tir_footprint):
        try:
            fp_for_lines = tir_footprint
            if tir_footprint.crs != structs.crs:
                fp_for_lines = tir_footprint.to_crs(structs.crs)
            structs_clipped = gpd.clip(structs, fp_for_lines)
        except Exception:
            structs_clipped = structs
        if len(structs_clipped):
            structs_clipped.plot(
                ax=ax, color="#d35400", linewidth=1.4, alpha=0.85, zorder=2,
            )
    elif has_structure:
        structs.plot(ax=ax, color="#d35400", linewidth=1.4, alpha=0.85, zorder=2)
    return has_structure, buffer_m


def _fig03_prepare_deposits(
    deposits: gpd.GeoDataFrame,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, int, str, bool]:
    outside = deposits[~deposits["inside_zone"]] if "inside_zone" in deposits.columns else deposits
    inside = deposits[deposits["inside_zone"]] if "inside_zone" in deposits.columns else gpd.GeoDataFrame()
    n_outside = len(outside)
    capped_note = ""
    show_cap = n_outside > 500
    if show_cap:
        outside = outside.iloc[:500]
        capped_note = f" (showing 500 of {n_outside})"
    return outside, inside, n_outside, capped_note, show_cap


def _fig03_draw_deposits(
    ax: Any,
    outside: gpd.GeoDataFrame,
    inside: gpd.GeoDataFrame,
    n_outside: int,
) -> None:
    ms = 10 if n_outside > 100 else 24
    alpha = 0.80 if n_outside > 100 else 0.9
    if len(outside):
        # Two-pass black halo: dark outline ring makes steelblue readable on
        # both the bright satellite and the gray hillshade panels.
        outside.plot(ax=ax, color="black", markersize=ms + 9, alpha=0.85, zorder=3)
        outside.plot(ax=ax, color="steelblue", markersize=ms, alpha=alpha, zorder=3)
    if len(inside):
        inside.plot(
            ax=ax, color="gold", markersize=110, marker="*", zorder=4,
            edgecolors="black", linewidths=0.7,
        )


def _fig03_draw_tir_footprint(
    ax: Any,
    site: SiteConfig,
    zones: gpd.GeoDataFrame,
    deposits: gpd.GeoDataFrame,
    tir_footprint: gpd.GeoDataFrame | None,
) -> float | None:
    if tir_footprint is None or not len(tir_footprint):
        return None
    plot_crs = zones.crs if len(zones) else (deposits.crs if len(deposits) else None)
    fp_plot = tir_footprint
    if plot_crs is not None and fp_plot.crs != plot_crs:
        fp_plot = fp_plot.to_crs(plot_crs)
    fp_plot.plot(
        ax=ax, facecolor="none", edgecolor="#1a252f",
        linewidth=2.2, linestyle="--", alpha=0.85, zorder=5,
    )
    try:
        from shapely.geometry import box as _box

        fp_wgs84 = tir_footprint.to_crs("EPSG:4326").iloc[0].geometry
        bbox_geom = _box(*site.bbox_wgs84)
        if bbox_geom.area > 0:
            return min(100.0, fp_wgs84.area / bbox_geom.area * 100)
    except Exception:
        pass
    return None


def save_deposit_overlay_figure(
    site: SiteConfig,
    paths: SitePaths,
    zones: gpd.GeoDataFrame,
    deposits: gpd.GeoDataFrame,
    repo_root: Path,
    hillshade: np.ndarray | None = None,
    tir_footprint: "gpd.GeoDataFrame | None" = None,
    structs: "gpd.GeoDataFrame | None" = None,
    n_on_structure: "int | None" = None,
    n_total_deposits: "int | None" = None,
    hs_transform: "rasterio.Affine | None" = None,
    hs_shape: "tuple[int, int] | None" = None,
    basemap_rgb: np.ndarray | None = None,
    basemap_source: str | None = None,
    basemap_cached: bool | None = None,
) -> dict[str, Any]:
    """Figure 03 — satellite context (left) and strong anomaly zones (right)."""
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches
    import matplotlib.ticker

    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    fig, (ax_ctx, ax_anom) = plt.subplots(
        1, 2, figsize=(20, 10), sharex=True, sharey=True,
        gridspec_kw={"wspace": 0.04},
        layout="constrained",
    )

    outside, inside, n_outside, capped_note, show_cap = _fig03_prepare_deposits(deposits)
    has_structure, buffer_m = _fig03_draw_structure(
        ax_anom, site, repo_root, zones, deposits, tir_footprint, structs,
    )
    _fig03_draw_structure(
        ax_ctx, site, repo_root, zones, deposits, tir_footprint, structs,
    )

    # Right panel: filled zones + outlines for legibility; outlines make edges crisp
    # even where small zones scatter over the hillshade.
    _fig03_draw_zones(ax_anom, zones)
    _fig03_draw_deposits(ax_ctx, outside, inside, n_outside)
    _fig03_draw_deposits(ax_anom, outside, inside, n_outside)

    xlim, ylim = _fig03_view_limits(hs_transform, hs_shape, zones, ax_anom)

    # Darker background gives hillshade more contrast (visible in no-data swath gaps).
    for ax in (ax_ctx, ax_anom):
        ax.set_facecolor("#d0d0d0")

    hs_extent: tuple[float, float, float, float] | None = None
    if hs_transform is not None and hs_shape is not None:
        hs_extent = _fig03_hs_extent(hs_transform, hs_shape)

    if basemap_rgb is not None and hs_extent is not None:
        ax_ctx.imshow(
            basemap_rgb, extent=hs_extent, origin="upper", zorder=0, interpolation="bilinear",
        )
    elif hillshade is not None and hs_extent is not None:
        hs_cmap = plt.cm.gray.copy()
        hs_cmap.set_bad(alpha=0.0)
        ax_ctx.imshow(
            hillshade, cmap=hs_cmap, alpha=0.75, vmin=0, vmax=1,
            extent=hs_extent, origin="upper", zorder=0,
        )

    if hillshade is not None and hs_extent is not None:
        hs_cmap = plt.cm.gray.copy()
        hs_cmap.set_bad(alpha=0.0)
        # Higher alpha (0.55) on right panel so topographic texture reads through zone fill.
        ax_anom.imshow(
            hillshade, cmap=hs_cmap, alpha=0.55, vmin=0, vmax=1,
            extent=hs_extent, origin="upper", zorder=0,
        )

    tir_coverage_pct = _fig03_draw_tir_footprint(ax_ctx, site, zones, deposits, tir_footprint)
    _fig03_draw_tir_footprint(ax_anom, site, zones, deposits, tir_footprint)

    for ax in (ax_ctx, ax_anom):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    km_fmt = matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1000:.0f}")
    for ax in (ax_ctx, ax_anom):
        ax.xaxis.set_major_formatter(km_fmt)
        ax.yaxis.set_major_formatter(km_fmt)

    crs_label = zones.crs if len(zones) else (deposits.crs if len(deposits) else None)
    epsg = crs_label.to_epsg() if crs_label is not None else None
    epsg_suffix = f" (EPSG:{epsg})" if epsg else ""
    for ax in (ax_ctx, ax_anom):
        ax.set_xlabel(f"Easting (km{epsg_suffix})")
    ax_ctx.set_ylabel("Northing (km)")
    ax_anom.set_ylabel("Northing (km)")

    ctx_title = "Satellite imagery"
    if basemap_rgb is None:
        ctx_title += " (unavailable — hillshade)"
    ax_ctx.set_title(ctx_title, fontsize=11)
    ax_anom.set_title("Strong alteration zones", fontsize=11)

    # Matplotlib legend fills COLUMN-MAJOR (col 0 top→bottom, then col 1, …).
    # To achieve the desired visual grouping with ncol=3:
    #   Visual col 0 (left):   fault trace, fault buffer
    #   Visual col 1 (center): MRDS outside, MRDS inside
    #   Visual col 2 (right):  TIR boundary, zone <10, zone ≥10
    #
    # Items must be ordered so that each column's items are contiguous:
    #   indices 0,1,2 → col 0;  indices 3,4,5 → col 1;  indices 6,7,8 → col 2
    #
    # Blank spacers at the foot of cols 0 and 1 pad those columns to 3 rows so
    # that all three zone entries stay in col 2 and no item bleeds left.
    _blank = mlines.Line2D([], [], linewidth=0, markersize=0, label='')
    _has_tir = tir_footprint is not None and len(tir_footprint) > 0

    _fault_trace = mlines.Line2D(
        [], [], color="#d35400", linewidth=2.0, alpha=0.85, label="Fault / structure trace"
    )
    _fault_buffer = mpatches.Patch(
        facecolor="#e67e22", edgecolor="none", alpha=0.45,
        label=f"Fault buffer (±{buffer_m:.0f} m, within TIR)"
    )
    _mrds_outside = mlines.Line2D(
        [], [], marker="o", color="w", markerfacecolor="steelblue",
        markeredgecolor="black", markeredgewidth=1.0, markersize=9,
        label=f"MRDS deposit (outside zone, n={n_outside}{capped_note})"
    )
    _mrds_inside = mlines.Line2D(
        [], [], marker="*", color="w", markerfacecolor="gold",
        markeredgecolor="black", markersize=16,
        label=f"MRDS deposit (inside zone, n={len(inside)})"
    )
    _tir_entry = mpatches.Patch(
        facecolor="none", edgecolor="#1a252f", linewidth=2.2,
        linestyle="--", alpha=0.85, label="TIR data boundary"
    )
    _zone_small = mpatches.Patch(
        facecolor="#922b21", alpha=0.70, label="Strong anomaly zone (< 10 km²)"
    )
    _zone_large = mpatches.Patch(
        facecolor="#4a0000", alpha=0.90, label="Strong anomaly zone (≥ 10 km²)"
    )

    if has_structure and _has_tir:
        # 9 items → clean 3×3 grid
        # Col 0: fault trace, fault buffer, blank
        # Col 1: MRDS outside, MRDS inside, blank
        # Col 2: TIR boundary, zone <10, zone ≥10
        legend_elements: list[Any] = [
            _fault_trace, _fault_buffer, _blank,
            _mrds_outside, _mrds_inside, _blank,
            _tir_entry, _zone_small, _zone_large,
        ]
    elif has_structure:
        # No TIR entry → 6 items, 2 rows × 3 cols
        # Col 0: fault trace, fault buffer
        # Col 1: MRDS outside, MRDS inside
        # Col 2: zone <10, zone ≥10
        legend_elements = [
            _fault_trace, _fault_buffer,
            _mrds_outside, _mrds_inside,
            _zone_small, _zone_large,
        ]
    elif _has_tir:
        # No structure → 5 items (+ 1 blank), 2 rows × 3 cols
        # Col 0: MRDS outside, MRDS inside
        # Col 1: TIR boundary, zone <10
        # Col 2: blank, zone ≥10
        legend_elements = [
            _mrds_outside, _mrds_inside,
            _tir_entry, _zone_small,
            _blank, _zone_large,
        ]
    else:
        # No structure, no TIR → 4 items, 2 rows × 2 cols
        legend_elements = [
            _mrds_outside, _mrds_inside,
            _zone_small, _zone_large,
        ]

    _ncol = 3 if (has_structure or _has_tir) else 2
    # Legend sits just below the suptitle inside the figure bbox.
    fig.legend(
        handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, 0.995),
        ncol=_ncol, framealpha=0.95, fontsize=8.5,
    )

    # Context annotations go on the LEFT panel so viewers reading A→B get
    # geographic/coverage info before interpreting the anomaly zones.
    if tir_coverage_pct is not None and tir_coverage_pct < 99.0:
        ax_ctx.text(
            0.02, 0.14,
            f"TIR data covers {tir_coverage_pct:.0f}% of site bbox",
            transform=ax_ctx.transAxes, fontsize=8, va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="none"),
            zorder=10,
        )

    if show_cap:
        ax_anom.text(
            0.98, 0.03,
            f"⚠ {n_outside:,} deposits total\n(displaying 500)",
            transform=ax_anom.transAxes, fontsize=7.5, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff3cd", alpha=0.92, edgecolor="#d4a017"),
            zorder=10,
        )

    if has_structure and n_on_structure is not None and n_total_deposits is not None:
        pct = (n_on_structure / n_total_deposits * 100) if n_total_deposits > 0 else 0.0
        ax_ctx.text(
            0.02, 0.09,
            f"{n_on_structure}/{n_total_deposits} deposits ({pct:.0f}%) within {buffer_m:.0f} m of structure",
            transform=ax_ctx.transAxes, fontsize=8.5, va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="none"),
            zorder=10,
        )

    fig.suptitle(
        f"Strong alteration zones vs MRDS deposits — {site.name}\n{_bbox_annotation(site)}",
        fontsize=12, y=1.01,
    )

    # N arrow on left panel; scale bar on both panels.
    ax_ctx.annotate(
        "N\n↑", xy=(0.04, 0.96), xycoords="axes fraction",
        ha="center", va="top", fontsize=13, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="gray"),
    )

    def _add_scalebar(ax: Any) -> None:
        try:
            from matplotlib_scalebar.scalebar import ScaleBar
            ax.add_artist(ScaleBar(1, "m", length_fraction=0.15, location="lower left",
                                   box_alpha=0.7, font_properties={"size": 9}))
        except Exception:
            xlim_bar = ax.get_xlim()
            ylim_bar = ax.get_ylim()
            bar_len_m = 10_000
            bar_x0 = xlim_bar[0] + (xlim_bar[1] - xlim_bar[0]) * 0.05
            bar_y = ylim_bar[0] + (ylim_bar[1] - ylim_bar[0]) * 0.04
            ax.plot([bar_x0, bar_x0 + bar_len_m], [bar_y, bar_y],
                    color="black", linewidth=3, solid_capstyle="butt", zorder=10)
            ax.text(bar_x0 + bar_len_m / 2, bar_y + (ylim_bar[1] - ylim_bar[0]) * 0.015,
                    "10 km", ha="center", va="bottom", fontsize=8, zorder=10)

    _add_scalebar(ax_ctx)
    _add_scalebar(ax_anom)

    plt.savefig(paths.figures_dir / "03_deposit_overlay.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    return {
        "fig03_version": FIG03_VERSION,
        "fig03_basemap_source": basemap_source,
        "fig03_basemap_cached": basemap_cached,
    }


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


def _max_mosaic(paths_list: list[Path], out_path: Path) -> None:
    """Merge uint8 classified maps by taking the elementwise max on a common grid.

    Reprojects each input to the union extent (CRS and resolution of the first
    file) using nearest-neighbour resampling, then writes max(cls_i) per pixel.
    Pixels covered by no granule stay 0 (nodata sentinel for classified maps).
    """
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.transform import from_bounds
    from rasterio.warp import reproject as warp_reproject, transform_bounds

    datasets = [rasterio.open(p) for p in paths_list]
    try:
        ref_ds = datasets[0]
        dst_crs = ref_ds.crs
        res_x = abs(ref_ds.transform.a)
        res_y = abs(ref_ds.transform.e)

        all_bounds = [transform_bounds(ds.crs, dst_crs, *ds.bounds) for ds in datasets]
        left   = min(b[0] for b in all_bounds)
        bottom = min(b[1] for b in all_bounds)
        right  = max(b[2] for b in all_bounds)
        top    = max(b[3] for b in all_bounds)

        out_width  = max(1, int(round((right - left) / res_x)))
        out_height = max(1, int(round((top - bottom) / res_y)))
        out_transform = from_bounds(left, bottom, right, top, out_width, out_height)

        acc = np.zeros((out_height, out_width), dtype=np.uint8)

        for ds in datasets:
            src_arr = ds.read(1).astype(np.uint8)
            dst_arr = np.zeros((out_height, out_width), dtype=np.uint8)
            warp_reproject(
                source=src_arr,
                destination=dst_arr,
                src_transform=ds.transform,
                src_crs=ds.crs,
                dst_transform=out_transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest,
                src_nodata=0,
                dst_nodata=0,
            )
            acc = np.maximum(acc, dst_arr)

        profile = ref_ds.profile.copy()
        profile.update(
            height=out_height, width=out_width, transform=out_transform,
            crs=dst_crs, dtype=rasterio.uint8, count=1, nodata=0,
        )
        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(acc, 1)
    finally:
        for ds in datasets:
            ds.close()


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


def _build_ratio_mosaic_from_paths(
    b10_files: list[Path],
    paths: SitePaths,
    mosaic_id: str,
) -> None:
    """Compute ratio mosaics from a list of local B10 TIF paths.

    Normalizes per-ratio distributions across granules before feather-blending,
    so radiometric differences between acquisition dates don't distort the ratios.
    Writes ``{mosaic_id}_ratio_{silica,carbonate,mafic}.tif`` and
    ``{mosaic_id}_TIR_B10.tif`` into ``paths.aster_dir``.
    """
    import tempfile

    import rasterio

    def _write_tmp(arr: np.ndarray, profile: dict, path: Path) -> None:
        p = profile.copy()
        p.update(dtype=rasterio.float32)
        out = arr.astype(np.float64).copy()
        out[~np.isfinite(out)] = 0.0
        with rasterio.open(path, "w", **p) as dst:
            dst.write(out.astype(np.float32), 1)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        granule_data: list[dict] = []
        for b10_file in b10_files:
            prefix = b10_file.stem.replace("_TIR_B10", "")
            try:
                with rasterio.open(b10_file) as ds:
                    b10 = ds.read(1).astype(np.float64)
                    profile = ds.profile.copy()
                b10[b10 == 0] = np.nan
                bands: dict[int, np.ndarray] = {}
                for bnum in (12, 13, 14):
                    bf = b10_file.parent / f"{prefix}_TIR_B{bnum}.tif"
                    with rasterio.open(bf) as ds:
                        arr = ds.read(1).astype(np.float64)
                    arr[arr == 0] = np.nan
                    bands[bnum] = arr
                granule_data.append({
                    "b10": b10,
                    "silica":    band_ratio(bands[13], bands[14]),
                    "carbonate": band_ratio(bands[13], bands[12]),
                    "mafic":     band_ratio(bands[12], bands[13]),
                    "profile": profile,
                })
            except Exception as exc:
                print(
                    f"  Warning: could not load granule {b10_file.name}: {exc}",
                    file=sys.stderr,
                )

        if not granule_data:
            return

        ref = granule_data[0]

        # --- Ratio mosaics (continuous values, used for visualization figures) ---
        if len(granule_data) == 1:
            # Single granule: write directly without reprojection to preserve the
            # original pixel grid.  Bilinear resampling in _feathered_mosaic would
            # shift edge-pixel values via interpolation with nodata neighbors,
            # distorting percentile thresholds and moving strong-anomaly zone boundaries.
            for ratio_name in ("silica", "carbonate", "mafic"):
                _write_tmp(
                    ref[ratio_name],
                    ref["profile"],
                    paths.aster_dir / f"{mosaic_id}_ratio_{ratio_name}.tif",
                )
            _write_tmp(
                ref["b10"],
                ref["profile"],
                paths.aster_dir / f"{mosaic_id}_TIR_B10.tif",
            )
        else:
            for ratio_name in ("silica", "carbonate", "mafic"):
                ref_arr = ref[ratio_name]
                ref_tmp = tmpdir / f"ref_{ratio_name}.tif"
                _write_tmp(ref_arr, ref["profile"], ref_tmp)
                corrected_paths: list[Path] = [ref_tmp]

                for idx, g in enumerate(granule_data[1:], 1):
                    normed = _normalise_to_reference(g[ratio_name], ref_arr)
                    sec_tmp = tmpdir / f"sec_{idx}_{ratio_name}.tif"
                    _write_tmp(normed, g["profile"], sec_tmp)
                    corrected_paths.append(sec_tmp)

                _feathered_mosaic(
                    corrected_paths,
                    paths.aster_dir / f"{mosaic_id}_ratio_{ratio_name}.tif",
                )

            b10_tmp_paths: list[Path] = []
            for idx, g in enumerate(granule_data):
                b10_tmp = tmpdir / f"b10_{idx}.tif"
                _write_tmp(g["b10"], g["profile"], b10_tmp)
                b10_tmp_paths.append(b10_tmp)
            _feathered_mosaic(
                b10_tmp_paths,
                paths.aster_dir / f"{mosaic_id}_TIR_B10.tif",
            )

        # --- Classification from the bbox-clipped ratio mosaic ---
        # Read the ratio mosaic written above, clip to the site bbox once, and
        # classify in a single pass.  Per-granule classify + max-merge inflates
        # the strong-anomaly fraction: with N independent granules each at the
        # 70th/90th percentile, P(max ≥ strong) = 1 − 0.9^N, so 3 granules →
        # ~27% strong and 6 granules → ~47% strong instead of the expected 10%.
        cp = paths.site.classification
        if cp is None:
            from critical_minerals_aster.config import ClassificationParams
            cp = ClassificationParams()

        with rasterio.open(paths.aster_dir / f"{mosaic_id}_TIR_B10.tif") as ds:
            mosaic_b10 = ds.read(1).astype(np.float64)
            mosaic_transform = ds.transform
            mosaic_crs = ds.crs
            mosaic_profile = ds.profile.copy()
        mosaic_b10[mosaic_b10 == 0] = np.nan

        mosaic_ratios: dict[str, np.ndarray] = {}
        for rname in ("silica", "carbonate", "mafic"):
            with rasterio.open(paths.aster_dir / f"{mosaic_id}_ratio_{rname}.tif") as ds:
                arr = ds.read(1).astype(np.float64)
            arr[arr == 0] = np.nan
            mosaic_ratios[rname] = arr

        (mosaic_b10_c, sil_c, carb_c, maf_c), cls_transform = clip_bands_to_bbox(
            [mosaic_b10, mosaic_ratios["silica"], mosaic_ratios["carbonate"], mosaic_ratios["mafic"]],
            mosaic_transform, mosaic_crs, paths.site.bbox_wgs84,
        )

        silica_cls, _, _ = classify_percentiles(sil_c, cp.low_pct, cp.high_pct)
        carbonate_cls, _, _ = classify_percentiles(carb_c, cp.low_pct, cp.high_pct)
        mafic_cls, _, _ = classify_percentiles(maf_c, cp.low_pct, cp.high_pct)
        valid_mask = np.isfinite(mosaic_b10_c).astype(np.uint8)

        rows, cols = silica_cls.shape
        cls_profile = mosaic_profile.copy()
        cls_profile.update(
            dtype="uint8", count=1, height=rows, width=cols,
            transform=cls_transform, crs=mosaic_crs, nodata=0,
        )
        for name, arr in [
            ("silica", silica_cls),
            ("carbonate", carbonate_cls),
            ("mafic", mafic_cls),
        ]:
            with rasterio.open(paths.aster_dir / f"{mosaic_id}_cls_{name}.tif", "w", **cls_profile) as dst:
                dst.write(arr, 1)
        with rasterio.open(paths.aster_dir / f"{mosaic_id}_cls_valid.tif", "w", **cls_profile) as dst:
            dst.write(valid_mask, 1)


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

    When local granule TIF files already exist in ``paths.aster_dir`` (non-mosaic
    B10 files), earthaccess download is skipped entirely — the local files are
    used directly to build the ratio mosaics.
    """
    import tempfile

    mosaic_id = f"{site.id}_mosaic"

    # Short-circuit: if both ratio mosaics (for figures) and classification mosaics
    # (for science) are already built, nothing to do.
    ratio_files = [paths.aster_dir / f"{mosaic_id}_ratio_{r}.tif" for r in ("silica", "carbonate", "mafic")]
    cls_files = [paths.aster_dir / f"{mosaic_id}_cls_{r}.tif" for r in ("silica", "carbonate", "mafic")]
    valid_file = paths.aster_dir / f"{mosaic_id}_cls_valid.tif"
    if all(p.is_file() for p in ratio_files) and all(p.is_file() for p in cls_files) and valid_file.is_file():
        return mosaic_id

    # Prefer local granule TIFs over re-downloading.  A granule is a set of
    # per-band files sharing the same prefix (everything before "_TIR_B10.tif").
    local_b10s = [
        p for p in sorted(paths.aster_dir.glob("*_TIR_B10.tif"))
        if "mosaic" not in p.name
    ]

    if local_b10s:
        # Build mosaic from local files — no earthaccess login needed.
        _build_ratio_mosaic_from_paths(local_b10s, paths, mosaic_id)
        return mosaic_id

    # No local files: fall back to earthaccess download.
    import earthaccess

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

        downloaded_b10s = sorted(Path(tmpdir).glob("*_TIR_B10.tif"))
        if not downloaded_b10s:
            return download_aster(site, paths, interactive_login=False)

        _build_ratio_mosaic_from_paths(downloaded_b10s, paths, mosaic_id)

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

    if skip_existing and fig03_outputs_current(repo_root, paths):
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
        raster_bbox,  # WGS84 bounding box of the analysed (clipped) ASTER data
        raster_transform,
        raster_shape,
        raster_crs,
        tir_footprint,  # valid-pixel polygon (diagonal scene shape, not bbox rectangle)
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
        local = filter_mrds_bbox(mrds, raster_bbox)  # fast bbox pre-filter
        _deposits_gdf = mrds_to_points_gdf(local, zones.crs)
        # Clip to actual valid-pixel footprint (handles diagonal scene edges).
        if tir_footprint is not None and len(tir_footprint):
            try:
                _deposits_gdf = gpd.clip(_deposits_gdf, tir_footprint)
            except Exception:
                pass  # fall back to bbox-only clipping on geometry errors
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
    if tir_footprint is not None:
        import json as _json
        try:
            provenance_extra["tir_footprint_wgs84"] = _json.loads(
                tir_footprint.to_crs("EPSG:4326").to_json()
            )
        except Exception:
            pass
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

    # Use tir_footprint (actual valid-pixel polygon) so deposits in diagonal
    # no-data corners of the ASTER scene are excluded.  mrds_bbox is the
    # rectangular pre-filter; tir_footprint is the precise polygon clip.
    summary = compute_site_summary(
        site, paths, zones, granule_id, mrds_bbox=raster_bbox,
        tir_footprint=tir_footprint,
        n_on_structure=n_on_structure, mean_nearest_m=mean_nearest_m,
        annotated_deposits=_annotated_gdf,
    )

    if not skip_figures and _deposits_gdf is not None:
        basemap_rgb: np.ndarray | None = None
        basemap_source: str | None = None
        basemap_cached: bool | None = None
        if hs_transform is not None and hs_shape is not None:
            try:
                from critical_minerals_aster.basemap import fetch_satellite_basemap_for_site

                _bm = fetch_satellite_basemap_for_site(
                    site, paths, hs_transform, hs_shape, raster_crs,
                )
                if _bm is not None:
                    basemap_rgb, basemap_source, basemap_cached = _bm
            except Exception as exc:
                print(f"  [basemap] Skipped for {site.id}: {exc}", file=sys.stderr)

        fig03_meta = save_deposit_overlay_figure(
            site, paths, zones, _deposits_gdf, repo_root,
            hillshade=hillshade,
            tir_footprint=tir_footprint,
            structs=_structs_gdf,
            n_on_structure=n_on_structure,
            n_total_deposits=len(_deposits_gdf),
            hs_transform=hs_transform,
            hs_shape=hs_shape,
            basemap_rgb=basemap_rgb,
            basemap_source=basemap_source,
            basemap_cached=basemap_cached,
        )
        provenance_extra.update(fig03_meta)
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
