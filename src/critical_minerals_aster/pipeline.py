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
from critical_minerals_aster.structure import annotate_deposits_with_structure

_GRANULE_ID_RE = re.compile(r"(AST_L1T_\d+_\d+)")


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


def save_band_ratio_figure(
    site: SiteConfig,
    paths: SitePaths,
    silica: np.ndarray,
    carbonate: np.ndarray,
    mafic: np.ndarray,
    hillshade: np.ndarray | None = None,
) -> None:
    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    ratios = [
        (silica, "Silica/quartz (B13/B14)", "RdYlGn_r"),
        (carbonate, "Carbonate/dolomite (B13/B12)", "YlOrBr"),
        (mafic, "Mafic minerals (B12/B13)", "PuBu"),
    ]
    for ax, (ratio, title, cmap) in zip(axes, ratios):
        vmin, vmax = _percentile_limits(ratio, 2, 98)
        im = ax.imshow(
            ratio,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        if hillshade is not None:
            ax.imshow(hillshade, cmap="gray", alpha=0.25, vmin=0, vmax=255)
        ax.set_title(title)
        ax.axis("off")
        plt.colorbar(im, ax=ax, shrink=0.8)
    plt.suptitle(f"ASTER TIR Band Ratios — {site.name}", fontsize=13)
    plt.tight_layout()
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
        im = ax.imshow(arr, cmap="YlOrRd" if title != "Combined score" else "RdYlGn")
        if hillshade is not None:
            ax.imshow(hillshade, cmap="gray", alpha=0.25, vmin=0, vmax=255)
        ax.set_title(title)
        ax.axis("off")
        plt.colorbar(im, ax=ax, shrink=0.8)
    plt.suptitle(f"Alteration classification — {site.name}", fontsize=13)
    plt.tight_layout()
    plt.savefig(paths.figures_dir / "02_classification.png", dpi=150, bbox_inches="tight")
    plt.close()


def save_deposit_overlay_figure(
    site: SiteConfig,
    paths: SitePaths,
    zones: gpd.GeoDataFrame,
    deposits: gpd.GeoDataFrame,
) -> None:
    """Figure 03 — spatial map of strong anomaly zones with MRDS deposit overlay."""
    import matplotlib.lines as mlines

    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 10))

    if len(zones) > 0:
        small = zones[zones["area_km2"] < 10]
        large = zones[zones["area_km2"] >= 10]
        if len(small):
            small.plot(ax=ax, color="#c0392b", alpha=0.35, linewidth=0)
        if len(large):
            large.plot(ax=ax, color="#4a0000", alpha=0.85, linewidth=0)

    outside = deposits[~deposits["inside_zone"]] if "inside_zone" in deposits.columns else deposits
    inside = deposits[deposits["inside_zone"]] if "inside_zone" in deposits.columns else gpd.GeoDataFrame()

    if len(outside):
        outside.plot(ax=ax, color="steelblue", markersize=20, alpha=0.8, zorder=3)
    if len(inside):
        inside.plot(ax=ax, color="gold", markersize=70, marker="*", zorder=4,
                    edgecolor="black", linewidth=0.5)

    # Add terrain basemap for structural geology context (semi-transparent).
    try:
        import contextily as cx

        crs_str = zones.crs.to_string() if len(zones) else deposits.crs.to_string()
        providers_to_try = [
            cx.providers.OpenTopoMap,
            cx.providers.Stadia.StamenTerrain,  # type: ignore[attr-defined]
            cx.providers.CartoDB.Positron,  # type: ignore[attr-defined]
        ]
        for provider in providers_to_try:
            try:
                cx.add_basemap(ax, crs=crs_str, source=provider, alpha=0.35, zoom=11)
                break
            except Exception:
                continue
    except Exception:
        pass  # basemap is fully optional

    legend_elements = [
        mlines.Line2D([], [], color="#c0392b", linewidth=6, alpha=0.5, label="Strong anomaly zone (< 10 km²)"),
        mlines.Line2D([], [], color="#4a0000", linewidth=6, alpha=0.9, label="Strong anomaly zone (≥ 10 km²)"),
        mlines.Line2D([], [], marker="o", color="w", markerfacecolor="steelblue",
                      markersize=9, label=f"MRDS deposit (outside zone, n={len(outside)})"),
        mlines.Line2D([], [], marker="*", color="w", markerfacecolor="gold",
                      markeredgecolor="black", markersize=13,
                      label=f"MRDS deposit (inside zone, n={len(inside)})"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", framealpha=0.95, fontsize=9)
    ax.set_title(f"Strong alteration zones vs MRDS deposits\n{site.name}", fontsize=13)
    ax.set_xlabel(f"Easting (m, {zones.crs})" if len(zones) else "Easting (m)")
    ax.set_ylabel("Northing (m)")
    plt.tight_layout()
    plt.savefig(paths.figures_dir / "03_deposit_overlay.png", dpi=150, bbox_inches="tight")
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


def download_and_mosaic_aster(
    site: SiteConfig,
    paths: SitePaths,
    interactive_login: bool = True,
) -> str:
    """Download all bbox-covering ASTER granules, merge per-band, return mosaic granule_id.

    If ``site.granule_id`` is pinned, falls back to :func:`download_aster` so
    single-granule overrides are honoured without change.  Otherwise, every
    granule that covers >5 % of the site bbox (and has ≥3 TIR bands) is
    downloaded to a temporary directory, then merged band-by-band with
    ``rasterio.merge`` and written as ``{site_id}_mosaic_TIR_B{n}.tif``.
    """
    import shutil
    import tempfile

    import earthaccess
    import rasterio
    from rasterio.merge import merge as rasterio_merge

    from critical_minerals_aster.spectral import score_granule

    if interactive_login:
        earthaccess.login(strategy="netrc")

    # Honour a pinned granule_id — single-granule path.
    if site.granule_id:
        return download_aster(site, paths, interactive_login=False)

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
            if coverage > 0.05 and band_count >= 3:
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
            band_files = list(Path(tmpdir).glob(pattern))
            if not band_files:
                continue
            if len(band_files) == 1:
                shutil.copy(
                    band_files[0],
                    paths.aster_dir / f"{mosaic_id}_TIR_B{band_num}.tif",
                )
            else:
                datasets = [rasterio.open(f) for f in band_files]
                try:
                    merged, merged_transform = rasterio_merge(datasets)
                    profile = datasets[0].profile.copy()
                    profile.update(
                        {
                            "height": merged.shape[1],
                            "width": merged.shape[2],
                            "transform": merged_transform,
                            "count": 1,
                        }
                    )
                    out_path = paths.aster_dir / f"{mosaic_id}_TIR_B{band_num}.tif"
                    with rasterio.open(out_path, "w", **profile) as dst:
                        dst.write(merged[0], 1)
                finally:
                    for ds in datasets:
                        ds.close()

    return mosaic_id


def run_site(
    site: SiteConfig,
    repo_root: Path,
    *,
    download: bool = False,
    skip_figures: bool = False,
) -> pd.DataFrame:
    """Run classification, write vectors/summary/provenance; optional EarthData download."""
    paths = site_paths_for(site, repo_root)
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
    if not skip_figures:
        try:
            from critical_minerals_aster.terrain import compute_hillshade_for_site

            hillshade = compute_hillshade_for_site(
                site, paths, raster_transform, raster_shape, raster_crs
            )
        except Exception as exc:
            print(f"  [terrain] Hillshade skipped for {site.id}: {exc}", file=sys.stderr)

    if not skip_figures:
        save_composite_figure(site, paths, silica, carbonate, mafic)
        save_band_ratio_figure(site, paths, silica, carbonate, mafic, hillshade=hillshade)
        save_classification_figure(
            site, paths, silica_cls, carbonate_cls, mafic_cls, combined, hillshade=hillshade
        )

    # Use raster_bbox (actual TIR coverage) instead of site.bbox_wgs84 so
    # MRDS deposits that lie outside the ASTER scene footprint are excluded.
    summary = compute_site_summary(site, paths, zones, granule_id, mrds_bbox=raster_bbox)

    provenance_extra: dict[str, Any] = {
        "n_zones": len(zones),
        "raster_bbox_wgs84": list(raster_bbox),
    }

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

    if not skip_figures and _deposits_gdf is not None:
        save_deposit_overlay_figure(site, paths, zones, _deposits_gdf)
        save_commodity_correlation_figure(site, paths, _deposits_gdf)

    if site.structure_layers and _deposits_gdf is not None:
        annotated = annotate_deposits_with_structure(_deposits_gdf, site, paths)
        provenance_extra["n_deposits_on_structure"] = int(annotated["on_structure"].sum())
        provenance_extra["mean_nearest_structure_m"] = float(
            annotated["nearest_structure_m"].mean()
        )

    write_site_summary(summary, paths.site_summary_csv)
    write_provenance(paths, granule_id, provenance_extra)
    return summary


def run_batch(
    site_ids: list[str],
    repo_root: Path,
    *,
    download: bool = False,
    skip_figures: bool = False,
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
                )
            )
        except FileNotFoundError as exc:
            print(f"Skipping {site_id}: {exc}", file=sys.stderr)
    return outputs
