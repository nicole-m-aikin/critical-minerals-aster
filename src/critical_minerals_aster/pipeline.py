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
from critical_minerals_aster.config import SiteConfig, search_bbox
from critical_minerals_aster.metrics import compute_site_summary, write_site_summary
from critical_minerals_aster.paths import SitePaths, site_paths_for
from critical_minerals_aster.spectral import (
    alteration_ratios,
    extract_granule_id,
    load_tir_bands_10_14,
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
    for name in sorted(aster_dir.iterdir()):
        match = _GRANULE_ID_RE.search(name.name)
        if match and "TIR_B10" in name.name:
            return match.group(1)
    raise ValueError(f"No ASTER granule TIR files found under {aster_dir}")


def run_classification(
    site: SiteConfig, paths: SitePaths, granule_id: str
) -> tuple[
    gpd.GeoDataFrame,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Classify, vectorize, return zones and class maps."""
    _, _, b12, b13, b14, _, transform, crs = load_tir_bands_10_14(
        paths.aster_dir, granule_id
    )
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
    )


def save_band_ratio_figure(
    site: SiteConfig,
    paths: SitePaths,
    silica: np.ndarray,
    carbonate: np.ndarray,
    mafic: np.ndarray,
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
) -> None:
    paths.figures_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    arrays = [silica_cls, carbonate_cls, mafic_cls, combined]
    titles = ["Silica classes", "Carbonate classes", "Mafic classes", "Combined score"]
    for ax, arr, title in zip(axes, arrays, titles):
        im = ax.imshow(arr, cmap="YlOrRd" if title != "Combined score" else "RdYlGn")
        ax.set_title(title)
        ax.axis("off")
        plt.colorbar(im, ax=ax, shrink=0.8)
    plt.suptitle(f"Alteration classification — {site.name}", fontsize=13)
    plt.tight_layout()
    plt.savefig(paths.figures_dir / "02_classification.png", dpi=150, bbox_inches="tight")
    plt.close()


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
        earthaccess.login(strategy="interactive")

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
    ) = run_classification(site, paths, granule_id)
    paths.vectors_dir.mkdir(parents=True, exist_ok=True)
    zones.to_file(paths.strong_zones_geojson, driver="GeoJSON")

    if not skip_figures:
        save_composite_figure(site, paths, silica, carbonate, mafic)
        save_band_ratio_figure(site, paths, silica, carbonate, mafic)
        save_classification_figure(
            site, paths, silica_cls, carbonate_cls, mafic_cls, combined
        )

    summary = compute_site_summary(site, paths, zones, granule_id)

    provenance_extra: dict[str, Any] = {"n_zones": len(zones)}
    if site.structure_layers:
        from critical_minerals_aster.metrics import filter_mrds_bbox, read_mrds_national
        from critical_minerals_aster.mrds import mrds_to_points_gdf

        mrds = read_mrds_national(paths)
        local = filter_mrds_bbox(mrds, site.bbox_wgs84)
        deposits = mrds_to_points_gdf(local, zones.crs)
        annotated = annotate_deposits_with_structure(deposits, site, paths)
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
