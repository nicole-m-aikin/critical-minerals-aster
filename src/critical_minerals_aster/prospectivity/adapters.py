"""Sensor adapters into a common ThermalScene / HyperspectralScene contract."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

from critical_minerals_aster.classification import combined_score
from critical_minerals_aster.config import SiteConfig, load_site_by_id
from critical_minerals_aster.paths import SitePaths, site_paths_for
from critical_minerals_aster.pipeline import resolve_granule_id
from critical_minerals_aster.spectral import load_classification_mosaic, load_ratio_mosaic


@dataclass
class ThermalScene:
    """Host-filter thermal features on a study grid."""

    sensor_id: str
    gsd_m: float
    qi: np.ndarray
    ci: np.ndarray
    mi: np.ndarray
    score: np.ndarray
    qa_mask: np.ndarray
    transform: Any
    crs: Any
    coverage_frac: float
    features_used: tuple[str, ...] = ("QI", "CI", "MI")


@dataclass
class HyperspectralScene:
    """Alteration-check SWIR features (Al-OH family), optionally warped."""

    sensor_id: str
    gsd_m: float
    aloh: np.ndarray
    kaolinite: np.ndarray
    alunite: np.ndarray
    qa_mask: np.ndarray
    transform: Any
    crs: Any
    coverage_frac: float
    wavelengths_nm: tuple[float, ...] = ()
    features_used_nm: tuple[float, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)


def _flatten_nm(used) -> tuple[float, ...]:
    """Normalize provenance wavelength metadata into a flat float tuple."""
    out: list[float] = []

    def _walk(obj) -> None:
        if obj is None:
            return
        if isinstance(obj, dict):
            for v in obj.values():
                _walk(v)
            return
        if isinstance(obj, (list, tuple)):
            for v in obj:
                _walk(v)
            return
        try:
            out.append(float(obj))
        except (TypeError, ValueError):
            pass

    _walk(used)
    return tuple(out)


def _percentile_rank(arr: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Within-scene equal-area percentile ranks in [0, 1]; NaN outside mask."""
    out = np.full(arr.shape, np.nan, dtype=np.float32)
    vals = arr[mask]
    if vals.size == 0:
        return out
    order = np.argsort(vals, kind="mergesort")
    ranks = np.empty(vals.size, dtype=np.float32)
    ranks[order] = (np.arange(vals.size, dtype=np.float32) + 0.5) / vals.size
    out[mask] = ranks
    return out


def load_aster_thermal_scene(
    site: SiteConfig,
    paths: SitePaths,
) -> ThermalScene:
    """Load AST_05 QI/CI/MI mosaic + combined score for a site."""
    granule_id = resolve_granule_id(site, paths)
    aster_dir = paths.aster_dir
    try:
        sil_cls, carb_cls, maf_cls, valid, transform, crs = load_classification_mosaic(
            aster_dir, granule_id
        )
        qi, ci, mi, _, _, _ = load_ratio_mosaic(aster_dir, granule_id)
        score = combined_score(sil_cls, carb_cls, maf_cls).astype(np.float32)
    except FileNotFoundError:
        qi, ci, mi, b10, transform, crs = load_ratio_mosaic(aster_dir, granule_id)
        from critical_minerals_aster.classification import classify_percentiles

        sil_cls, _, _ = classify_percentiles(qi)
        carb_cls, _, _ = classify_percentiles(ci)
        maf_cls, _, _ = classify_percentiles(mi)
        score = combined_score(sil_cls, carb_cls, maf_cls).astype(np.float32)
        valid = np.where(np.isfinite(b10), 1.0, np.nan)

    qa = np.isfinite(valid) & np.isfinite(qi)
    cov = float(qa.mean()) if qa.size else 0.0
    return ThermalScene(
        sensor_id="ASTER_AST05",
        gsd_m=90.0,
        qi=qi.astype(np.float32),
        ci=ci.astype(np.float32),
        mi=mi.astype(np.float32),
        score=score,
        qa_mask=qa,
        transform=transform,
        crs=crs,
        coverage_frac=cov,
    )


def _warp_to_grid(
    src_path: Path,
    *,
    shape: tuple[int, int],
    transform,
    crs,
) -> np.ndarray:
    with rasterio.open(src_path) as src:
        dst = np.full(shape, np.nan, dtype=np.float32)
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=crs,
            resampling=Resampling.bilinear,
            src_nodata=src.nodata,
            dst_nodata=np.nan,
        )
    return dst


def load_aviris_hyperspectral_scene(
    site_id: str,
    *,
    repo_root: Path,
    target_shape: tuple[int, int] | None = None,
    target_transform=None,
    target_crs=None,
) -> HyperspectralScene | None:
    """Load GEMx AVIRIS Al-OH / kaolinite / alunite TIFFs for a site.

    When a target grid is supplied, warps onto it. Returns None if TIFFs missing.
    """
    out_dir = repo_root / "data" / "gemx" / "aviris" / site_id
    paths = {k: out_dir / f"{site_id}_{k}.tif" for k in ("aloh", "kaolinite", "alunite")}
    if not all(p.is_file() for p in paths.values()):
        return None

    prov_path = repo_root / "results" / f"gemx_{site_id}_swir_provenance.json"
    sensor_id = "AVIRIS"
    gsd_m = 5.0
    features_used_nm: tuple[float, ...] = ()
    wavelengths: tuple[float, ...] = ()
    if prov_path.is_file():
        meta = json.loads(prov_path.read_text())
        collection = meta.get("collection", "")
        sensor_id = meta.get("sensor") or collection or sensor_id
        gsd_m = float(meta.get("mosaic_resolution_m") or gsd_m)
        used = meta.get("features_used_nm") or meta.get("wavelengths_nm") or []
        features_used_nm = _flatten_nm(used)

    if target_shape is not None and target_transform is not None and target_crs is not None:
        aloh = _warp_to_grid(paths["aloh"], shape=target_shape, transform=target_transform, crs=target_crs)
        kaol = _warp_to_grid(paths["kaolinite"], shape=target_shape, transform=target_transform, crs=target_crs)
        alun = _warp_to_grid(paths["alunite"], shape=target_shape, transform=target_transform, crs=target_crs)
        transform, crs = target_transform, target_crs
    else:
        with rasterio.open(paths["aloh"]) as src:
            aloh = src.read(1).astype(np.float32)
            transform, crs = src.transform, src.crs
            if src.nodata is not None:
                aloh[aloh == src.nodata] = np.nan
        with rasterio.open(paths["kaolinite"]) as src:
            kaol = src.read(1).astype(np.float32)
        with rasterio.open(paths["alunite"]) as src:
            alun = src.read(1).astype(np.float32)

    qa = np.isfinite(aloh)
    cov = float(qa.mean()) if qa.size else 0.0
    return HyperspectralScene(
        sensor_id=str(sensor_id),
        gsd_m=gsd_m,
        aloh=aloh,
        kaolinite=kaol,
        alunite=alun,
        qa_mask=qa,
        transform=transform,
        crs=crs,
        coverage_frac=cov,
        wavelengths_nm=wavelengths,
        features_used_nm=features_used_nm,
    )


def load_emit_hyperspectral_scene(
    site_id: str,
    *,
    repo_root: Path,
    target_shape: tuple[int, int] | None = None,
    target_transform=None,
    target_crs=None,
) -> HyperspectralScene | None:
    """Load EMIT pilot Al-OH products into the same contract as AVIRIS.

    Expects ``data/emit/{site_id}/{site_id}_{aloh,kaolinite,alunite}.tif``
    written by the EMIT adapter pilot.
    """
    out_dir = repo_root / "data" / "emit" / site_id
    paths = {k: out_dir / f"{site_id}_{k}.tif" for k in ("aloh", "kaolinite", "alunite")}
    if not all(p.is_file() for p in paths.values()):
        return None

    prov_path = repo_root / "results" / f"emit_{site_id}_swir_provenance.json"
    features_used_nm: tuple[float, ...] = (2165.0, 2200.0, 2210.0)
    if prov_path.is_file():
        meta = json.loads(prov_path.read_text())
        used = meta.get("features_used_nm") or []
        if used:
            features_used_nm = tuple(float(x) for x in used)

    if target_shape is not None and target_transform is not None and target_crs is not None:
        aloh = _warp_to_grid(paths["aloh"], shape=target_shape, transform=target_transform, crs=target_crs)
        kaol = _warp_to_grid(paths["kaolinite"], shape=target_shape, transform=target_transform, crs=target_crs)
        alun = _warp_to_grid(paths["alunite"], shape=target_shape, transform=target_transform, crs=target_crs)
        transform, crs = target_transform, target_crs
    else:
        with rasterio.open(paths["aloh"]) as src:
            aloh = src.read(1).astype(np.float32)
            transform, crs = src.transform, src.crs
        with rasterio.open(paths["kaolinite"]) as src:
            kaol = src.read(1).astype(np.float32)
        with rasterio.open(paths["alunite"]) as src:
            alun = src.read(1).astype(np.float32)

    qa = np.isfinite(aloh)
    return HyperspectralScene(
        sensor_id="EMIT_L2A",
        gsd_m=60.0,
        aloh=aloh,
        kaolinite=kaol,
        alunite=alun,
        qa_mask=qa,
        transform=transform,
        crs=crs,
        coverage_frac=float(qa.mean()) if qa.size else 0.0,
        features_used_nm=features_used_nm,
        metadata={"source": "emit_adapter"},
    )


def site_thermal_for_id(site_id: str, repo_root: Path) -> ThermalScene:
    site = load_site_by_id(site_id, repo_root / "sites")
    paths = site_paths_for(site, repo_root)
    return load_aster_thermal_scene(site, paths)


# Re-export rank helper for feature builder
percentile_rank = _percentile_rank
