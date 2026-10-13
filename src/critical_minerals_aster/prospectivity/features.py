"""Per-site feature/label store builder (parquet) under the common contract."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import xy
from scipy import ndimage

from critical_minerals_aster.config import list_site_ids, load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.prospectivity.adapters import (
    load_aster_thermal_scene,
    load_aviris_hyperspectral_scene,
    load_emit_hyperspectral_scene,
    percentile_rank,
)
from critical_minerals_aster.prospectivity.labels import (
    DEFAULT_POS_RADIUS_M,
    build_stratum_grid,
    label_grid,
    load_critical_deposits,
)
from critical_minerals_aster.prospectivity.metrics import FEATURE_COLS_ALL


def _dist_to_strong_m(score: np.ndarray, transform, min_score: int = 3) -> np.ndarray:
    """Euclidean distance (m) from each cell to nearest strong-anomaly cell."""
    strong = score >= min_score
    if not np.any(strong):
        return np.full(score.shape, np.nan, dtype=np.float32)
    # Pixel distance × |pixel size|
    px = float(abs(transform.a))
    py = float(abs(transform.e)) if transform.e != 0 else px
    dist_px = ndimage.distance_transform_edt(~strong)
    return (dist_px * ((px + py) / 2.0)).astype(np.float32)


def _equal_area_clay_mask(aloh: np.ndarray, overlap: np.ndarray, target_frac: float) -> np.ndarray:
    """Boolean clay mask matching ``target_frac`` of overlap cells (high Al-OH)."""
    out = np.zeros(aloh.shape, dtype=np.uint8)
    vals = aloh[overlap & np.isfinite(aloh)]
    if vals.size == 0 or target_frac <= 0:
        return out
    thr = float(np.quantile(vals, 1.0 - min(target_frac, 0.99)))
    out[(overlap & np.isfinite(aloh) & (aloh >= thr))] = 1
    return out


def build_site_feature_table(
    site_id: str,
    *,
    repo_root: Path,
    pos_radius_m: float = DEFAULT_POS_RADIUS_M,
    n_neg_per_pos: int = 5,
    random_state: int = 0,
    prefer_emit: bool = False,
) -> pd.DataFrame | None:
    """Build sampled cell feature/label rows for one site.

    Returns None when TIR mosaic is missing. Sets ``hsi_missing=1`` when no
    AVIRIS/EMIT SWIR product is available.
    """
    repo_root = Path(repo_root)
    site = load_site_by_id(site_id, repo_root / "sites")
    paths = site_paths_for(site, repo_root)

    try:
        thermal = load_aster_thermal_scene(site, paths)
    except FileNotFoundError:
        return None

    shape = thermal.qi.shape
    transform, crs = thermal.transform, thermal.crs
    qa = thermal.qa_mask

    hsi = None
    if prefer_emit:
        hsi = load_emit_hyperspectral_scene(
            site_id,
            repo_root=repo_root,
            target_shape=shape,
            target_transform=transform,
            target_crs=crs,
        )
    if hsi is None:
        hsi = load_aviris_hyperspectral_scene(
            site_id,
            repo_root=repo_root,
            target_shape=shape,
            target_transform=transform,
            target_crs=crs,
        )
    if hsi is None and not prefer_emit:
        hsi = load_emit_hyperspectral_scene(
            site_id,
            repo_root=repo_root,
            target_shape=shape,
            target_transform=transform,
            target_crs=crs,
        )

    hsi_missing = hsi is None
    if hsi is None:
        aloh = np.full(shape, np.nan, dtype=np.float32)
        kaol = np.full(shape, np.nan, dtype=np.float32)
        alun = np.full(shape, np.nan, dtype=np.float32)
        hsi_qa = np.zeros(shape, dtype=bool)
        hsi_sensor = ""
        hsi_gsd = np.nan
        features_used_nm = ()
    else:
        aloh, kaol, alun = hsi.aloh, hsi.kaolinite, hsi.alunite
        hsi_qa = hsi.qa_mask
        hsi_sensor = hsi.sensor_id
        hsi_gsd = hsi.gsd_m
        features_used_nm = hsi.features_used_nm

    # Strong-zone fraction for equal-area clay mask
    strong_frac = float((thermal.score >= 3)[qa].mean()) if qa.any() else 0.1
    clay = _equal_area_clay_mask(aloh, qa & hsi_qa, strong_frac)

    stratum = build_stratum_grid(
        site_id,
        shape=shape,
        transform=transform,
        crs=crs,
        bbox_wgs84=tuple(site.bbox_wgs84),
        repo_root=repo_root,
    )

    deposits = load_critical_deposits(site, paths, target_crs=crs)
    rng = np.random.default_rng(random_state)
    label, disc_yr, sample = label_grid(
        shape=shape,
        transform=transform,
        crs=crs,
        deposits=deposits,
        stratum_grid=stratum,
        pos_radius_m=pos_radius_m,
        n_neg_per_pos=n_neg_per_pos,
        rng=rng,
    )

    # Always keep valid TIR cells that are in the sample mask
    keep = sample & qa
    if not np.any(keep):
        # Fallback: random valid cells so site still appears in store
        valid_idx = np.flatnonzero(qa.ravel())
        if valid_idx.size == 0:
            return None
        pick = rng.choice(valid_idx, size=min(200, valid_idx.size), replace=False)
        keep = np.zeros(shape, dtype=bool)
        keep.ravel()[pick] = True

    rows_idx, cols_idx = np.where(keep)
    xs, ys = xy(transform, rows_idx, cols_idx, offset="center")

    qi_rank = percentile_rank(thermal.qi, qa)
    ci_rank = percentile_rank(thermal.ci, qa)
    mi_rank = percentile_rank(thermal.mi, qa)
    aloh_rank = percentile_rank(aloh, qa & hsi_qa) if not hsi_missing else np.full(shape, np.nan)
    kaol_rank = percentile_rank(kaol, qa & hsi_qa) if not hsi_missing else np.full(shape, np.nan)
    alun_rank = percentile_rank(alun, qa & hsi_qa) if not hsi_missing else np.full(shape, np.nan)
    dist_strong = _dist_to_strong_m(thermal.score, transform)

    df = pd.DataFrame({
        "site_id": site_id,
        "row": rows_idx.astype(np.int32),
        "col": cols_idx.astype(np.int32),
        "x": np.asarray(xs, dtype=np.float64),
        "y": np.asarray(ys, dtype=np.float64),
        "qi": thermal.qi[keep],
        "ci": thermal.ci[keep],
        "mi": thermal.mi[keep],
        "qi_rank": qi_rank[keep],
        "ci_rank": ci_rank[keep],
        "mi_rank": mi_rank[keep],
        "score": thermal.score[keep],
        "dist_strong_m": dist_strong[keep],
        "aloh": aloh[keep],
        "kaolinite": kaol[keep],
        "alunite": alun[keep],
        "aloh_rank": aloh_rank[keep],
        "kaolinite_rank": kaol_rank[keep],
        "alunite_rank": alun_rank[keep],
        "clay_mask": clay[keep].astype(np.uint8),
        "stratum_code": stratum[keep].astype(np.uint8),
        "hsi_missing": np.uint8(1 if hsi_missing else 0),
        "label": label[keep].astype(np.int8),
        "disc_yr": disc_yr[keep],
        "thermal_sensor": thermal.sensor_id,
        "hsi_sensor": hsi_sensor,
        "hsi_gsd_m": hsi_gsd,
        "features_used_nm": json.dumps(list(features_used_nm)),
    })
    return df


def feature_store_dir(repo_root: Path) -> Path:
    return Path(repo_root) / "data" / "prospectivity" / "features"


def write_site_parquet(df: pd.DataFrame, repo_root: Path) -> Path:
    out_dir = feature_store_dir(repo_root)
    out_dir.mkdir(parents=True, exist_ok=True)
    site_id = str(df["site_id"].iloc[0])
    path = out_dir / f"{site_id}.parquet"
    df.to_parquet(path, index=False)
    return path


def load_feature_store(
    repo_root: Path,
    site_ids: list[str] | None = None,
) -> pd.DataFrame:
    """Concatenate per-site parquet tables."""
    out_dir = feature_store_dir(repo_root)
    paths = sorted(out_dir.glob("*.parquet"))
    if site_ids is not None:
        want = set(site_ids)
        paths = [p for p in paths if p.stem in want]
    if not paths:
        return pd.DataFrame()
    return pd.concat([pd.read_parquet(p) for p in paths], ignore_index=True)


def build_feature_store(
    repo_root: Path,
    *,
    site_ids: list[str] | None = None,
    feature_ready_csv: Path | None = None,
    max_sites: int | None = None,
) -> pd.DataFrame:
    """Build parquet feature store for feature-ready (or all listed) sites."""
    repo_root = Path(repo_root)
    if site_ids is None:
        if feature_ready_csv and Path(feature_ready_csv).is_file():
            fr = pd.read_csv(feature_ready_csv)
            if "feature_ready" in fr.columns:
                mask = fr["feature_ready"].astype(str).str.lower().isin(["1", "true", "yes"])
                site_ids = fr.loc[mask, "site_id"].tolist()
            elif "tir_ok" in fr.columns:
                site_ids = fr.loc[fr["tir_ok"].astype(bool), "site_id"].tolist()
            # Prefer HSI-complete sites first so the MVP train cohort is richer
            if "hsi_ok" in fr.columns and site_ids:
                order = fr.set_index("site_id").reindex(site_ids)["hsi_ok"].fillna(0)
                site_ids = order.sort_values(ascending=False).index.tolist()
        if not site_ids:
            site_ids = list_site_ids(repo_root / "sites")

    rows = []
    written = []
    for i, sid in enumerate(site_ids):
        if max_sites is not None and i >= max_sites:
            break
        print(f"[{i+1}/{len(site_ids)}] features {sid}", flush=True)
        df = build_site_feature_table(sid, repo_root=repo_root)
        if df is None or df.empty:
            print(f"  skip {sid}: no TIR or empty sample")
            continue
        path = write_site_parquet(df, repo_root)
        written.append(path)
        rows.append(df)
        print(f"  wrote {path.name} n={len(df)} pos={int(df['label'].sum())} hsi_missing={int(df['hsi_missing'].iloc[0])}")

    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def model_matrix(df: pd.DataFrame, *, cols: list[str] | None = None) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return X, y, column names with NaN filled for HSI-missing cells."""
    cols = list(cols or FEATURE_COLS_ALL)
    X = df[cols].copy()
    # Missing HSI → fill ranks/depths with 0; flag already present
    for c in cols:
        if c.startswith(("aloh", "kaol", "alun", "clay")):
            X[c] = X[c].fillna(0.0)
        else:
            X[c] = X[c].fillna(X[c].median() if X[c].notna().any() else 0.0)
    y = df["label"].to_numpy(dtype=int)
    return X.to_numpy(dtype=np.float32), y, cols
