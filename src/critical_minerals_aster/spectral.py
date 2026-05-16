"""ASTER TIR band I/O and alteration ratios."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import rasterio
from rasterio.profiles import Profile


def band_ratio(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(b != 0, a / b, np.nan)


def load_tir_band(
    aster_dir: os.PathLike[str] | str,
    granule_id: str,
    band_num: int,
) -> tuple[np.ndarray, Profile, rasterio.Affine, rasterio.crs.CRS]:
    """Load one ASTER L1T TIR band (float32), masking zero to NaN."""
    path = Path(aster_dir) / f"{granule_id}_TIR_B{band_num}.tif"
    with rasterio.open(path) as src:
        data = src.read(1).astype(float)
        profile = src.profile
        transform = src.transform
        crs = src.crs
    data[data == 0] = np.nan
    return data, profile, transform, crs


def load_tir_bands_10_14(
    aster_dir: os.PathLike[str] | str,
    granule_id: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    Profile,
    rasterio.Affine,
    rasterio.crs.CRS,
]:
    """Load TIR bands B10–B14; profile/transform/crs come from B10."""
    b10, profile, transform, crs = load_tir_band(aster_dir, granule_id, 10)
    b11, _, _, _ = load_tir_band(aster_dir, granule_id, 11)
    b12, _, _, _ = load_tir_band(aster_dir, granule_id, 12)
    b13, _, _, _ = load_tir_band(aster_dir, granule_id, 13)
    b14, _, _, _ = load_tir_band(aster_dir, granule_id, 14)
    return b10, b11, b12, b13, b14, profile, transform, crs


def alteration_ratios(
    b12: np.ndarray, b13: np.ndarray, b14: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Silica, carbonate, and mafic ratios (same definitions as notebooks)."""
    silica = band_ratio(b13, b14)
    carbonate = band_ratio(b13, b12)
    mafic = band_ratio(b12, b13)
    return silica, carbonate, mafic
