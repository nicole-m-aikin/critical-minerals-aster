"""
Pixel-level prospectivity probability and uncertainty estimation.

The deterministic pipeline produces hard binary zones (strong / not-strong)
from fixed percentile thresholds (70th / 90th by default).  These thresholds
are informed choices, not physical constants — different but equally defensible
values change which pixels cross the "strong anomaly" boundary.

This module implements a threshold-perturbation ensemble:

    for each scenario in [(65, 85), (70, 90), (75, 95)]:
        classify silica, carbonate, mafic → combined score
        binary strong = combined_score >= strong_score_min

    P(strong | pixel) = fraction of scenarios where pixel is strong
    σ(strong | pixel) = std across scenario binary maps

P = 1.0 means the pixel is strong in every scenario (high confidence).
P = 0.0 means it is never strong (clearly background).
P = 0.33 or 0.67 means the pixel sits near the classification boundary
        and its label is threshold-sensitive (uncertain).

This is a model-uncertainty estimate — it captures sensitivity to the
main modeling assumption (the percentile thresholds) rather than
observational noise.  For observational uncertainty (acquisition-date
variability) multiple temporal granules would be required.

Outputs per site
----------------
    data/sites/{id}/rasters/prospectivity_mean.tif   — P(strong), float32 [0–1]
    data/sites/{id}/rasters/prospectivity_std.tif    — σ(strong), float32 [0–0.5]
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine

from critical_minerals_aster.classification import classify_percentiles, combined_score

# Default threshold scenarios: (low_pct, high_pct)
# Chosen to bracket the nominal (70, 90) setting symmetrically.
DEFAULT_SCENARIOS: tuple[tuple[float, float], ...] = (
    (65.0, 85.0),  # loose  — more pixels classified as strong/moderate
    (70.0, 90.0),  # nominal — matches pipeline default
    (75.0, 95.0),  # tight  — only the most extreme pixels classified strong
)


def threshold_ensemble(
    silica: np.ndarray,
    carbonate: np.ndarray,
    mafic: np.ndarray,
    scenarios: tuple[tuple[float, float], ...] = DEFAULT_SCENARIOS,
    strong_score_min: int = 3,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-pixel prospectivity probability and uncertainty via
    threshold-perturbation ensemble.

    Parameters
    ----------
    silica, carbonate, mafic :
        Ratio arrays (float, NaN = nodata) from load_ratio_mosaic().
    scenarios :
        Sequence of (low_pct, high_pct) threshold pairs.  Each defines one
        classify_percentiles() call.
    strong_score_min :
        Combined score threshold for "strong anomaly" label (default 3).

    Returns
    -------
    mean : np.ndarray (float32)
        Fraction of scenarios in which this pixel is classified as strong.
        Range [0, 1].  Pixels with all-NaN ratios remain NaN.
    std : np.ndarray (float32)
        Sample standard deviation across scenario binary maps.
        Maximum possible value ≈ 0.5 (at mean = 0.5, N = 3).

    Notes
    -----
    With N=3 scenarios the possible mean values are {0, 0.33, 0.67, 1.0}.
    Adding more scenarios would produce a finer-grained probability surface
    but the four levels are already interpretable:
        0.00 → never strong (background)
        0.33 → strong only in the loosest scenario (marginal)
        0.67 → strong in loose + nominal but not tight (moderate confidence)
        1.00 → strong in all scenarios (high confidence target)
    """
    binary_maps: list[np.ndarray] = []

    for low_pct, high_pct in scenarios:
        s_cls, _, _ = classify_percentiles(silica, low_pct, high_pct)
        c_cls, _, _ = classify_percentiles(carbonate, low_pct, high_pct)
        m_cls, _, _ = classify_percentiles(mafic, low_pct, high_pct)

        score = combined_score(s_cls, c_cls, m_cls)
        strong = (score >= strong_score_min).astype(np.float32)

        # Propagate NaN from any nodata pixel
        nodata_mask = np.isnan(silica) | np.isnan(carbonate) | np.isnan(mafic)
        strong[nodata_mask] = np.nan
        binary_maps.append(strong)

    stack = np.stack(binary_maps, axis=0)  # shape (N, H, W)
    mean = np.nanmean(stack, axis=0).astype(np.float32)
    std = np.nanstd(stack, axis=0, ddof=0).astype(np.float32)

    # Restore NaN where all inputs are nodata
    all_nodata = np.all(np.isnan(stack), axis=0)
    mean[all_nodata] = np.nan
    std[all_nodata] = np.nan

    return mean, std


def save_prospectivity_rasters(
    mean: np.ndarray,
    std: np.ndarray,
    transform: Affine,
    crs,
    out_dir: Path,
) -> tuple[Path, Path]:
    """Write prospectivity_mean.tif and prospectivity_std.tif to *out_dir*.

    Returns (mean_path, std_path).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "width": mean.shape[1],
        "height": mean.shape[0],
        "count": 1,
        "crs": crs,
        "transform": transform,
        "nodata": float("nan"),
        "compress": "lzw",
    }

    mean_path = out_dir / "prospectivity_mean.tif"
    std_path = out_dir / "prospectivity_std.tif"

    for path, arr in [(mean_path, mean), (std_path, std)]:
        with rasterio.open(path, "w", **profile) as dst:
            dst.write(arr, 1)

    return mean_path, std_path
