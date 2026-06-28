"""Tests for uncertainty.py — threshold_ensemble and save_prospectivity_rasters."""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.uncertainty import (  # noqa: E402
    DEFAULT_SCENARIOS,
    threshold_ensemble,
    save_prospectivity_rasters,
)


def _synthetic_ratios(shape=(20, 20), seed=0):
    """Uniform random ratio arrays with some NaN nodata pixels."""
    rng = np.random.default_rng(seed)
    silica = rng.uniform(0.8, 1.2, shape)
    carbonate = rng.uniform(0.9, 1.3, shape)
    mafic = rng.uniform(0.7, 1.1, shape)
    # Add a nodata border
    for arr in (silica, carbonate, mafic):
        arr[0, :] = np.nan
        arr[-1, :] = np.nan
    return silica, carbonate, mafic


# ---------------------------------------------------------------------------
# threshold_ensemble
# ---------------------------------------------------------------------------


def test_output_shape_matches_input():
    silica, carbonate, mafic = _synthetic_ratios((30, 40))
    mean, std = threshold_ensemble(silica, carbonate, mafic)
    assert mean.shape == (30, 40)
    assert std.shape == (30, 40)


def test_mean_in_range():
    silica, carbonate, mafic = _synthetic_ratios()
    mean, _ = threshold_ensemble(silica, carbonate, mafic)
    valid = mean[np.isfinite(mean)]
    assert valid.min() >= 0.0
    assert valid.max() <= 1.0


def test_std_non_negative():
    silica, carbonate, mafic = _synthetic_ratios()
    _, std = threshold_ensemble(silica, carbonate, mafic)
    assert np.all(std[np.isfinite(std)] >= 0.0)


def test_nodata_propagates():
    """NaN in any input ratio → NaN in both output arrays at that pixel."""
    silica, carbonate, mafic = _synthetic_ratios()
    # Force a specific pixel to NaN in silica only
    silica[5, 5] = np.nan
    mean, std = threshold_ensemble(silica, carbonate, mafic)
    assert np.isnan(mean[5, 5])
    assert np.isnan(std[5, 5])


def test_mean_discrete_values_with_3_scenarios():
    """With 3 scenarios, valid mean values can only be 0, 1/3, 2/3, or 1."""
    silica, carbonate, mafic = _synthetic_ratios()
    mean, _ = threshold_ensemble(silica, carbonate, mafic)
    valid = mean[np.isfinite(mean)]
    allowed = {0.0, 1 / 3, 2 / 3, 1.0}
    # Use approximate comparison
    for v in valid:
        assert any(abs(v - a) < 1e-5 for a in allowed), f"Unexpected mean value: {v}"


def test_uniform_high_ratios_produce_high_mean():
    """Pixels with very high silica and carbonate ratios should be strong in all scenarios."""
    shape = (10, 10)
    # Extremely high values → percentile thresholds will all classify these as strong
    silica = np.full(shape, 10.0)
    carbonate = np.full(shape, 10.0)
    mafic = np.full(shape, 0.1)  # low mafic (mafic ratio = B12/B13 — low is non-mafic)
    mean, std = threshold_ensemble(silica, carbonate, mafic)
    # At least some pixels should be strongly anomalous (mean=1.0 for silica+carbonate)
    assert (mean == 1.0).any()


def test_std_zero_when_all_scenarios_agree():
    """If all pixels are unanimously strong or not-strong, std should be 0."""
    shape = (10, 10)
    # All zeros → never strong in any scenario → mean=0, std=0
    silica = np.full(shape, 0.5)
    carbonate = np.full(shape, 0.5)
    mafic = np.full(shape, 0.5)
    mean, std = threshold_ensemble(silica, carbonate, mafic)
    # Background pixels should all agree (std=0)
    assert np.all(std[mean == 0.0] == pytest.approx(0.0, abs=1e-5))


def test_custom_scenarios():
    """Accepts a single custom scenario without error."""
    silica, carbonate, mafic = _synthetic_ratios()
    mean, std = threshold_ensemble(
        silica, carbonate, mafic, scenarios=((70.0, 90.0),)
    )
    # Single scenario: mean is binary (0 or 1), std is always 0
    valid = mean[np.isfinite(mean)]
    assert set(np.unique(np.round(valid, 3))).issubset({0.0, 1.0})
    assert np.all(std[np.isfinite(std)] == pytest.approx(0.0, abs=1e-5))


# ---------------------------------------------------------------------------
# save_prospectivity_rasters
# ---------------------------------------------------------------------------


def test_save_creates_two_tifs():
    import rasterio
    from rasterio.transform import from_bounds

    silica, carbonate, mafic = _synthetic_ratios((20, 20))
    mean, std = threshold_ensemble(silica, carbonate, mafic)
    transform = from_bounds(-118, 41, -117, 42, 20, 20)
    crs = rasterio.crs.CRS.from_epsg(4326)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        mean_path, std_path = save_prospectivity_rasters(mean, std, transform, crs, out_dir)

        assert mean_path.exists()
        assert std_path.exists()

        with rasterio.open(mean_path) as src:
            arr = src.read(1)
            assert arr.shape == (20, 20)
            assert src.crs == crs
