"""Moment-matching destriping for ASTER TIR detector striping."""

from __future__ import annotations

import numpy as np
import pytest

from critical_minerals_aster.spectral import (
    destripe_amplitude,
    destripe_moment_match,
    destripe_scan_axis,
)


def _striped(shape=(40, 50), n_detectors=10, axis=0, amp=0.05, base=0.95):
    """Synthetic image: global base plus a per-detector mean offset."""
    arr = np.full(shape, base, dtype=np.float64)
    n = shape[axis]
    for d in range(n_detectors):
        idx = np.arange(d, n, n_detectors)
        offset = amp * (d - (n_detectors - 1) / 2) / ((n_detectors - 1) / 2)
        if axis == 0:
            arr[idx, :] += offset
        else:
            arr[:, idx] += offset
    return arr


def test_destripe_equalizes_detector_means():
    arr = _striped(amp=0.08)
    out = destripe_moment_match(arr, n_detectors=10, axis=0)
    means = [
        np.nanmean(out[d::10, :])
        for d in range(10)
    ]
    assert max(means) - min(means) < 1e-6


def test_destripe_preserves_global_mean():
    arr = _striped(amp=0.06)
    out = destripe_moment_match(arr, n_detectors=10, axis=0)
    assert np.nanmean(out) == pytest.approx(np.nanmean(arr), abs=1e-9)


def test_destripe_reduces_amplitude():
    arr = _striped(amp=0.1)
    before = destripe_amplitude(arr, n_detectors=10, axis=0)
    after = destripe_amplitude(
        destripe_moment_match(arr, n_detectors=10, axis=0), n_detectors=10, axis=0
    )
    assert before > 0.05
    assert after < 1e-6


def test_destripe_preserves_nans():
    arr = _striped()
    arr[5:8, 10:15] = np.nan
    out = destripe_moment_match(arr, n_detectors=10, axis=0)
    assert np.isnan(out[5:8, 10:15]).all()
    assert np.isfinite(out[0, 0])


def test_destripe_noop_on_empty_or_constant():
    empty = np.full((20, 20), np.nan)
    assert np.isnan(destripe_moment_match(empty)).all()
    const = np.ones((20, 20))
    out = destripe_moment_match(const)
    np.testing.assert_allclose(out, const)


def test_destripe_scan_axis_picks_stronger_modulation():
    row_striped = _striped(axis=0, amp=0.1)
    col_striped = _striped(axis=1, amp=0.1)
    assert destripe_scan_axis(row_striped, n_detectors=10) == 0
    assert destripe_scan_axis(col_striped, n_detectors=10) == 1
