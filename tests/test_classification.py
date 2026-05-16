import sys
from pathlib import Path

import numpy as np
import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.classification import (  # noqa: E402
    classify_percentiles,
    combined_score,
)
from critical_minerals_aster.spectral import alteration_ratios, band_ratio  # noqa: E402


def test_band_ratio_divide_by_zero():
    a = np.array([1.0, 1.0])
    b = np.array([0.0, 2.0])
    r = band_ratio(a, b)
    assert np.isnan(r[0])
    assert r[1] == pytest.approx(0.5)


def test_alteration_ratios_shape():
    n = 4
    b12 = np.ones((n, n))
    b13 = np.full((n, n), 2.0)
    b14 = np.full((n, n), 2.0)
    silica, carb, mafic = alteration_ratios(b12, b13, b14)
    assert silica.shape == (n, n)
    assert np.allclose(silica, 1.0)


def test_combined_score_range():
    s = np.zeros((2, 2), dtype=np.uint8)
    c = np.array([[1, 2], [0, 1]], dtype=np.uint8)
    m = np.array([[0, 1], [2, 0]], dtype=np.uint8)
    out = combined_score(s, c, m)
    assert out.max() <= 6
    assert out.min() >= 0


def test_classify_percentiles_nan():
    ratio = np.array([[np.nan, 1.0], [2.0, 3.0]])
    classes, lo, hi = classify_percentiles(ratio, low_pct=50, high_pct=90)
    assert lo <= hi
    assert classes[0, 0] == 0
