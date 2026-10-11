import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from gemx_goldfield_correspondence import (  # noqa: E402
    binary_lift,
    cliffs_delta,
    jaccard,
    median_enrichment,
)


def test_median_enrichment_independent_is_one():
    rng = np.random.default_rng(0)
    x = rng.random(200)
    assert median_enrichment(x, x) == 1.0


def test_median_enrichment_nested_is_above_one():
    inside = np.array([0.2, 0.2, 0.2])
    outside = np.array([0.1, 0.1, 0.1])
    assert median_enrichment(inside, outside) == 2.0


def test_binary_lift_independent_is_one():
    event = np.array([True, False, True, False])
    group = np.array([True, True, False, False])
    assert binary_lift(event, group) == 1.0


def test_binary_lift_nested_is_above_one():
    event = np.array([True, True, False, False])
    group = np.array([True, True, False, False])
    assert binary_lift(event, group) == 2.0


def test_jaccard_identical_is_one():
    a = np.array([True, True, False])
    assert jaccard(a, a) == 1.0


def test_jaccard_disjoint_is_zero():
    a = np.array([True, False])
    b = np.array([False, True])
    assert jaccard(a, b) == 0.0


def test_cliffs_delta_identical_near_zero():
    x = np.arange(20, dtype=float)
    assert abs(cliffs_delta(x, x)) < 1e-12
