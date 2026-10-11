"""Unit tests for the n=2 cutoff-sensitivity helpers.  No cubes."""

import sys
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from gemx_swir_sensitivity import (  # noqa: E402
    JACCARD_MAX,
    PUBLISHED,
    cutoff_cell,
    evaluate_predictions,
    pin_cell,
)


def _disjoint_arrays():
    """5 cells: TIR-strong on the left, SWIR-high on the right, no overlap."""
    overlap = np.ones((1, 5), dtype=bool)
    combined = np.array([[4, 4, 0, 0, 0]], dtype=np.int16)
    depth = np.array([[0.01, 0.01, 0.01, 0.20, 0.20]], dtype=np.float32)
    return overlap, combined, depth


def test_cutoff_cell_disjoint_jaccard_is_zero():
    overlap, combined, depth = _disjoint_arrays()
    cell = cutoff_cell(overlap, combined, depth, score_min=3, pct=80)
    assert cell["jaccard"] == 0.0
    assert cell["n_both"] == 0
    assert cell["n_strong"] == 2
    assert cell["n_high"] == 2


def test_cutoff_cell_identical_jaccard_is_one():
    overlap = np.ones((1, 4), dtype=bool)
    combined = np.array([[4, 4, 0, 0]], dtype=np.int16)
    depth = np.array([[0.20, 0.20, 0.01, 0.01]], dtype=np.float32)
    cell = cutoff_cell(overlap, combined, depth, score_min=3, pct=50)
    assert cell["jaccard"] == 1.0
    assert cell["n_both"] == 2


def test_published_cell_is_in_the_grid():
    assert PUBLISHED == (3, 90)
    assert JACCARD_MAX == 0.20


def test_bodie_pins_are_wired_as_ls_control():
    sys.path.insert(0, str(SCRIPTS))
    from gemx_swir_pilot import SITE_PINS

    assert "bodie" in SITE_PINS
    assert SITE_PINS["bodie"]["town_label"] == "Bodie"
    assert SITE_PINS["bodie"]["target_label"] == "Bodie Bluff"


def test_pin_both_is_coincidence():
    overlap = np.ones((1, 2), dtype=bool)
    strong = np.array([[True, False]])
    high = np.array([[True, True]])
    depth = np.array([[0.2, 0.2]])
    pin = pin_cell(overlap, strong, high, depth, 0, 0, "town")
    assert pin["both"] is True
    pin2 = pin_cell(overlap, strong, high, depth, 0, 1, "ridge")
    assert pin2["both"] is False


def test_evaluate_predictions_passes_on_complementary_grid():
    cell = {
        "score_min": 3,
        "swir_pct": 90,
        "jaccard": 0.05,
        "lift": 1.2,
        "pins": {
            "town": {
                "in_strong_zone": False,
                "alunite_high": True,
                "both": False,
            },
            "target": {
                "in_strong_zone": False,
                "alunite_high": True,
                "both": False,
            },
        },
    }
    goldfield = {
        "site_id": "goldfield",
        "grid": [cell],
    }
    mary = {
        "site_id": "marysvale",
        "grid": [{
            **cell,
            "pins": {
                "town": {
                    "in_strong_zone": True,
                    "alunite_high": False,
                    "both": False,
                },
                "target": {
                    "in_strong_zone": False,
                    "alunite_high": True,
                    "both": False,
                },
            },
        }],
    }
    verdict = evaluate_predictions([goldfield, mary])
    assert verdict["jaccard_stays_below_0.20"] is True
    assert verdict["no_coincidence_flip"] is True
    assert verdict["published_pin_story_holds"] is True
    assert verdict["cutoff_fragile"] is False


def test_evaluate_predictions_fails_on_high_jaccard():
    cell = {
        "score_min": 3,
        "swir_pct": 90,
        "jaccard": 0.25,
        "lift": 3.0,
        "pins": {
            "town": {
                "in_strong_zone": False,
                "alunite_high": True,
                "both": False,
            },
            "target": {
                "in_strong_zone": False,
                "alunite_high": True,
                "both": False,
            },
        },
    }
    goldfield = {"site_id": "goldfield", "grid": [cell]}
    mary = {
        "site_id": "marysvale",
        "grid": [{
            **cell,
            "jaccard": 0.01,
            "pins": {
                "town": {
                    "in_strong_zone": True,
                    "alunite_high": False,
                    "both": False,
                },
                "target": {
                    "in_strong_zone": False,
                    "alunite_high": True,
                    "both": False,
                },
            },
        }],
    }
    verdict = evaluate_predictions([goldfield, mary])
    assert verdict["jaccard_stays_below_0.20"] is False
    assert verdict["cutoff_fragile"] is True
