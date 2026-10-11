import numpy as np
import pandas as pd

from critical_minerals_aster.mrds import validate_mrds


def _frame(rows):
    cols = ["dep_id", "site_name", "latitude", "longitude", "dev_stat", "score", "commod1"]
    return pd.DataFrame(rows, columns=cols)


def test_drops_missing_and_out_of_range_coords():
    df = _frame(
        [
            (1, "A Mine", 40.0, -115.0, "Producer", "B", "Gold"),
            (2, "B Mine", np.nan, -115.0, "Producer", "B", "Gold"),
            (3, "C Mine", 95.0, -115.0, "Producer", "B", "Gold"),
            (4, "D Mine", 40.0, "not-a-number", "Producer", "B", "Gold"),
        ]
    )
    out, rep = validate_mrds(df)
    assert list(out["dep_id"]) == [1]
    assert rep["missing_coords"] == 2
    assert rep["out_of_range_coords"] == 1


def test_drops_processing_plants():
    df = _frame(
        [
            (1, "Smelter", 40.0, -115.0, "Plant", "C", "Copper"),
            (2, "Pit", 40.1, -115.1, "Past Producer", "C", "Copper"),
        ]
    )
    out, rep = validate_mrds(df)
    assert list(out["dep_id"]) == [2]
    assert rep["non_deposit_plant"] == 1


def test_duplicate_keeps_best_score_then_dev_stat():
    df = _frame(
        [
            (1, "Big Mine ", 40.0, -115.0, "Occurrence", "D", "Gold"),
            (2, "big mine", 40.0, -115.0, "Producer", "B", "Gold"),
            (3, "Big Mine", 40.0, -115.0, "Past Producer", "B", "Silver"),
        ]
    )
    out, rep = validate_mrds(df)
    assert list(out["dep_id"]) == [2]
    assert rep["duplicate_record"] == 2


def test_shared_coords_flagged_not_dropped_and_index_preserved():
    df = _frame(
        [
            (1, "North Adit", 40.0, -115.0, "Prospect", "A", "Gold"),
            (2, "South Adit", 40.0, -115.0, "Occurrence", "E", "Gold"),
            (3, "Elsewhere", 41.0, -116.0, "Producer", "C", "Gold"),
        ],
    )
    df.index = [10, 20, 30]
    out, rep = validate_mrds(df)
    assert list(out.index) == [10, 20, 30]
    assert out.loc[[10, 20], "qc_shared_coords"].all()
    assert not out.loc[30, "qc_shared_coords"]
    assert out["qc_score_high"].tolist() == [True, False, True]
    assert out["qc_producer"].tolist() == [False, False, True]
    assert rep["n_output"] == 3
    assert rep["n_input"] == 3
