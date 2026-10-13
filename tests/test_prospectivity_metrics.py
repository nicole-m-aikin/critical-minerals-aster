"""Site-mean metrics never silently pool across sites as the headline."""

from __future__ import annotations

import numpy as np

from critical_minerals_aster.prospectivity.metrics import (
    evaluate_site,
    site_mean_metrics,
)


def test_site_mean_averages_per_site_not_pixels():
    # Site A: perfect separation on few cells; site B: chance on many cells
    a = evaluate_site(
        np.array([1, 0, 0, 0]),
        np.array([0.99, 0.1, 0.05, 0.02]),
        site_id="a",
    )
    b = evaluate_site(
        np.array([1, 0] * 50),
        np.array([0.5] * 100),
        site_id="b",
    )
    summary = site_mean_metrics([a, b])
    assert summary["n_sites"] == 2
    # Headline is mean of site metrics, not pixel-weighted
    expected = np.nanmean([a["pr_auc"], b["pr_auc"]])
    assert abs(summary["site_mean_pr_auc"] - expected) < 1e-9
