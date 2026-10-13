"""Unit tests for GEMx SWIR batch site selection and verdicts.  No cubes."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(ROOT / "src"))

from gemx_site_select import (  # noqa: E402
    select_census_collection,
    census_collection_by_site,
    JACCARD_OVERLAP,
    LIFT_COMPLEMENTARY,
    LIFT_OVERLAP,
    SITE_PINS,
    aviris5_coverage_by_site,
    bbox_centroid,
    correspondence_verdict,
    primary_fdr_sites,
    resolve_pins,
    select_sites,
    swir_covered_sites,
)


def test_select_sites_single_and_list():
    assert select_sites(repo_root=ROOT, site="eureka") == ["eureka"]
    assert select_sites(repo_root=ROOT, sites=["a", "b"]) == ["a", "b"]


def test_select_sites_rejects_mixed_modes():
    with pytest.raises(ValueError, match="choose only one"):
        select_sites(repo_root=ROOT, site="eureka", primary_fdr=True)


def test_select_sites_primary_fdr_matches_significance():
    ids = select_sites(repo_root=ROOT, primary_fdr=True)
    expected = primary_fdr_sites(ROOT)
    assert ids == expected
    assert len(ids) == 19


def test_select_sites_default_is_primary_fdr():
    ids = select_sites(repo_root=ROOT)
    assert ids == primary_fdr_sites(ROOT)


def test_select_sites_all_covered_uses_has_swir_csv():
    ids = select_sites(repo_root=ROOT, all_covered=True)
    from_csv = swir_covered_sites(ROOT)
    assert ids == from_csv
    assert len(ids) > 19  # many more than FDR-only
    assert "mcdermitt" in ids


def test_aviris5_coverage_positive_for_known_pilots():
    cov = aviris5_coverage_by_site(ROOT)
    assert cov.get("goldfield", 0) > 0
    assert cov.get("marysvale", 0) > 0
    assert cov.get("bodie", 0) > 0


def test_bbox_centroid():
    assert bbox_centroid((-2.0, 10.0, 2.0, 20.0)) == (0.0, 15.0)


def test_resolve_pins_uses_site_pins_table():
    site = SimpleNamespace(
        id="bodie",
        name="Bodie, CA",
        bbox_wgs84=(-119.1, 38.1, -118.9, 38.3),
    )
    pins = resolve_pins(site)
    assert pins["town"] == SITE_PINS["bodie"]["town"]
    assert pins["target_label"] == "Bodie Bluff"
    assert pins.get("pins_source") == "site_pins"


def test_resolve_pins_centroid_fallback_when_unpinned():
    site = SimpleNamespace(
        id="not_a_real_pin_site",
        name="Fake District, NV",
        bbox_wgs84=(-118.0, 38.0, -117.0, 39.0),
    )
    pins = resolve_pins(site, pins_table={})
    lon, lat = bbox_centroid(site.bbox_wgs84)
    assert pins["town"] == (lon, lat)
    assert pins["target"] == (lon, lat)
    assert pins["pins_source"] == "bbox_centroid"
    assert pins["town_label"] == "Fake District"


def test_resolve_pins_yaml_center_overrides_centroid(tmp_path):
    yaml_path = tmp_path / "fake_site.yaml"
    yaml_path.write_text(
        "id: fake_site\n"
        "name: Fake\n"
        "bbox_wgs84: [-1.0, 0.0, 1.0, 2.0]\n"
        "center_wgs84: [-0.25, 1.5]\n"
    )
    site = SimpleNamespace(
        id="fake_site",
        name="Fake",
        bbox_wgs84=(-1.0, 0.0, 1.0, 2.0),
    )
    pins = resolve_pins(site, pins_table={}, sites_dir=tmp_path)
    assert pins["town"] == (-0.25, 1.5)
    assert pins["pins_source"] == "yaml_center"


def test_verdict_no_coverage_and_error():
    assert correspondence_verdict(status="no_coverage") == "no_coverage"
    assert correspondence_verdict(swir_coverage=0.0) == "no_coverage"
    assert correspondence_verdict(status="error") == "error"


def test_verdict_insufficient_small_overlap_or_nan():
    assert (
        correspondence_verdict(
            jaccard_aloh=0.01,
            lift_aloh=0.5,
            n_overlap=10,
            n_tir_strong=5,
        )
        == "insufficient"
    )
    assert (
        correspondence_verdict(
            jaccard_aloh=None,
            lift_aloh=1.0,
            n_overlap=500,
            n_tir_strong=10,
        )
        == "insufficient"
    )


def test_verdict_overlapping_by_jaccard_or_lift():
    assert (
        correspondence_verdict(
            jaccard_aloh=JACCARD_OVERLAP,
            lift_aloh=1.0,
            n_overlap=500,
            n_tir_strong=20,
        )
        == "overlapping"
    )
    assert (
        correspondence_verdict(
            jaccard_aloh=0.05,
            lift_aloh=LIFT_OVERLAP,
            n_overlap=500,
            n_tir_strong=20,
        )
        == "overlapping"
    )


def test_verdict_complementary():
    assert (
        correspondence_verdict(
            jaccard_aloh=0.05,
            lift_aloh=LIFT_COMPLEMENTARY - 0.01,
            n_overlap=500,
            n_tir_strong=20,
        )
        == "complementary"
    )


def test_verdict_mixed_band_is_insufficient():
    """Jaccard low but lift in [1.1, 1.3) → inconclusive → insufficient."""
    assert (
        correspondence_verdict(
            jaccard_aloh=0.05,
            lift_aloh=1.2,
            n_overlap=500,
            n_tir_strong=20,
        )
        == "insufficient"
    )


def test_pilot_reexports_site_pins():
    from gemx_swir_pilot import SITE_PINS as pilot_pins

    assert "marysvale" in pilot_pins
    assert pilot_pins is SITE_PINS or pilot_pins == SITE_PINS


def test_select_census_collection_national_counts():
    """Catches drifting away from the 71-site full-SWIR census rule."""
    picked = census_collection_by_site(ROOT)
    assert len(picked) == 71
    n_av5 = sum(1 for k, _, _ in picked.values() if k == "aviris5_l2a_rfl")
    n_classic = sum(1 for k, _, _ in picked.values() if k == "avirisc_l2_rfl")
    assert n_av5 + n_classic == 71
    assert n_av5 >= 25
    assert n_classic >= 11
    # Courtland is Classic-only in CMR inventory.
    assert picked["courtland_gleeson"][0] == "avirisc_l2_rfl"
