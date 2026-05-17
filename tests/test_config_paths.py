import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.config import (  # noqa: E402
    SiteConfig,
    load_site_by_id,
    load_site_config,
    list_site_ids,
    search_bbox,
)
from critical_minerals_aster.paths import site_paths_for  # noqa: E402


REPO = Path(__file__).resolve().parents[1]


def test_list_site_ids():
    ids = list_site_ids(REPO / "sites")
    assert "mcdermitt" in ids
    assert "silver_peak" in ids


def test_flat_layout_paths():
    # No site yaml uses layout="flat" anymore, so synthesize a config
    # to keep coverage of the flat branch in SitePaths.
    site = SiteConfig(
        id="dummy_flat",
        name="Dummy Flat Site",
        bbox_wgs84=(-120.0, 40.0, -119.0, 41.0),
        layout="flat",
    )
    paths = site_paths_for(site, REPO)
    assert paths.aster_dir == REPO / "data" / "aster"
    assert paths.vectors_dir == REPO / "data" / "vectors"
    assert paths.figures_dir == REPO / "figures"


def test_silver_peak_nested_paths():
    site = load_site_config(REPO / "sites" / "silver_peak.yaml")
    paths = site_paths_for(site, REPO)
    assert paths.aster_dir == REPO / "data" / "sites" / "silver_peak" / "aster"
    assert paths.figures_dir == REPO / "figures" / "sites" / "silver_peak"


def test_search_bbox_buffer():
    site = load_site_by_id("mcdermitt", REPO / "sites")
    site.buffer_deg = 0.1
    bbox = search_bbox(site)
    assert bbox[0] < site.bbox_wgs84[0]
