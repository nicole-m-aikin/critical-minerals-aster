"""Tests for finer-scale geology loaders (offline; mocks / synthetic GeoDataFrames)."""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch

import geopandas as gpd
import pandas as pd
import pytest
from shapely.geometry import box

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.fine_geology import (  # noqa: E402
    fine_units_cache_path,
    load_fine_units,
    map_row_to_generalized_lith,
    normalize_units_lithology,
)

BBOX = (-115.0, 37.0, -114.5, 37.5)


def test_fine_units_cache_path_naming(tmp_path: Path):
    p = fine_units_cache_path("pioche", "nbmg_100k", tmp_path)
    assert p == tmp_path / "data" / "lithology" / "pioche_nbmg_100k_units.geojson"
    assert fine_units_cache_path("karnes_county", "tx_250k", tmp_path).name == (
        "karnes_county_tx_250k_units.geojson"
    )
    assert fine_units_cache_path("courtland_gleeson", "sgmc", tmp_path).name == (
        "courtland_gleeson_sgmc_units.geojson"
    )


@pytest.mark.parametrize(
    "props,expected",
    [
        # SGMC-style already present
        ({"GENERALIZED_LITH": "Sedimentary, carbonate"}, "Sedimentary, carbonate"),
        # NBMG OneGeology lithology strings
        ({"lithology": "limestone"}, "Sedimentary, carbonate"),
        ({"lithology": "dolomite"}, "Sedimentary, carbonate"),
        ({"lithology": "sandstone"}, "Sedimentary, clastic"),
        ({"lithology": "alluvium"}, "Unconsolidated, undifferentiated"),
        ({"lithology": "granite"}, "Igneous, intrusive"),
        ({"lithology": "rhyolite"}, "Igneous, volcanic"),
        ({"lithology": "schist"}, "Metamorphic, schist"),
        # Texas description / long name
        ({"description": "Clay, silt, sand, and gravel"}, "Unconsolidated, undifferentiated"),
        ({"long_nm": "Edwards Limestone"}, "Sedimentary, carbonate"),
        ({"rockunit_cd": "Qal", "description": "Alluvium"}, "Unconsolidated, undifferentiated"),
        # AZGS / generic local map fields
        ({"LITH": "limestone and dolomite"}, "Sedimentary, carbonate"),
        ({"UNIT_NAME": "Escabrosa Limestone"}, "Sedimentary, carbonate"),
        ({"name": "Quaternary alluvium"}, "Unconsolidated, undifferentiated"),
        # Unknown → Unknown (maps to "other" via generalize_lithology)
        ({"lithology": "mystery rock"}, "Unknown"),
        ({}, "Unknown"),
    ],
)
def test_map_row_to_generalized_lith(props, expected):
    assert map_row_to_generalized_lith(props) == expected


def test_normalize_units_adds_generalized_lith_and_lith_class():
    gdf = gpd.GeoDataFrame(
        {"lithology": ["limestone", "granite"]},
        geometry=[box(0, 0, 1, 1), box(1, 0, 2, 1)],
        crs="EPSG:4326",
    )
    out = normalize_units_lithology(gdf)
    assert list(out["GENERALIZED_LITH"]) == [
        "Sedimentary, carbonate",
        "Igneous, intrusive",
    ]
    assert list(out["lith_class"]) == ["carbonate_chemical", "intrusive"]


def test_auto_falls_back_to_sgmc_when_fine_source_empty(tmp_path: Path):
    """NV site prefers nbmg_100k; empty fine fetch → SGMC cache path / loader."""
    sgmc = gpd.GeoDataFrame(
        {"GENERALIZED_LITH": ["Sedimentary, carbonate"], "STATE": ["NV"], "SGMC_LABEL": ["Pl"]},
        geometry=[box(*BBOX)],
        crs="EPSG:4326",
    )
    empty = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")

    with patch(
        "critical_minerals_aster.fine_geology.fetch_nbmg_units",
        return_value=empty,
    ) as mock_nbmg, patch(
        "critical_minerals_aster.fine_geology.load_sgmc_units",
        return_value=sgmc,
    ) as mock_sgmc:
        out = load_fine_units(
            "pioche",
            BBOX,
            "auto",
            tmp_path,
            refresh=True,
            site_name="Pioche, NV",
        )
    mock_nbmg.assert_called_once()
    mock_sgmc.assert_called_once()
    assert len(out) == 1
    assert out.iloc[0]["GENERALIZED_LITH"] == "Sedimentary, carbonate"
    assert out.attrs.get("geology_source") == "sgmc_fallback_from_nbmg"
    # Empty fine fetch is not cached; fallback returns SGMC units.
    assert not fine_units_cache_path("pioche", "nbmg_100k", tmp_path).is_file()


def test_explicit_nbmg_empty_falls_back_to_sgmc(tmp_path: Path):
    sgmc = gpd.GeoDataFrame(
        {"GENERALIZED_LITH": ["Sedimentary, carbonate"], "STATE": ["NV"], "SGMC_LABEL": ["Pl"]},
        geometry=[box(*BBOX)],
        crs="EPSG:4326",
    )
    empty = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    with patch(
        "critical_minerals_aster.fine_geology.fetch_nbmg_units",
        return_value=empty,
    ), patch(
        "critical_minerals_aster.fine_geology.load_sgmc_units",
        return_value=sgmc,
    ):
        out = load_fine_units("pioche", BBOX, "nbmg_100k", tmp_path, refresh=True)
    assert out.attrs.get("geology_source") == "sgmc_fallback_from_nbmg"
    assert len(out) == 1


def test_load_fine_units_sgmc_delegates_and_caches_compatibly(tmp_path: Path):
    sgmc = gpd.GeoDataFrame(
        {"GENERALIZED_LITH": ["Igneous, volcanic"], "STATE": ["AZ"], "SGMC_LABEL": ["Tv"]},
        geometry=[box(-110.1, 31.7, -109.9, 31.9)],
        crs="EPSG:4326",
    )
    with patch(
        "critical_minerals_aster.fine_geology.load_sgmc_units",
        return_value=sgmc,
    ) as mock_sgmc:
        out = load_fine_units(
            "courtland_gleeson",
            (-110.1, 31.7, -109.9, 31.9),
            "sgmc",
            tmp_path,
        )
    mock_sgmc.assert_called_once()
    assert "lith_class" in out.columns
    assert out.iloc[0]["lith_class"] == "volcanic"


def test_azgs_local_reads_path_and_caches(tmp_path: Path):
    src = tmp_path / "CourtlandGleeson.gpkg"
    # Write a tiny GeoJSON and treat it as the local geology path (gpd reads both).
    local = tmp_path / "azgs_units.geojson"
    gpd.GeoDataFrame(
        {"UNIT_NAME": ["Escabrosa Limestone"], "LITH": ["limestone"]},
        geometry=[box(-110.05, 31.8, -109.95, 31.9)],
        crs="EPSG:4326",
    ).to_file(local, driver="GeoJSON")

    out = load_fine_units(
        "courtland_gleeson",
        (-110.08, 31.67, -109.58, 32.11),
        "azgs_local",
        tmp_path,
        geology_path=local,
        refresh=True,
    )
    cache = fine_units_cache_path("courtland_gleeson", "azgs_local", tmp_path)
    assert cache.is_file()
    assert len(out) == 1
    assert out.iloc[0]["GENERALIZED_LITH"] == "Sedimentary, carbonate"
    assert out.iloc[0]["lith_class"] == "carbonate_chemical"
    # unused but documents the expected download name
    assert src.name == "CourtlandGleeson.gpkg"


def test_preferred_auto_source_by_site():
    from critical_minerals_aster.fine_geology import preferred_fine_source

    assert preferred_fine_source("pioche", "Pioche, NV") == "nbmg_100k"
    assert preferred_fine_source("karnes_county", "Karnes County, TX") == "tx_250k"
    assert preferred_fine_source("courtland_gleeson", "Courtland-Gleeson, AZ") == "sgmc"
    assert preferred_fine_source("mcdermitt", "McDermitt Caldera, NV/OR") == "nbmg_100k"
