import sys
from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import LineString, Point

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.structure import (  # noqa: E402
    nearest_structure_distance_m,
    points_on_structure,
)


def test_nearest_structure_distance():
    structures = gpd.GeoDataFrame(
        geometry=[LineString([(0, 0), (10, 0)])],
        crs="EPSG:32611",
    )
    points = gpd.GeoDataFrame(
        geometry=[Point(5, 3)],
        crs="EPSG:32611",
    )
    dist = nearest_structure_distance_m(points, structures)
    assert dist.iloc[0] == pytest.approx(3.0)


def test_points_on_structure_within_buffer():
    structures = gpd.GeoDataFrame(
        geometry=[LineString([(0, 0), (100, 0)])],
        crs="EPSG:32611",
    )
    points = gpd.GeoDataFrame(
        geometry=[Point(50, 40), Point(50, 200)],
        crs="EPSG:32611",
    )
    on = points_on_structure(points, structures, buffer_m=50)
    assert on.iloc[0]
    assert not on.iloc[1]
