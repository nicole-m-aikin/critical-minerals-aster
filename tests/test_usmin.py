"""Offline tests for critical_minerals_aster.usmin (synthetic geometries, no network)."""

import math

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Point, Polygon, box

from critical_minerals_aster.usmin import (
    disturbed_mask,
    nearest_feature_distance,
    osm_elements_to_gdf,
    overlap_fraction,
    parse_wfs_gml,
    points_in_mask,
    quad_cells,
    quad_coverage,
    subtract_mask,
)

CRS = "EPSG:32612"

GML = """<?xml version='1.0' encoding="UTF-8" ?>
<wfs:FeatureCollection xmlns:ms="http://mapserver.gis.umn.edu/mapserver"
   xmlns:gml="http://www.opengis.net/gml" xmlns:wfs="http://www.opengis.net/wfs">
  <gml:featureMember>
    <ms:polygons fid="polygons.1">
      <ms:geometry>
        <gml:Polygon srsName="EPSG:4267">
          <gml:outerBoundaryIs><gml:LinearRing>
            <gml:coordinates>-110.0,31.0 -109.9,31.0 -109.9,31.1 -110.0,31.1 -110.0,31.0</gml:coordinates>
          </gml:LinearRing></gml:outerBoundaryIs>
          <gml:innerBoundaryIs><gml:LinearRing>
            <gml:coordinates>-109.98,31.02 -109.92,31.02 -109.92,31.08 -109.98,31.08 -109.98,31.02</gml:coordinates>
          </gml:LinearRing></gml:innerBoundaryIs>
        </gml:Polygon>
      </ms:geometry>
      <ms:ftr_type>Tailings - Undifferentiated</ms:ftr_type>
      <ms:topo_name>Tombstone</ms:topo_name>
    </ms:polygons>
  </gml:featureMember>
  <gml:featureMember>
    <ms:points-high fid="points-high.2">
      <ms:geometry><gml:Point srsName="EPSG:4267"><gml:coordinates>-110.05,31.05</gml:coordinates></gml:Point></ms:geometry>
      <ms:ftr_type>Adit</ms:ftr_type>
    </ms:points-high>
  </gml:featureMember>
</wfs:FeatureCollection>
"""


def test_parse_wfs_gml_polygon_with_hole_and_point():
    gdf = parse_wfs_gml(GML)
    assert len(gdf) == 2
    assert gdf.crs.to_epsg() == 4267
    poly = gdf.geometry.iloc[0]
    assert poly.geom_type == "Polygon" and len(poly.interiors) == 1
    assert poly.area == pytest.approx(0.01 - 0.06 * 0.06)
    assert list(gdf["ftr_type"]) == ["Tailings - Undifferentiated", "Adit"]
    assert gdf.geometry.iloc[1].equals(Point(-110.05, 31.05))


def test_parse_wfs_gml_raises_on_service_exception():
    xml = ('<?xml version="1.0"?><ServiceExceptionReport version="1.2.0">'
           '<ServiceException code="InvalidParameterValue">bad</ServiceException>'
           '</ServiceExceptionReport>')
    with pytest.raises(RuntimeError, match="service exception"):
        parse_wfs_gml(xml)


def test_parse_wfs_gml_empty_collection_is_empty_not_error():
    xml = ('<?xml version="1.0"?><wfs:FeatureCollection xmlns:wfs="http://www.opengis.net/wfs" '
           'xmlns:gml="http://www.opengis.net/gml"/>')
    assert parse_wfs_gml(xml).empty


def test_disturbed_mask_polygons_only_ignores_points_at_zero_buffer():
    polys = [box(0, 0, 100, 100), box(50, 0, 150, 100)]
    pts = [Point(1000, 1000)]
    m = disturbed_mask(polys, pts, 0)
    assert m.area == pytest.approx(150 * 100)


def test_disturbed_mask_adds_point_buffers():
    m = disturbed_mask([box(0, 0, 10, 10)], [Point(1000, 1000)], 100)
    assert m.area == pytest.approx(100 + math.pi * 100**2, rel=0.01)


def test_disturbed_mask_empty_inputs():
    assert disturbed_mask(None, None, 100).is_empty
    assert disturbed_mask(gpd.GeoSeries([], crs=CRS), gpd.GeoSeries([], crs=CRS), 50).is_empty


def test_overlap_fraction_and_subtract_mask():
    zone = box(0, 0, 100, 100)
    mask = box(50, 0, 200, 100)
    assert overlap_fraction(zone, mask) == pytest.approx(0.5)
    assert overlap_fraction(zone, Polygon()) == 0.0
    assert overlap_fraction(Polygon(), mask) == 0.0
    assert subtract_mask(zone, mask).area == pytest.approx(5000)
    assert subtract_mask(zone, Polygon()).equals(zone)


def test_masked_null_stays_consistent():
    """Removing a mask that lies entirely outside zones raises p0 by the
    footprint shrinkage only: p0' = zone / (footprint - mask)."""
    footprint = box(0, 0, 1000, 1000)
    zones = box(0, 0, 100, 1000)
    mask = box(900, 0, 1000, 1000)
    p0 = zones.area / footprint.area
    p0m = subtract_mask(zones, mask).area / subtract_mask(footprint, mask).area
    assert p0 == pytest.approx(0.1)
    assert p0m == pytest.approx(100 / 900)


def test_nearest_feature_distance_points_and_polygons():
    feats = [Point(0, 0), box(1000, 0, 1100, 100)]
    targets = [Point(3, 4), Point(1050, 50), Point(1200, 50)]
    d = nearest_feature_distance(targets, feats)
    np.testing.assert_allclose(d, [5.0, 0.0, 100.0])


def test_nearest_feature_distance_no_reference_is_nan():
    d = nearest_feature_distance([Point(0, 0), Point(1, 1)], [])
    assert d.shape == (2,) and np.isnan(d).all()
    assert nearest_feature_distance([], [Point(0, 0)]).size == 0


def test_nearest_feature_distance_accepts_geodataframes():
    t = gpd.GeoDataFrame(geometry=[Point(0, 0)], crs=CRS)
    f = gpd.GeoDataFrame(geometry=[Point(0, 250)], crs=CRS)
    assert nearest_feature_distance(t, f)[0] == pytest.approx(250)


def test_points_in_mask():
    pts = gpd.GeoSeries([Point(5, 5), Point(50, 50), Point(10, 5)], crs=CRS)
    got = points_in_mask(pts, box(0, 0, 10, 10))
    assert got.tolist() == [True, False, True]
    assert points_in_mask(pts, Polygon()).tolist() == [False, False, False]


def test_quad_cells_and_coverage():
    cells = quad_cells((-110.2, 31.0, -110.0, 31.2))
    assert len(cells) == 4
    assert all(c.bounds[2] - c.bounds[0] == pytest.approx(0.125) for c in cells)
    fp = box(-110.2, 31.0, -110.0, 31.2)
    feats = gpd.GeoSeries([Point(-110.1, 31.1)], crs="EPSG:4326")
    n, n_hit = quad_coverage(fp, feats)
    assert n == 4 and n_hit == 1
    assert quad_coverage(fp, gpd.GeoSeries([], crs="EPSG:4326")) == (4, 0)


def test_osm_elements_to_gdf():
    elements = [
        {"type": "node", "id": 1, "lat": 31.0, "lon": -110.0, "tags": {"historic": "mine"}},
        {"type": "way", "id": 2, "tags": {"landuse": "quarry"},
         "geometry": [{"lat": 0, "lon": 0}, {"lat": 0, "lon": 1}, {"lat": 1, "lon": 1},
                      {"lat": 1, "lon": 0}, {"lat": 0, "lon": 0}]},
        {"type": "relation", "id": 3, "tags": {"man_made": "spoil_heap"},
         "members": [
             {"type": "way", "role": "outer",
              "geometry": [{"lat": 0, "lon": 10}, {"lat": 0, "lon": 14}, {"lat": 4, "lon": 14}]},
             {"type": "way", "role": "outer",
              "geometry": [{"lat": 4, "lon": 14}, {"lat": 4, "lon": 10}, {"lat": 0, "lon": 10}]},
             {"type": "way", "role": "inner",
              "geometry": [{"lat": 1, "lon": 11}, {"lat": 1, "lon": 12}, {"lat": 2, "lon": 12},
                           {"lat": 2, "lon": 11}, {"lat": 1, "lon": 11}]},
         ]},
    ]
    gdf = osm_elements_to_gdf(elements)
    assert list(gdf.geom_type) == ["Point", "Polygon", "Polygon"]
    assert list(gdf["ftr_type"]) == ["historic=mine", "landuse=quarry", "man_made=spoil_heap"]
    assert gdf.geometry.iloc[2].area == pytest.approx(16 - 1)
