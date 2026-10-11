import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.spectral import (  # noqa: E402
    CLOUD_SCORE_WEIGHT,
    extract_granule_id,
    granule_cloud_cover,
    granule_day_night,
    score_granule,
    select_granule,
    select_mosaic_granules,
)


class MockGranule:
    def __init__(
        self, gid: str, lons, lats, bands=(10, 11, 12, 13, 14), cloud=None,
        day_night=None, size_mb=10.0, month=None,
    ):
        self._gid = gid
        self._size = size_mb
        self._urls = [f"https://x/{gid}_TIR_B{b}.tif" for b in bands]
        self._umm = {
            "SpatialExtent": {
                "HorizontalSpatialDomain": {
                    "Geometry": {
                        "GPolygons": [
                            {
                                "Boundary": {
                                    "Points": [
                                        {"Longitude": lo, "Latitude": la}
                                        for lo, la in zip(lons, lats)
                                    ]
                                }
                            }
                        ]
                    }
                }
            }
        }
        if cloud is not None:
            self._umm["CloudCover"] = cloud
        if day_night is not None:
            self._umm["DataGranule"] = {"DayNightFlag": day_night}
        if month is not None:
            self._umm["TemporalExtent"] = {"SingleDateTime": f"2015-{month:02d}-15T18:00:00Z"}

    def data_links(self):
        return self._urls

    def size(self):
        return self._size

    def __getitem__(self, key):
        if key == "umm":
            return self._umm
        raise KeyError(key)


BBOX = (-118.1, 41.8, -117.3, 42.4)
GID_A = "AST_L1T_00407232010184946_20250705074029"
GID_B = "AST_L1T_00407232010184954_20250705074032"


def test_override_selects_exact_granule():
    g_a = MockGranule(GID_A, [-118.5, -117.0, -117.0, -118.5], [41.5, 41.5, 42.5, 42.5])
    g_b = MockGranule(GID_B, [-118.0, -117.2, -117.2, -118.0], [41.0, 41.0, 42.0, 42.0])
    picked = select_granule([g_a, g_b], BBOX, granule_id_override=GID_A)
    assert extract_granule_id(picked) == GID_A


def test_auto_pick_higher_coverage():
    g_a = MockGranule(
        GID_A,
        [-118.0, -117.5, -117.5, -118.0],
        [41.9, 41.9, 42.1, 42.1],
        bands=(10, 11),
    )
    g_b = MockGranule(
        GID_B,
        [-118.5, -117.0, -117.0, -118.5],
        [41.5, 41.5, 42.5, 42.5],
    )
    picked = select_granule([g_a, g_b], BBOX, granule_id_override=None)
    assert extract_granule_id(picked) == GID_B


FULL_LONS = [-118.5, -117.0, -117.0, -118.5]
FULL_LATS = [41.5, 41.5, 42.5, 42.5]
GID_C = "AST_L1T_00408082010184950_20250705080000"


def test_cloud_cover_parsing_is_defensive():
    assert granule_cloud_cover(MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=12)) == 12.0
    assert granule_cloud_cover(MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud="7")) == 7.0
    assert granule_cloud_cover(MockGranule(GID_A, FULL_LONS, FULL_LATS)) is None
    for bad in (-1, 150, "n/a", float("nan")):
        assert granule_cloud_cover(MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=bad)) is None
    assert granule_cloud_cover(object()) is None


def test_score_penalizes_cloud_but_keeps_signature():
    clear = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=0)
    cloudy = MockGranule(GID_B, FULL_LONS, FULL_LATS, cloud=100)
    unknown = MockGranule(GID_C, FULL_LONS, FULL_LATS)
    cov_c, comp_c, bands_c = score_granule(clear, BBOX)
    cov_x, comp_x, bands_x = score_granule(cloudy, BBOX)
    _, comp_u, _ = score_granule(unknown, BBOX)
    assert cov_c == pytest.approx(cov_x) == pytest.approx(1.0)
    assert bands_c == bands_x == 5
    assert comp_c == pytest.approx(15.0)
    assert comp_c - comp_x == pytest.approx(CLOUD_SCORE_WEIGHT)
    assert comp_u == pytest.approx(comp_c)
    assert score_granule(cloudy, BBOX, cloud_weight=0.0)[1] == pytest.approx(comp_c)


def test_lower_cloud_wins_at_similar_coverage():
    # A covers the full bbox at 25% cloud; B covers ~96% of it with 2% cloud.
    g_a = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=25)
    g_b = MockGranule(GID_B, [-118.5, -117.33, -117.33, -118.5], [41.5, 41.5, 42.5, 42.5], cloud=2)
    assert score_granule(g_b, BBOX)[0] == pytest.approx(0.9625, abs=1e-3)
    picked = select_granule([g_a, g_b], BBOX)
    assert extract_granule_id(picked) == GID_B


def test_coverage_still_dominates_large_gaps():
    g_a = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=20)
    g_b = MockGranule(GID_B, [-118.5, -117.7, -117.7, -118.5], [41.5, 41.5, 42.5, 42.5], cloud=0)
    picked = select_granule([g_a, g_b], BBOX)
    assert extract_granule_id(picked) == GID_A


def test_granules_above_max_cloud_are_skipped():
    g_cloudy = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=45)
    g_partial = MockGranule(
        GID_B, [-118.5, -117.7, -117.7, -118.5], [41.5, 41.5, 42.5, 42.5], cloud=5
    )
    assert score_granule(g_cloudy, BBOX)[1] > score_granule(g_partial, BBOX)[1]
    picked = select_granule([g_cloudy, g_partial], BBOX)
    assert extract_granule_id(picked) == GID_B
    picked = select_granule([g_cloudy, g_partial], BBOX, max_cloud_pct=50)
    assert extract_granule_id(picked) == GID_A
    picked = select_granule([g_cloudy, g_partial], BBOX, max_cloud_pct=None)
    assert extract_granule_id(picked) == GID_A


def test_unknown_cloud_is_not_skipped():
    g_cloudy = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=80)
    g_unknown = MockGranule(GID_B, [-118.5, -117.7, -117.7, -118.5], [41.5, 41.5, 42.5, 42.5])
    picked = select_granule([g_cloudy, g_unknown], BBOX)
    assert extract_granule_id(picked) == GID_B


def test_falls_back_to_least_bad_cloudy_granule_with_warning():
    g_a = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=90)
    g_b = MockGranule(GID_B, FULL_LONS, FULL_LATS, cloud=40)
    with pytest.warns(UserWarning, match="falling back"):
        picked = select_granule([g_a, g_b], BBOX)
    assert extract_granule_id(picked) == GID_B


def test_override_ignores_cloud_cutoff():
    g_a = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=95)
    g_b = MockGranule(GID_B, FULL_LONS, FULL_LATS, cloud=0)
    picked = select_granule([g_a, g_b], BBOX, granule_id_override=GID_A)
    assert extract_granule_id(picked) == GID_A


GID_D = "AST_L1T_00408092010184950_20250705080001"
WEST_LONS = [-118.5, -117.7, -117.7, -118.5]
EAST_LONS = [-117.7, -117.0, -117.0, -117.7]


def test_day_night_parsing():
    assert granule_day_night(MockGranule(GID_A, FULL_LONS, FULL_LATS, day_night="Day")) == "Day"
    assert granule_day_night(MockGranule(GID_A, FULL_LONS, FULL_LATS, day_night="NIGHT")) == "Night"
    assert granule_day_night(MockGranule(GID_A, FULL_LONS, FULL_LATS, day_night="Both")) is None
    assert granule_day_night(MockGranule(GID_A, FULL_LONS, FULL_LATS)) is None
    assert granule_day_night(object()) is None


def test_clear_day_beats_clearer_night():
    night = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=0, day_night="Night")
    day = MockGranule(GID_B, WEST_LONS, FULL_LATS, cloud=20, day_night="Day")
    assert extract_granule_id(select_granule([night, day], BBOX)) == GID_B
    assert extract_granule_id(select_granule([night, day], BBOX, prefer_day=False)) == GID_A


def test_clear_night_beats_cloudy_day_with_warning():
    night = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=5, day_night="Night")
    day = MockGranule(GID_B, FULL_LONS, FULL_LATS, cloud=60, day_night="Day")
    with pytest.warns(UserWarning, match="night_clear"):
        picked = select_granule([night, day], BBOX)
    assert extract_granule_id(picked) == GID_A


def test_mosaic_covers_bbox_with_clear_day_pair_reference_first():
    west = MockGranule(GID_A, WEST_LONS, FULL_LATS, cloud=2, day_night="Day")
    east = MockGranule(GID_B, EAST_LONS, FULL_LATS, cloud=10, day_night="Day")
    full_night = MockGranule(GID_C, FULL_LONS, FULL_LATS, cloud=0, day_night="Night")
    duplicate_west = MockGranule(GID_D, WEST_LONS, FULL_LATS, cloud=15, day_night="Day")
    sel = select_mosaic_granules([full_night, east, duplicate_west, west], BBOX)
    assert sel.tier == "day_clear"
    assert [r["granule_id"] for r in sel.records] == [GID_A, GID_B]
    assert sel.union_coverage == pytest.approx(1.0)


def test_mosaic_never_mixes_day_and_night():
    west_day = MockGranule(GID_A, WEST_LONS, FULL_LATS, cloud=0, day_night="Day")
    east_night = MockGranule(GID_B, EAST_LONS, FULL_LATS, cloud=0, day_night="Night")
    full_night = MockGranule(GID_C, FULL_LONS, FULL_LATS, cloud=10, day_night="Night")
    with pytest.warns(UserWarning, match="night_clear"):
        sel = select_mosaic_granules([west_day, east_night, full_night], BBOX)
    assert {r["day_night"] for r in sel.records} == {"Night"}
    assert sel.union_coverage == pytest.approx(1.0)


def test_mosaic_excludes_cloudy_and_oversize():
    cloudy = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=70, day_night="Day")
    big = MockGranule(GID_B, FULL_LONS, FULL_LATS, cloud=0, day_night="Day", size_mb=500)
    clear = MockGranule(GID_C, WEST_LONS, FULL_LATS, cloud=0, day_night="Day")
    with pytest.warns(UserWarning, match="coverage"):
        sel = select_mosaic_granules([cloudy, big, clear], BBOX, max_size_mb=100)
    assert sel.tier == "day_cloudy"
    assert [r["granule_id"] for r in sel.records] == [GID_A]
    sel = select_mosaic_granules([cloudy, big, clear], BBOX, min_acceptable_coverage=0.5)
    assert sel.tier == "day_clear"
    assert [r["granule_id"] for r in sel.records][0] == GID_B


def test_mosaic_prefers_summer_then_falls_back_to_any_month():
    winter = MockGranule(GID_A, FULL_LONS, FULL_LATS, cloud=0, day_night="Day", month=1)
    summer = MockGranule(GID_B, FULL_LONS, FULL_LATS, cloud=15, day_night="Day", month=7)
    sel = select_mosaic_granules([winter, summer], BBOX)
    assert sel.tier == "day_clear"
    assert sel.records[0]["granule_id"] == GID_B
    assert sel.records[0]["month"] == 7
    with pytest.warns(UserWarning, match="day_clear_offseason"):
        sel = select_mosaic_granules([winter], BBOX)
    assert sel.records[0]["granule_id"] == GID_A
    sel = select_mosaic_granules([winter, summer], BBOX, preferred_months=None)
    assert sel.records[0]["granule_id"] == GID_A


def test_fallback_keeps_day_when_night_barely_better():
    day = MockGranule(GID_A, [-118.5, -117.4, -117.4, -118.5], FULL_LATS, cloud=0, day_night="Day")
    night = MockGranule(GID_B, [-118.5, -117.37, -117.37, -118.5], FULL_LATS, cloud=0, day_night="Night")
    with pytest.warns(UserWarning, match="day_clear"):
        sel = select_mosaic_granules([night, day], BBOX, min_acceptable_coverage=0.99)
    assert sel.tier == "day_clear"
    with pytest.warns(UserWarning, match="night_clear"):
        sel = select_mosaic_granules(
            [night, day], BBOX, min_acceptable_coverage=0.99, fallback_coverage_tolerance=0.0
        )
    assert sel.tier == "night_clear"
