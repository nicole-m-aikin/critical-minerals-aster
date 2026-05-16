import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.spectral import (  # noqa: E402
    extract_granule_id,
    score_granule,
    select_granule,
)


class MockGranule:
    def __init__(self, gid: str, lons, lats, bands=(10, 11, 12, 13, 14)):
        self._gid = gid
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

    def data_links(self):
        return self._urls

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
