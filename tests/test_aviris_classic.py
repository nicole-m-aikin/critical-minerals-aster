"""Tests for dual-sensor AVIRIS collection selection and Classic ENVI I/O."""

from __future__ import annotations

import json
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from critical_minerals_aster.aviris import (  # noqa: E402
    MAX_WAVELENGTH_DELTA_NM,
    parse_envi_header,
    read_envi_rfl_subset,
    wavelength_targets_ok,
)
from gemx_site_select import select_census_collection  # noqa: E402


def test_select_census_collection_prefers_higher_coverage():
    """Catches choosing AVIRIS-5 when Classic covers more of the bbox."""
    long_rows = [
        {"site_id": "a", "collection": "aviris5_l2a_rfl", "coverage_frac": 0.4},
        {"site_id": "a", "collection": "avirisc_l2_rfl", "coverage_frac": 0.9},
        {"site_id": "b", "collection": "aviris5_l2a_rfl", "coverage_frac": 0.0},
        {"site_id": "b", "collection": "avirisc_l2_rfl", "coverage_frac": 0.0},
    ]
    assert select_census_collection("a", long_rows) == (
        "avirisc_l2_rfl",
        0.9,
        "AVIRIS-Classic",
    )
    assert select_census_collection("b", long_rows) is None


def test_select_census_collection_tiebreaks_to_aviris5():
    """Catches Classic winning when coverages are equal."""
    long_rows = [
        {"site_id": "x", "collection": "avirisc_l2_rfl", "coverage_frac": 1.0},
        {"site_id": "x", "collection": "aviris5_l2a_rfl", "coverage_frac": 1.0},
    ]
    assert select_census_collection("x", long_rows)[0] == "aviris5_l2a_rfl"


def test_wavelength_targets_ok_rejects_sparse_sampling():
    """Catches inventing band depths when nearest bands are too far."""
    good = np.array([2120.0, 2165.0, 2175.0, 2200.0, 2206.0, 2210.0, 2270.0])
    assert wavelength_targets_ok(good, (2120, 2165, 2175, 2200, 2206, 2210, 2270))
    bad = np.array([2100.0, 2150.0, 2250.0])
    assert not wavelength_targets_ok(bad, (2120, 2165, 2175, 2200, 2206, 2210, 2270))
    assert MAX_WAVELENGTH_DELTA_NM == 8.0


def _write_envi_pair(tmp_path: Path, *, bands: int = 7, lines: int = 4, samples: int = 5):
    """Tiny BSQ float32 ENVI cube with map info and wavelength list."""
    waves = [2120.0, 2165.0, 2175.0, 2200.0, 2206.0, 2210.0, 2270.0][:bands]
    data = np.arange(bands * lines * samples, dtype=np.float32).reshape(bands, lines, samples)
    data[0, 0, 0] = -9999.0
    bin_path = tmp_path / "scene.rfl.bin"
    hdr_path = tmp_path / "scene.rfl.hdr"
    data.tofile(bin_path)
    hdr = (
        "ENVI\n"
        "description = {test}\n"
        f"samples = {samples}\n"
        f"lines = {lines}\n"
        f"bands = {bands}\n"
        "header offset = 0\n"
        "file type = ENVI Standard\n"
        "data type = 4\n"
        "interleave = bsq\n"
        "byte order = 0\n"
        "map info = {UTM, 1.000, 1.000, 500000.000, 4200000.000, 20.000, 20.000, "
        "11, North, WGS-84, units=Meters}\n"
        "wavelength units = Nanometers\n"
        "wavelength = {" + ", ".join(f"{w:.1f}" for w in waves) + "}\n"
        "data ignore value = -9999\n"
    )
    hdr_path.write_text(hdr)
    return hdr_path, bin_path, waves, data


def test_parse_envi_header_reads_map_info_and_wavelengths(tmp_path):
    """Catches broken ENVI map-info / wavelength parsing."""
    hdr_path, _, waves, _ = _write_envi_pair(tmp_path)
    meta = parse_envi_header(hdr_path)
    assert meta["samples"] == 5
    assert meta["lines"] == 4
    assert meta["bands"] == 7
    assert meta["interleave"] == "bsq"
    assert meta["data_type"] == 4
    assert meta["wavelengths_nm"] == pytest.approx(waves)
    assert meta["transform"].c == pytest.approx(500000.0 - 10.0)  # pixel center → corner
    assert "32611" in meta["crs"] or "UTM" in meta["crs"]


def test_read_envi_rfl_subset_returns_nearest_bands(tmp_path):
    """Catches wrong band indexing or nodata handling on Classic cubes."""
    hdr_path, bin_path, waves, data = _write_envi_pair(tmp_path)
    targets = (2120.0, 2200.0, 2270.0)
    stack, actual, transform, crs = read_envi_rfl_subset(bin_path, targets, hdr_path=hdr_path)
    assert stack.shape == (3, 4, 5)
    assert actual.tolist() == [2120.0, 2200.0, 2270.0]
    assert np.isnan(stack[0, 0, 0])
    assert stack[1, 0, 1] == pytest.approx(data[3, 0, 1])
    assert transform.a == pytest.approx(20.0)
    assert "32611" in crs or "UTM" in crs


def test_read_envi_rfl_subset_requires_matching_hdr(tmp_path):
    """Catches silent read when the .hdr sibling is missing."""
    _, bin_path, _, _ = _write_envi_pair(tmp_path)
    with pytest.raises(FileNotFoundError, match="hdr"):
        read_envi_rfl_subset(bin_path, (2200.0,), hdr_path=tmp_path / "missing.hdr")


def test_select_census_collection_good_tier_prefers_aviris5():
    """Catches Classic winning when both sensors already have good coverage."""
    long_rows = [
        {"site_id": "a", "collection": "aviris5_l2a_rfl", "coverage_frac": 0.83},
        {"site_id": "a", "collection": "avirisc_l2_rfl", "coverage_frac": 1.0},
    ]
    assert select_census_collection("a", long_rows)[0] == "aviris5_l2a_rfl"


def test_read_envi_rfl_subset_bil(tmp_path):
    """Catches rejecting Classic BIL cubes that ORNL actually ships."""
    waves = [2120.0, 2165.0, 2175.0, 2200.0, 2206.0, 2210.0, 2270.0]
    bands, lines, samples = 7, 3, 4
    bil = np.zeros((lines, bands, samples), dtype=np.float32)
    for b in range(bands):
        bil[:, b, :] = b + 0.1
    bil[0, 0, 0] = -9999.0
    bin_path = tmp_path / "scene_rfl.bin"
    hdr_path = tmp_path / "scene_rfl.hdr"
    bil.tofile(bin_path)
    hdr = [
        "ENVI",
        f"samples = {samples}",
        f"lines = {lines}",
        f"bands = {bands}",
        "header offset = 0",
        "data type = 4",
        "interleave = bil",
        "byte order = 0",
        "map info = {UTM, 1.000, 1.000, 500000.000, 4200000.000, 20.000, 20.000, 11, North, WGS-84, units=Meters}",
        "wavelength = {" + ", ".join(f"{w:.1f}" for w in waves) + "}",
        "data ignore value = -9999",
    ]
    hdr_path.write_text("\n".join(hdr) + "\n")
    stack, actual, transform, crs = read_envi_rfl_subset(
        bin_path, (2120.0, 2200.0, 2270.0), hdr_path=hdr_path
    )
    assert stack.shape == (3, 3, 4)
    assert np.isnan(stack[0, 0, 0])
    assert stack[1, 1, 1] == pytest.approx(3.1)
