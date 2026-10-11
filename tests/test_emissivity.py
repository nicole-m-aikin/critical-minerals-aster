"""AST_05 emissivity inputs: QA mask, Ninomiya indices, granule loading, priority stacking."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from critical_minerals_aster.config import ClassificationParams, SiteConfig
from critical_minerals_aster.pipeline import _match_to_reference, build_emissivity_index_mosaic
from critical_minerals_aster.spectral import (
    acquisition_token,
    alteration_indices,
    ast05_qa_good,
    load_ast05_granule,
)


def test_qa_good_bits():
    qa = np.array(
        [
            0b0000_00_00,  # good, clear, far
            0b0000_00_01,  # good, clear, near -> kept at default adjacency
            0b0000_00_10,  # adjacency 2 -> rejected
            0b0000_01_00,  # thin cloud
            0b0000_10_00,  # thick cloud
            0b0001_00_00,  # suspect quality
            144,  # 1001 0000: bad (border / fill)
        ],
        dtype=np.uint8,
    )
    assert ast05_qa_good(qa).tolist() == [True, True, False, False, False, False, False]
    assert ast05_qa_good(qa, max_adjacency=0).tolist()[:2] == [True, False]


def test_alteration_indices_definitions():
    b10, b11, b12, b13, b14 = (np.array([v]) for v in (0.90, 0.95, 0.92, 0.96, 0.97))
    qi, ci, mi = alteration_indices(b10, b11, b12, b13, b14)
    assert qi[0] == pytest.approx(0.95**2 / (0.90 * 0.92))
    assert ci[0] == pytest.approx(0.96 / 0.97)
    assert mi[0] == pytest.approx(0.92 / 0.96)


def test_alteration_indices_zero_denominator_is_nan():
    z = np.array([0.0])
    o = np.array([1.0])
    qi, ci, mi = alteration_indices(z, o, o, z, z)
    assert np.isnan(qi[0]) and np.isnan(ci[0]) and np.isnan(mi[0])


def test_acquisition_token_matches_across_products():
    l1t = "AST_L1T_00405302014181013_20250705120000"
    ast05 = "AST_05_00405302014181013_20250728080701"
    assert acquisition_token(l1t) == acquisition_token(ast05) == "05302014181013"
    assert acquisition_token("bisbee_mosaic") is None
    assert acquisition_token(None) is None


CRS = "EPSG:32612"
RES = 90.0


def _write_granule(d: Path, gid: str, origin: tuple[float, float], shape, emis, qa_value=0):
    transform = from_origin(origin[0], origin[1], RES, RES)
    profile = dict(driver="GTiff", height=shape[0], width=shape[1], count=1, crs=CRS,
                   transform=transform)
    for b, val in emis.items():
        arr = np.full(shape, int(round(val * 1000)), dtype=np.uint16)
        with rasterio.open(d / f"{gid}_SRE_TIR_B{b}.tif", "w", dtype="uint16", **profile) as dst:
            dst.write(arr, 1)
    qa = np.full(shape, qa_value, dtype=np.uint8)
    with rasterio.open(d / f"{gid}_SRE_TIR_QA_DataPlane.tif", "w", dtype="uint8", **profile) as dst:
        dst.write(qa, 1)
    return transform


def test_load_ast05_granule_masks_qa_and_range(tmp_path):
    gid = "AST_05_00405302014181013_20250728080701"
    _write_granule(tmp_path, gid, (500000, 3500000), (4, 4),
                   {10: 0.9, 11: 0.95, 12: 0.92, 13: 0.96, 14: 0.97})
    with rasterio.open(tmp_path / f"{gid}_SRE_TIR_QA_DataPlane.tif", "r+") as dst:
        qa = dst.read(1)
        qa[0, 0] = 0b0000_10_00
        dst.write(qa, 1)
    with rasterio.open(tmp_path / f"{gid}_SRE_TIR_B11.tif", "r+") as dst:
        b = dst.read(1)
        b[1, 1] = 0
        b[2, 2] = 300  # emissivity 0.3: outside valid range
        dst.write(b, 1)
    bands, qa_good, *_ = load_ast05_granule(tmp_path, gid)
    assert not qa_good[0, 0]
    for arr in bands.values():
        assert np.isnan(arr[0, 0]) and np.isnan(arr[1, 1]) and np.isnan(arr[2, 2])
        assert np.isfinite(arr[3, 3])
    assert bands[13][3, 3] == pytest.approx(0.96)


def test_match_to_reference_uses_overlap():
    rng = np.random.default_rng(0)
    ref = rng.normal(1.0, 0.01, (40, 40))
    sec = ref * 2 + 5
    ref_partial = ref.copy()
    ref_partial[:, 30:] = np.nan
    matched, how = _match_to_reference(sec, ref_partial)
    assert how == "overlap"
    np.testing.assert_allclose(matched, ref, atol=1e-9)


def _site(bbox) -> SiteConfig:
    return SiteConfig(
        id="synthetic", name="synthetic", bbox_wgs84=bbox, granule_id=None,
        temporal_start="2010-01-01", temporal_end="2020-12-31",
        classification=ClassificationParams(),
    )


def test_priority_stack_fills_gaps_only(tmp_path):
    from rasterio.warp import transform_bounds

    x0, y0 = 500000.0, 3500000.0
    shape = (40, 40)
    emis_a = {10: 0.90, 11: 0.95, 12: 0.92, 13: 0.96, 14: 0.97}
    emis_b = {10: 0.85, 11: 0.95, 12: 0.90, 13: 0.94, 14: 0.97}
    ga = "AST_05_00405302014181013_20250728080701"
    gb = "AST_05_00406012014181013_20250728080702"
    _write_granule(tmp_path, ga, (x0, y0), shape, emis_a)
    _write_granule(tmp_path, gb, (x0, y0), shape, emis_b)
    with rasterio.open(tmp_path / f"{ga}_SRE_TIR_QA_DataPlane.tif", "r+") as dst:
        qa = dst.read(1)
        qa[:, :10] = 0b0000_10_00  # cloud over the left quarter of the reference
        dst.write(qa, 1)
    bbox = transform_bounds(CRS, "EPSG:4326", x0 + 5 * RES, y0 - 35 * RES, x0 + 35 * RES, y0 - 5 * RES)
    out = tmp_path / "out"
    stats = build_emissivity_index_mosaic(
        _site(tuple(bbox)), [ga, gb], tmp_path, out, "synthetic_mosaic",
        repo_root=tmp_path, apply_surface_masks=False,
    )
    with rasterio.open(out / "synthetic_mosaic_ratio_carbonate.tif") as src:
        ci = src.read(1)
    with rasterio.open(out / "synthetic_mosaic_cls_valid.tif") as src:
        valid = src.read(1).astype(bool)
    assert valid.mean() > 0.95
    # Reference wherever it is clear; the second granule only fills the cloud gap
    # (constant images have no variance, so they are left unmatched).
    ci_a, ci_b = 0.96 / 0.97, 0.94 / 0.97
    from_a = np.isclose(ci, ci_a, rtol=1e-5) & valid
    from_b = np.isclose(ci, ci_b, rtol=1e-5) & valid
    assert (from_a | from_b)[valid].all()
    assert 0.1 < from_b.sum() / valid.sum() < 0.25
    assert stats["valid_frac_of_bbox"] > 0.95
    assert {m["granule_id"] for m in stats["matching"]} == {gb}
    for name in ("silica", "mafic"):
        assert (out / f"synthetic_mosaic_cls_{name}.tif").is_file()
