"""Snapshot test for ``pipeline.run_site()`` — pins the output bytes.

Two tests:

* ``test_pipeline_synthetic_snapshot`` (always-on): runs ``run_site`` against
  committed minimal synthetic fixtures in ``tests/fixtures/synthetic_aster/``
  (5 tiny 16x16 GeoTIFFs for ASTER TIR bands 10-14, plus a 4-row MRDS CSV
  chosen to exercise all three groupby loops in ``compute_site_summary``)
  and byte-matches the three always-written outputs (strong-zones GeoJSON,
  summary CSV, normalised provenance JSON) against committed snapshots in
  ``tests/fixtures/snapshots/``.
* ``test_pipeline_mcdermitt_real`` (``@pytest.mark.slow``, default-skipped):
  smoke-runs the pipeline against the real mcdermitt ASTER scenes if they're
  already on disk.  Sanity assertions only (file exists, summary not empty);
  not a byte-equality check.

To accept an intentional output change, re-run with ``UPDATE_SNAPSHOTS=1``:

    UPDATE_SNAPSHOTS=1 .venv/bin/pytest tests/test_pipeline_snapshot.py

then ``git diff tests/fixtures/snapshots/`` to review the deltas before
committing.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from difflib import unified_diff
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC))

from critical_minerals_aster.config import (  # noqa: E402
    ClassificationParams,
    SiteConfig,
    load_site_by_id,
)
from critical_minerals_aster.paths import site_paths_for  # noqa: E402
from critical_minerals_aster.pipeline import run_site  # noqa: E402


TESTS_DIR = Path(__file__).resolve().parent
SYNTH_FIXTURE = TESTS_DIR / "fixtures" / "synthetic_aster"
SNAPSHOTS = TESTS_DIR / "fixtures" / "snapshots"


# --- snapshot helpers -------------------------------------------------------

# Provenance fields that vary per-run (timestamp) or per-environment
# (git_commit, python build, package versions).  Replaced with this sentinel
# before snapshot diff so the diff still catches changes to *structural*
# fields like site_id, granule_id, n_zones, raster_bbox_wgs84.
_NORMALIZED_KEYS = {"timestamp_utc", "git_commit", "python"}
_NORMALIZED_SENTINEL = "<normalized>"


def _update_mode() -> bool:
    return os.environ.get("UPDATE_SNAPSHOTS") == "1"


def _diff_capped(expected: str, actual: str, expected_name: str, actual_name: str) -> str:
    diff_lines = list(
        unified_diff(
            expected.splitlines(keepends=True),
            actual.splitlines(keepends=True),
            fromfile=f"expected/{expected_name}",
            tofile=f"actual/{actual_name}",
            n=3,
        )
    )
    capped = diff_lines[:60]
    if len(diff_lines) > 60:
        capped.append(f"... ({len(diff_lines) - 60} more lines)\n")
    return "".join(capped)


def assert_text_snapshot(actual_path: Path, expected_path: Path) -> None:
    actual = actual_path.read_text()
    if _update_mode():
        expected_path.parent.mkdir(parents=True, exist_ok=True)
        expected_path.write_text(actual)
        return
    expected = expected_path.read_text()
    if actual != expected:
        diff = _diff_capped(expected, actual, expected_path.name, actual_path.name)
        pytest.fail(
            f"Snapshot mismatch for {expected_path.name}.\n{diff}\n"
            "Re-run with UPDATE_SNAPSHOTS=1 to accept the new output."
        )


def _normalize_provenance(text: str) -> str:
    data = json.loads(text)
    for key in _NORMALIZED_KEYS:
        if key in data:
            data[key] = _NORMALIZED_SENTINEL
    if isinstance(data.get("packages"), dict):
        data["packages"] = {k: _NORMALIZED_SENTINEL for k in data["packages"]}
    return json.dumps(data, indent=2)


def assert_provenance_snapshot(actual_path: Path, expected_path: Path) -> None:
    normalized = _normalize_provenance(actual_path.read_text())
    if _update_mode():
        expected_path.parent.mkdir(parents=True, exist_ok=True)
        expected_path.write_text(normalized)
        return
    expected = expected_path.read_text()
    if normalized != expected:
        diff = _diff_capped(expected, normalized, expected_path.name, actual_path.name)
        pytest.fail(
            f"Snapshot mismatch for {expected_path.name}.\n{diff}\n"
            "Re-run with UPDATE_SNAPSHOTS=1 to accept the new output."
        )


# --- tests ------------------------------------------------------------------


def test_pipeline_synthetic_snapshot(tmp_path: Path) -> None:
    """End-to-end pipeline run on committed minimal synthetic fixtures."""
    shutil.copytree(SYNTH_FIXTURE / "data", tmp_path / "data")

    site = SiteConfig(
        id="synth",
        name="Synthetic Test Site",
        bbox_wgs84=(-120.0, 35.0, -119.9, 35.1),
        granule_id="AST_L1T_99999999_99999999",
        layout="flat",
        classification=ClassificationParams(),
    )
    run_site(site, tmp_path, skip_figures=True)

    paths = site_paths_for(site, tmp_path)
    assert_text_snapshot(
        paths.strong_zones_geojson, SNAPSHOTS / "synth_strong_zones.geojson"
    )
    assert_text_snapshot(
        paths.site_summary_csv, SNAPSHOTS / "synth_summary.csv"
    )
    assert_provenance_snapshot(
        paths.site_provenance_json, SNAPSHOTS / "synth_provenance.json"
    )


@pytest.mark.slow
def test_pipeline_mcdermitt_real(tmp_path: Path) -> None:
    """Smoke test against real mcdermitt ASTER scenes when present on disk.

    Sandbox approach: symlink the real inputs (ASTER bands, MRDS, DEM,
    structure layers, the sites/ dir) from the live repo into ``tmp_path``,
    then point ``run_site`` at ``tmp_path``.  Outputs land under ``tmp_path``
    and never overwrite the live repo's results.
    """
    real_repo = Path(__file__).resolve().parents[1]
    real_aster_dir = real_repo / "data" / "sites" / "mcdermitt" / "aster"
    if not list(real_aster_dir.glob("*_TIR_B10.tif")):
        pytest.skip(
            f"no real ASTER scenes at {real_aster_dir}; download mcdermitt to run"
        )

    sandbox = tmp_path
    (sandbox / "sites").symlink_to(real_repo / "sites")
    (sandbox / "data" / "sites" / "mcdermitt").mkdir(parents=True)
    (sandbox / "data" / "sites" / "mcdermitt" / "aster").symlink_to(real_aster_dir)
    for candidate in (
        real_repo / "data" / "mrds.csv",
        real_repo / "data" / "mrds" / "mrds.csv",
    ):
        if candidate.is_file():
            (sandbox / candidate.relative_to(real_repo).parent).mkdir(
                parents=True, exist_ok=True
            )
            (sandbox / candidate.relative_to(real_repo)).symlink_to(candidate)
            break
    real_dem = real_repo / "data" / "dem" / "mcdermitt"
    if real_dem.is_dir():
        (sandbox / "data" / "dem").mkdir(parents=True, exist_ok=True)
        (sandbox / "data" / "dem" / "mcdermitt").symlink_to(real_dem)
    real_structures = real_repo / "data" / "structures"
    if real_structures.is_dir():
        (sandbox / "data" / "structures").symlink_to(real_structures)

    site = load_site_by_id("mcdermitt", sandbox / "sites")
    summary = run_site(site, sandbox, skip_figures=True)

    paths = site_paths_for(site, sandbox)
    assert paths.strong_zones_geojson.exists()
    assert paths.site_summary_csv.exists()
    assert paths.site_provenance_json.exists()
    assert not summary.empty
    assert {"site_id", "granule_id"} <= set(summary.columns)
