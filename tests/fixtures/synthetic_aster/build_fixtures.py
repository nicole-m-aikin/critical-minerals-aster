"""Generate synthetic ASTER TIR fixtures for the pipeline snapshot test.

This script is committed for transparency and reproducibility; it is NOT
auto-run by the test suite. The generated ``.tif`` files (and the synthetic
``data/mrds/mrds.csv``) are also committed.

Re-run only when you intentionally want to change the synthetic fixture
(e.g., adjust pixel patterns to exercise a new pipeline branch). After
re-running this, regenerate the committed snapshots:

    UPDATE_SNAPSHOTS=1 pytest tests/test_pipeline_snapshot.py

Run from this directory:

    python build_fixtures.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_bounds


GRANULE = "AST_L1T_99999999_99999999"
HEIGHT = WIDTH = 16
# (west, south, east, north) in WGS84 degrees.  Tiny patch matches the
# synthetic SiteConfig in the test.
BBOX = (-120.0, 35.0, -119.9, 35.1)


def main() -> None:
    here = Path(__file__).resolve().parent
    aster_dir = here / "data" / "aster"
    aster_dir.mkdir(parents=True, exist_ok=True)
    mrds_dir = here / "data" / "mrds"
    mrds_dir.mkdir(parents=True, exist_ok=True)

    # Distinct deterministic linear gradients per band so the three band
    # ratios (B13/B14, B13/B12, B12/B13) all vary across the grid and the
    # 70th/90th percentile thresholds land on real values instead of degenerate
    # constants.  Pipeline produces a non-empty mix of strong/moderate cells.
    i, j = np.indices((HEIGHT, WIDTH), dtype=np.float32)
    bands = {
        10: 100.0 + 1.0 * i + 1.0 * j,
        11: 100.0 + 2.0 * i + 1.0 * j,
        12: 100.0 + 5.0 * i + 2.0 * j,
        13: 100.0 + 3.0 * i + 7.0 * j,
        14: 100.0 + 4.0 * i + 4.0 * j,
    }
    transform = from_bounds(*BBOX, WIDTH, HEIGHT)

    for band_num, data in bands.items():
        out_path = aster_dir / f"{GRANULE}_TIR_B{band_num}.tif"
        with rasterio.open(
            out_path,
            "w",
            driver="GTiff",
            height=HEIGHT,
            width=WIDTH,
            count=1,
            dtype=np.float32,
            transform=transform,
            crs="EPSG:4326",
        ) as dst:
            dst.write(data, 1)
        print(f"wrote {out_path}")

    # Synthetic MRDS rows chosen to exercise all three groupby loops in
    # compute_site_summary (commodity_group, earth_mri_category,
    # mineral_system) AND both branches of is_critical_mineral.
    #
    # Each row's classification (with the rules as of this commit):
    #   1. lithium       → Lithium       / Battery Metals – Li/Brine / Lacustrine Evaporite   (commod1 tier)
    #   2. gold          → Gold/Silver   / Gold/Silver               / Orogenic               (commod1 tier)
    #   3. copper+model  → Other         / Base Metals               / Porphyry Cu-Mo-Au      (model tier — most-specific cascade)
    #   4. building stone → Stone        / Non-Critical              / Unknown System         (fallthrough; is_critical_mineral=False)
    #
    # The `model` column is required to exercise the model-tier branch of
    # reclassify_mrds_mineral_system's three-field cascade.
    mrds_csv = mrds_dir / "mrds.csv"
    mrds_csv.write_text(
        "longitude,latitude,commod1,model,dep_type\n"
        "-119.95,35.05,lithium,,\n"
        "-119.92,35.02,gold,,\n"
        "-119.98,35.08,copper,60c: porphyry cu,\n"
        "-119.94,35.06,building stone,,\n"
    )
    print(f"wrote {mrds_csv}")


if __name__ == "__main__":
    main()
