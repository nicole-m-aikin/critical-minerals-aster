#!/usr/bin/env python
"""
Phase 8b -- fix stale cached zone polygons found during threshold-sensitivity
testing.

Phase 8a's sensitivity sweep reclassifies fresh from each site's cached
ratio mosaic (data/sites/{id}/aster/{id}_mosaic_ratio_*.tif) at the nominal
70/90/>=3 settings and, as a side effect, gives an independent check on
whether the CACHED zone polygons (data/sites/{id}/vectors/
strong_anomaly_zones.geojson) still match what that same ratio mosaic
would produce today. For most sites they agree closely (small pixel-level
noise). Two sites relevant to this project's reported significant-site
list disagree by far more than noise:

    site        fresh px hits   cached poly hits   ratio
    yerington   49              18                 2.6x
    iron_hill   33              20                 1.6x

Both are in the 6-site "continuous-only FDR-significant" group (Phase 3/4)
-- meaning their BINARY test result (not significant, per Phase 2) may have
been computed against an undercount of true hits. This script regenerates
ONLY these two sites' zone polygons from their current ratio mosaics (same
method as the rest of the pipeline: classify_percentiles -> combined_score
-> vectorize_strong_zones, using each site's own ClassificationParams), so
every phase from here on is self-consistent with the ratio mosaic actually
on disk. The old files are git-tracked, so this is a recoverable, reviewed
change (see `git diff` after running), not a destructive edit.

After regenerating, re-run in order:
    scripts/regenerate_site_summaries.py
    scripts/site_specific_null_significance.py
    scripts/phase3_monte_carlo_and_continuous_score.py
    scripts/phase4_fdr_correction.py
to propagate the correction through every downstream table.

Usage:
    conda run -n aster-minerals python scripts/phase8b_fix_stale_zones.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.classification import (
    classify_percentiles,
    combined_score,
    vectorize_strong_zones,
)
from critical_minerals_aster.spectral import load_ratio_mosaic

REPO_ROOT = Path(__file__).parent.parent
SITES_TO_FIX = ["yerington", "iron_hill"]


def main() -> None:
    for site_id in SITES_TO_FIX:
        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        paths = site_paths_for(site, REPO_ROOT)
        prefix = f"{site_id}_mosaic"

        silica, carbonate, mafic, _, transform, crs = load_ratio_mosaic(paths.aster_dir, prefix)
        cp = site.classification
        s_cls, _, _ = classify_percentiles(silica, cp.low_pct, cp.high_pct)
        c_cls, _, _ = classify_percentiles(carbonate, cp.low_pct, cp.high_pct)
        m_cls, _, _ = classify_percentiles(mafic, cp.low_pct, cp.high_pct)
        combined = combined_score(s_cls, c_cls, m_cls)

        old_path = paths.strong_zones_geojson
        import geopandas as gpd
        old_zones = gpd.read_file(old_path)
        old_area = float(old_zones["area_km2"].sum()) if len(old_zones) else 0.0

        new_zones = vectorize_strong_zones(combined, transform, crs, min_score=cp.strong_score_min)
        new_area = float(new_zones["area_km2"].sum()) if len(new_zones) else 0.0

        new_zones.to_file(old_path, driver="GeoJSON")
        print(f"{site_id:15s} zones regenerated: {len(old_zones)} polygons ({old_area:.2f} km2) "
              f"-> {len(new_zones)} polygons ({new_area:.2f} km2)")


if __name__ == "__main__":
    main()
