"""Inventory EarthAccess archive for all sites — shows what's blocked by filters.

For each site prints:
  - granules already passing both filters (would be mosaicked today)
  - granules blocked by size only (large bundles, good coverage)
  - granules blocked by coverage only or both
"""
from __future__ import annotations
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import earthaccess
earthaccess.login(strategy="netrc")

from critical_minerals_aster.config import load_site_by_id, list_site_ids, search_bbox
from critical_minerals_aster.spectral import score_granule

sites_dir = REPO_ROOT / "sites"
site_ids = list_site_ids(sites_dir)

print(f"{'SITE':<22} {'passing':>7} {'blocked-size':>12} {'max_cov_blocked':>15} {'current_cap':>11}")
print("-" * 72)

for sid in site_ids:
    site = load_site_by_id(sid, sites_dir)
    bbox = search_bbox(site)
    try:
        results = earthaccess.search_data(
            short_name="AST_L1T",
            bounding_box=bbox,
            temporal=(site.temporal_start, site.temporal_end),
            count=20,
        )
    except Exception as e:
        print(f"{sid:<22}  ERROR: {e}")
        continue

    passing = []       # pass both filters
    size_blocked = []  # good coverage but too large
    other = []         # low coverage (with or without size issue)

    for g in results:
        try:
            cov, _, bands = score_granule(g, site.bbox_wgs84)
            sz = g.size()
            if bands < 3:
                continue
            if cov > 0.30 and sz < site.max_bundle_mb:
                passing.append((cov, sz))
            elif cov > 0.30 and sz >= site.max_bundle_mb:
                size_blocked.append((cov, sz))
            else:
                other.append((cov, sz))
        except Exception:
            continue

    max_blocked_cov = max((c for c, _ in size_blocked), default=0.0)
    print(
        f"{sid:<22} {len(passing):>7} {len(size_blocked):>12} {max_blocked_cov:>15.2f} {site.max_bundle_mb:>11.0f} MB"
    )
