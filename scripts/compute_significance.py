#!/usr/bin/env python
"""
Binomial and spatial permutation significance tests for all 15 sites.

Usage:
    conda run -n aster-minerals python scripts/compute_significance.py

Writes: results/significance.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import yaml
import geopandas as gpd
import pandas as pd

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import SitePaths
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    mrds_to_points_gdf,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.significance import (
    coverage_fraction,
    run_binomial,
    run_permutation,
)

REPO_ROOT = Path(__file__).parent.parent
N_ITER = 10_000


def load_site_ids() -> list[str]:
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        return yaml.safe_load(f)["sites"]


def main() -> None:
    site_ids = load_site_ids()
    rows: list[dict] = []

    for site_id in site_ids:
        print(f"{site_id:25s}", end=" ", flush=True)

        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        paths = SitePaths(repo_root=REPO_ROOT, site=site)

        zones_path = paths.strong_zones_geojson
        if not zones_path.exists():
            print("SKIP — no zones GeoJSON")
            continue

        zones = gpd.read_file(zones_path)
        if zones.empty:
            print("SKIP — empty zones")
            continue

        # Load deposits within the bbox
        try:
            mrds_df = read_mrds_national(paths)
        except FileNotFoundError as exc:
            print(f"SKIP — {exc}")
            continue

        local = filter_mrds_bbox(mrds_df, site.bbox_wgs84)
        deposits = mrds_to_points_gdf(local, zones.crs)
        joined, _, _ = spatial_join_deposits_zones(deposits, zones)
        inside = joined[joined["index_right"].notna()]

        n_dep = len(deposits)
        n_hit = inside.index.nunique()

        if n_dep == 0:
            print("SKIP — no MRDS deposits in bbox")
            continue

        p_cov = coverage_fraction(zones, site.bbox_wgs84)
        p_binom, expected = run_binomial(n_hit, n_dep, p_cov)
        p_perm = run_permutation(zones, site.bbox_wgs84, n_dep, n_hit, n_iter=N_ITER)

        sig_marker = "**" if p_binom < 0.01 else ("*" if p_binom < 0.05 else "  ")
        print(
            f"n={n_dep:4d}  hits={n_hit:3d}  cov={p_cov:5.1%}"
            f"  expected={expected:5.1f}  p_binom={p_binom:.4f}"
            f"  p_perm={p_perm:.4f}  {sig_marker}"
        )

        rows.append(
            {
                "site_id": site_id,
                "n_deposits": n_dep,
                "n_hits": n_hit,
                "hit_rate_pct": round(n_hit / n_dep * 100, 1),
                "zone_coverage_pct": round(p_cov * 100, 2),
                "expected_hits_null": round(expected, 1),
                "p_binom": round(p_binom, 6),
                "p_perm": round(p_perm, 4),
                "significant_binom_05": p_binom < 0.05,
                "significant_binom_01": p_binom < 0.01,
            }
        )

    df = pd.DataFrame(rows).sort_values("p_binom")
    out = REPO_ROOT / "results" / "significance.csv"
    df.to_csv(out, index=False)

    print(f"\n{'─'*80}")
    print(f"Results written to {out}\n")
    print(df.to_string(index=False))
    print("\n* p < 0.05   ** p < 0.01   (one-sided binomial, H₀: hit rate ≤ zone coverage)")


if __name__ == "__main__":
    main()
