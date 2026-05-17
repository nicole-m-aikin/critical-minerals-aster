#!/usr/bin/env python
"""
Filtered significance tests: TIR-detectable deposit types only.

Compares unfiltered vs. TIR-relevant-only MRDS significance to answer
whether the whole MRDS catalog dilutes the signal.

TIR-detectable systems = those with surface alteration mineralogy visible in
ASTER TIR (silica, carbonate, mafic alteration zones):
  - Porphyry Cu-Mo-Au, Climax-type     → potassic / argillic / propylitic alteration
  - Magmatic REE, Hybrid Magmatic REE  → carbonate signature (B13/B12)
  - Mafic Magmatic                     → mafic ratio (B12/B13)
  - Reduced Intrusion-Related          → silicic / argillic halo
  - Porphyry Sn                        → greisen / argillic

TIR-invisible (excluded from filtered test):
  - Orogenic, Carlin-type              → no or buried alteration
  - Placer, Basin Brine, Evaporites    → detrital or brine-hosted
  - Meteoric Recharge/Convection       → no alteration signature
  - Marine Chemocline, Seafloor        → sediment-hosted
  - Unknown System                     → unclassified
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
    reclassify_mrds_mineral_system,
)
from critical_minerals_aster.significance import (
    coverage_fraction,
    run_binomial,
    run_permutation,
)

REPO_ROOT = Path(__file__).parent.parent
N_ITER = 10_000

TIR_DETECTABLE = {
    "Porphyry Cu-Mo-Au",
    "Climax-type",
    "Porphyry Sn",
    "Magmatic REE",
    "Hybrid Magmatic REE",
    "Mafic Magmatic",
    "Reduced Intrusion-Related",
}


def load_site_ids() -> list[str]:
    with open(REPO_ROOT / "sites" / "index.yaml") as f:
        return yaml.safe_load(f)["sites"]


def run_site(site_id: str, deposit_filter=None) -> dict | None:
    site = load_site_by_id(site_id, REPO_ROOT / "sites")
    paths = SitePaths(repo_root=REPO_ROOT, site=site)

    zones_path = paths.strong_zones_geojson
    if not zones_path.exists():
        return None

    zones = gpd.read_file(zones_path)
    if zones.empty:
        return None

    try:
        mrds_df = read_mrds_national(paths)
    except FileNotFoundError:
        return None

    local = filter_mrds_bbox(mrds_df, site.bbox_wgs84)
    if deposit_filter is not None:
        local = reclassify_mrds_mineral_system(local)
        local = local[local["mineral_system"].isin(deposit_filter)]

    if local.empty:
        return None

    deposits = mrds_to_points_gdf(local, zones.crs)
    joined, _, _ = spatial_join_deposits_zones(deposits, zones)
    inside = joined[joined["index_right"].notna()]

    n_dep = len(deposits)
    n_hit = inside.index.nunique()

    if n_dep == 0:
        return None

    p_cov = coverage_fraction(zones, site.bbox_wgs84)
    p_binom, expected = run_binomial(n_hit, n_dep, p_cov)
    p_perm = run_permutation(zones, site.bbox_wgs84, n_dep, n_hit, n_iter=N_ITER)

    return {
        "site_id": site_id,
        "n_deposits": n_dep,
        "n_hits": n_hit,
        "hit_rate_pct": round(n_hit / n_dep * 100, 1),
        "zone_coverage_pct": round(p_cov * 100, 2),
        "expected_hits_null": round(expected, 1),
        "p_binom": round(p_binom, 6),
        "p_perm": round(p_perm, 4),
    }


def sig_label(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "** "
    if p < 0.05:
        return "*  "
    return "   "


def main() -> None:
    site_ids = load_site_ids()
    unfiltered, filtered = [], []

    print(f"{'Site':25s}  {'─── All MRDS ───':>30s}    {'─── TIR-detectable only ───':>35s}")
    print(f"{'':25s}  {'n':>5} {'hits':>5} {'p_binom':>9}     {'n':>5} {'hits':>5} {'p_binom':>9}")
    print("─" * 90)

    for site_id in site_ids:
        r_all = run_site(site_id, deposit_filter=None)
        r_flt = run_site(site_id, deposit_filter=TIR_DETECTABLE)

        if r_all:
            unfiltered.append(r_all)
        if r_flt:
            filtered.append({**r_flt, "filter": "tir_detectable"})

        n_all = r_all["n_deposits"] if r_all else 0
        h_all = r_all["n_hits"] if r_all else 0
        p_all = r_all["p_binom"] if r_all else float("nan")

        n_flt = r_flt["n_deposits"] if r_flt else 0
        h_flt = r_flt["n_hits"] if r_flt else 0
        p_flt = r_flt["p_binom"] if r_flt else float("nan")

        arrow = ""
        if r_all and r_flt:
            if p_flt < p_all * 0.5:
                arrow = " ← stronger"
            elif p_flt > p_all * 2 and p_flt > 0.1:
                arrow = " ← weaker"

        print(
            f"{site_id:25s}  {n_all:5d} {h_all:5d} {p_all:8.4f}{sig_label(p_all)}"
            f"  {n_flt:5d} {h_flt:5d} {p_flt:8.4f}{sig_label(p_flt)}{arrow}"
        )

    # Save filtered results
    df_flt = pd.DataFrame(filtered)
    if not df_flt.empty:
        out = REPO_ROOT / "results" / "significance_tir_filtered.csv"
        df_flt.to_csv(out, index=False)
        print(f"\nFiltered results written to {out}")

    print("\n* p<0.05  ** p<0.01  *** p<0.001")
    print(f"TIR-detectable systems: {sorted(TIR_DETECTABLE)}")


if __name__ == "__main__":
    main()
