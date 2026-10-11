#!/usr/bin/env python
"""
Phase 5 — MRDS clustering sensitivity (DBSCAN at 250/500/1000 m).

Motivation (docs/results.md): scipy.stats.binomtest treats
every MRDS row as an independent Bernoulli trial, but the Phase 1 spatial
audit found this is badly violated at exactly the sites carrying the
strongest claims -- e.g. Eureka's median nearest-neighbor distance between
"independent" critical-mineral deposits is 4.2 m. Counting the same mine
several times as several "independent" hits inflates the effective sample
size the binomial test sees, which mechanically shrinks p-values and can
manufacture significance from spatial autocorrelation rather than real
deposit-alteration co-location.

Method: DBSCAN (density-based clustering)
------------------------------------------
A point is a *core point* if at least `min_samples` other points lie within
Euclidean distance `eps` of it. Core points within `eps` of each other join
into one cluster; points reachable through a chain of core points join too;
points that are neither core points nor reachable from one are labelled
"noise" (cluster label -1) and are NOT merged with anything.

Why DBSCAN is suitable here: it needs no pre-specified number of clusters
(unlike k-means) and treats isolated, genuinely independent occurrences as
singletons rather than forcing them into a cluster (unlike single-linkage
hierarchical clustering with a fixed cut height, which absorbs stragglers).
That "leave true singletons alone" behavior matters because most MRDS
occurrences in this dataset *are* independent -- only a minority are
tightly duplicated records of the same mine.

Parameters:
  eps         -- 250, 500, 1000 m, matching the existing buffer-sensitivity
                 convention in scripts/mrds_location_sensitivity.py, so the
                 two analyses are directly comparable radius-for-radius.
  min_samples -- 2. MRDS routinely splits one historical mine into several
                 rows (different commodities, workings, report vintages),
                 so even a *pair* of coincident points indicates
                 non-independence. A higher min_samples (e.g. 5) would
                 leave common pairwise duplicates labelled "noise" and
                 uncorrected -- min_samples=2 is the smallest value that
                 fully captures pairwise clustering.

CRS: coordinates are projected to each site's own local UTM zone
(zones.crs -- already the CRS the pipeline stores zone polygons in) BEFORE
clustering. DBSCAN's eps is a raw Euclidean distance with no notion of
units; in geographic coordinates (EPSG:4326) one degree of longitude spans
a different physical distance depending on latitude (~111 km * cos(lat)),
so an eps expressed in degrees would mean a different real-world radius at
every site and even at different latitudes within one site's bbox. Local
UTM keeps distance distortion under ~0.04% within a single site's ~50x50 km
extent -- far below the precision this analysis needs.

Effective sample size: each DBSCAN cluster (label >= 0) collapses to ONE
effective unit; each noise point (label -1) remains its own unit.

A unit's hit status is the in-zone status of its MEDOID member (the point
with the smallest total distance to every other point in its cluster) --
NOT "any point in the cluster hits" (OR-aggregation). This distinction
matters a great deal and was caught empirically before trusting any output
of this script: an initial OR-rule implementation caused 8 of 45 sites to
flip from non-significant to significant purely from clustering, because
OR-ing k correlated Bernoulli(p0) draws has null-expected hit probability
1-(1-p0)^k > p0 -- strictly higher than the single-point null p0 this
script still compares against. Comparing an OR-inflated observed rate to an
uncorrected p0 null systematically manufactures significance as radius (and
therefore typical cluster size) grows, exactly backwards from clustering
sensitivity's purpose. Picking one representative point's true hit status
per cluster is a fair single Bernoulli trial under the same p0 as any other
point, so it does not have this bias -- this is genuine deduplication (in
the sense used throughout docs/results.md: "collapsing
clustered points into one representative point"), not a merged-evidence
rule. See docs/results.md for the full derivation.

Null after declustering (bedrock-vs-cover stratified, the Phase 2 primary)
---------------------------------------------------------------------------
Every deposit carries its own null probability p_i = zone coverage fraction
of its stratum (bedrock / cover / unmapped; lithology.cover_strata). After
declustering, each effective unit keeps the p_i of the point whose hit
status it uses (the medoid for a cluster, itself for a noise point), and
E = Σ p_i over effective units; p_binomial is the exact Poisson-binomial
tail. Re-using the site-level p0 would be wrong here: declustering changes
the bedrock/cover mix of the deposit set (dense mining camps sit in bedrock),
so the expected hits must be recomputed from the surviving units.
p0_null = E / n_effective; enrichment = hits / E. The site-level geometric
null is kept as p0_geometric / p_binomial_geometric for comparison.

Outputs
-------
    results/phase5_clustering_sensitivity.csv     -- per site x radius table
    figures/phase5_clustering_stability.png       -- p-value stability figure

Usage:
    conda run -n aster-minerals python scripts/phase5_clustering_sensitivity.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import box
from shapely.ops import unary_union
from shapely.validation import make_valid
from sklearn.cluster import DBSCAN

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.lithology import (
    assign_cover_stratum,
    cover_strata,
    stratified_null,
    stratum_zone_fractions,
)
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.significance import run_binomial

REPO_ROOT = Path(__file__).parent.parent
RADII_M = (0, 250, 500, 1000)  # 0 = no clustering (raw), for comparison
MIN_SAMPLES = 2


def _medoid_index(coords: np.ndarray) -> int:
    """Index (within *coords*) of the point with smallest summed distance
    to every other point in the group -- a deterministic, unbiased choice
    of "the" representative location for a cluster."""
    if len(coords) == 1:
        return 0
    diffs = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diffs ** 2).sum(axis=-1))
    return int(dist.sum(axis=1).argmin())


def effective_units(deposits: gpd.GeoDataFrame, eps_m: float) -> pd.DataFrame:
    """Collapse deposits to DBSCAN clusters + noise singletons.

    Returns one row per effective unit with the ``in_zone`` and ``p_null``
    of its representative point -- the cluster MEDOID, or the noise point
    itself (see module docstring for why OR-aggregation is biased and was
    rejected).
    """
    cols = ["in_zone", "p_null"]
    if eps_m == 0 or len(deposits) <= 1:
        return deposits[cols].reset_index(drop=True)

    coords = np.column_stack([deposits.geometry.x.values, deposits.geometry.y.values])
    labels = DBSCAN(eps=eps_m, min_samples=MIN_SAMPLES).fit_predict(coords)

    df = deposits.copy()
    df["_cluster"] = labels
    units = []
    for label, grp in df.groupby("_cluster"):
        if label == -1:
            # noise points are NOT merged -- each is its own independent unit
            units.append(grp[cols])
        else:
            grp_coords = np.column_stack([grp.geometry.x.values, grp.geometry.y.values])
            units.append(grp[cols].iloc[[_medoid_index(grp_coords)]])
    return pd.concat(units, ignore_index=True)


def effective_n_hits(deposits: gpd.GeoDataFrame, eps_m: float) -> tuple[int, int]:
    """(effective_n, effective_hits) after medoid declustering."""
    d = deposits if "p_null" in deposits.columns else deposits.assign(p_null=0.0)
    units = effective_units(d, eps_m)
    return len(units), int(units["in_zone"].sum())


def main() -> None:
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    rows = []

    for _, ss_row in ss.iterrows():
        site_id = ss_row["site_id"]
        p0_geo = float(ss_row["p0_geometric"])
        print(f"{site_id:22s}", end=" ", flush=True)

        site = load_site_by_id(site_id, REPO_ROOT / "sites")
        paths = site_paths_for(site, REPO_ROOT)
        zones_path = paths.strong_zones_geojson
        prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
        if not zones_path.exists() or not prov_path.exists():
            print("SKIP")
            continue
        zones = gpd.read_file(zones_path)
        with open(prov_path) as f:
            prov = json.load(f)
        raster_bbox = tuple(prov["raster_bbox_wgs84"])
        footprint = None
        if prov.get("tir_footprint_wgs84") is not None:
            try:
                footprint = gpd.GeoDataFrame.from_features(
                    prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
                ).to_crs(zones.crs)
            except Exception:
                footprint = None

        mrds = read_mrds_national(paths)
        local = filter_mrds_bbox(mrds, raster_bbox)
        deposits = mrds_to_points_gdf(local, zones.crs)
        if footprint is not None and len(footprint):
            try:
                deposits = gpd.clip(deposits, footprint)
            except Exception:
                pass
        deposits = reclassify_mrds_earth_mri(deposits)
        deposits = deposits[deposits["earth_mri_category"] != "Non-Critical"].copy()

        if deposits.empty:
            print("SKIP — no critical deposits")
            continue

        joined, hits, _ = spatial_join_deposits_zones(deposits, zones)
        hit_ids = set(joined[joined["index_right"].notna()].index.unique())
        deposits["in_zone"] = deposits.index.isin(hit_ids)

        if footprint is not None and len(footprint):
            fp_geom = make_valid(unary_union(footprint.geometry.values))
        else:
            fp_geom = gpd.GeoSeries([box(*raster_bbox)], crs="EPSG:4326").to_crs(zones.crs).iloc[0]
        strata = cover_strata(site_id, fp_geom, zones.crs, bbox_wgs84=raster_bbox)
        frac = stratum_zone_fractions(strata, make_valid(unary_union(zones.geometry.values)))
        deposits["stratum"] = assign_cover_stratum(deposits, strata)
        deposits["p_null"] = deposits["stratum"].map(frac["p"]).astype(float)

        for radius in RADII_M:
            units = effective_units(deposits, radius)
            n_eff, hits_eff = len(units), int(units["in_zone"].sum())
            strat = stratified_null(units["p_null"].to_numpy(), hits_eff)
            p_geo, _ = run_binomial(hits_eff, n_eff, p0_geo)
            hit_rate = round(hits_eff / n_eff * 100, 2) if n_eff else 0.0
            rows.append({
                "site_id": site_id,
                "radius_m": radius,
                "n_effective": n_eff,
                "hits_effective": hits_eff,
                "hit_rate_pct": hit_rate,
                "p0_null": round(strat["p0_eff"], 4),
                "expected_hits_null": round(strat["expected"], 2),
                "enrichment": round(strat["enrichment"], 2) if strat["expected"] > 0 else float("nan"),
                "p_binomial": strat["p"],
                "sig_05": strat["p"] < 0.05,
                "p0_geometric": p0_geo,
                "enrichment_geometric": round((hit_rate / 100) / p0_geo, 2) if p0_geo > 0 and n_eff else float("nan"),
                "p_binomial_geometric": p_geo,
                "sig_05_geometric": p_geo < 0.05,
            })

        r0 = [r for r in rows if r["site_id"] == site_id and r["radius_m"] == 0][0]
        r1000 = [r for r in rows if r["site_id"] == site_id and r["radius_m"] == 1000][0]
        print(f"n(0m)={r0['n_effective']:4d} -> n(1000m)={r1000['n_effective']:4d}  "
              f"p(0m)={r0['p_binomial']:.4g} -> p(1000m)={r1000['p_binomial']:.4g}")

    out = pd.DataFrame(rows)
    out_path = REPO_ROOT / "results" / "phase5_clustering_sensitivity.csv"
    out.to_csv(out_path, index=False)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    main()
