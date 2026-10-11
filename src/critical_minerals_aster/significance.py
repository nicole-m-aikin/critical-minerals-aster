"""
Statistical significance tests for ASTER anomaly-zone hit rates.

Null hypothesis (both tests): MRDS deposits are spatially uniform within the
site bounding box.  Under H₀ each deposit has probability p = (zone area /
bbox area) of landing inside a strong-anomaly zone.

Functions
---------
coverage_fraction   — zone area / bbox area in the zones' projected CRS
run_binomial        — exact one-sided p-value via scipy.stats.binomtest
run_permutation     — Monte Carlo p-value via rasterised grid sampling
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.stats import binomtest
from shapely.geometry import box
from shapely.ops import unary_union


def _bbox_in_crs(bbox_wgs84: tuple, crs) -> tuple[float, float, float, float]:
    """Return (minx, miny, maxx, maxy) for bbox_wgs84 reprojected to crs."""
    minlon, minlat, maxlon, maxlat = bbox_wgs84
    gdf = gpd.GeoDataFrame(geometry=[box(minlon, minlat, maxlon, maxlat)], crs="EPSG:4326")
    return gdf.to_crs(crs).geometry[0].bounds


def coverage_fraction(zones_gdf: gpd.GeoDataFrame, bbox_wgs84: tuple) -> float:
    """Fraction of projected bbox area covered by strong-anomaly zones (0–1)."""
    minx, miny, maxx, maxy = _bbox_in_crs(bbox_wgs84, zones_gdf.crs)
    bbox_area = (maxx - minx) * (maxy - miny)
    if bbox_area == 0:
        return 0.0
    union = unary_union(zones_gdf.geometry)
    return min(union.area / bbox_area, 1.0)


def run_binomial(n_hits: int, n_deposits: int, p: float) -> tuple[float, float]:
    """
    One-sided exact binomial test: H₀ hit_rate ≤ p.

    Returns
    -------
    p_value : float
    expected_hits : float  — n_deposits × p
    """
    if n_deposits == 0 or p <= 0:
        return 1.0, 0.0
    result = binomtest(n_hits, n_deposits, p, alternative="greater")
    return float(result.pvalue), n_deposits * p


def run_permutation(
    zones_gdf: gpd.GeoDataFrame,
    bbox_wgs84: tuple,
    n_deposits: int,
    n_hits: int,
    n_iter: int = 10_000,
    seed: int = 42,
    grid_res: int = 1_000,
    footprint: gpd.GeoDataFrame | None = None,
) -> float:
    """
    Monte Carlo spatial permutation p-value.

    Rasterises the zone union onto a ``grid_res × grid_res`` grid covering
    the bbox, then for each iteration samples ``n_deposits`` grid cells
    uniformly at random and counts how many are inside a zone.

    Returns P(random_hits ≥ n_hits) over ``n_iter`` iterations.

    Mathematically equivalent to placing n_deposits random points uniformly
    in the bbox and checking containment, but runs in milliseconds via numpy
    array indexing rather than repeated point-in-polygon tests.

    Parameters
    ----------
    footprint :
        Optional TIR valid-pixel footprint polygon (in ``zones_gdf.crs``).
        When supplied, random points are drawn only from grid cells inside
        the footprint rather than the full rectangular bbox — this matches
        the denominator used by ``coverage_fraction_footprint`` / the
        site-specific binomial null, so the Monte Carlo and analytical
        p-values are directly comparable rather than differing by a
        footprint-vs-bbox area mismatch (see docs/results.md, Phase 3).
        When omitted, falls back to the original bbox-rectangle behaviour.
    """
    if n_deposits == 0:
        return 1.0

    minx, miny, maxx, maxy = _bbox_in_crs(bbox_wgs84, zones_gdf.crs)

    import rasterio.features
    from rasterio.transform import from_bounds as rio_from_bounds

    union = unary_union(zones_gdf.geometry)
    transform = rio_from_bounds(minx, miny, maxx, maxy, grid_res, grid_res)
    mask = rasterio.features.rasterize(
        [(union, 1)],
        out_shape=(grid_res, grid_res),
        transform=transform,
        fill=0,
        dtype=np.uint8,
    ).ravel()

    rng = np.random.default_rng(seed)

    if footprint is not None and len(footprint):
        fp_geom = unary_union(footprint.geometry)
        valid = rasterio.features.rasterize(
            [(fp_geom, 1)],
            out_shape=(grid_res, grid_res),
            transform=transform,
            fill=0,
            dtype=np.uint8,
        ).ravel()
        cell_pool = np.flatnonzero(valid)
        if len(cell_pool) == 0:
            return 1.0
        # Sample n_deposits *positions within the pool* per iteration, then
        # map back to grid-cell ids — restricts every draw to the footprint.
        pool_positions = rng.integers(0, len(cell_pool), size=(n_iter, n_deposits))
        indices = cell_pool[pool_positions]
    else:
        n_cells = len(mask)
        indices = rng.integers(0, n_cells, size=(n_iter, n_deposits))

    null_hits = mask[indices].sum(axis=1)
    return float((null_hits >= n_hits).mean())


def run_permutation_stratified(
    zones_gdf: gpd.GeoDataFrame,
    bbox_wgs84: tuple,
    strata: dict,
    n_per_stratum: dict[str, int],
    n_hits: int,
    n_iter: int = 10_000,
    seed: int = 42,
    grid_res: int = 1_000,
    chunk: int = 1_000,
) -> float:
    """Monte Carlo p-value with random deposits drawn WITHIN their own stratum.

    Same grid machinery as ``run_permutation``, but each iteration places
    ``n_per_stratum[s]`` points uniformly in stratum ``s`` (geometries in
    ``zones_gdf.crs``), so stratum counts are preserved exactly — the
    simulation counterpart of the Poisson-binomial stratified null. Deposits
    whose stratum has no grid cell contribute no random hits (their
    analytical p_i is ~0 too). Returns P(random_hits ≥ n_hits).
    """
    n_total = int(sum(n_per_stratum.values()))
    if n_total == 0:
        return 1.0

    import rasterio.features
    from rasterio.transform import from_bounds as rio_from_bounds

    minx, miny, maxx, maxy = _bbox_in_crs(bbox_wgs84, zones_gdf.crs)
    transform = rio_from_bounds(minx, miny, maxx, maxy, grid_res, grid_res)
    shape = (grid_res, grid_res)
    mask = rasterio.features.rasterize(
        [(unary_union(zones_gdf.geometry), 1)], out_shape=shape, transform=transform,
        fill=0, dtype=np.uint8,
    ).ravel()

    pools = []
    for name, n_s in n_per_stratum.items():
        geom = strata.get(name)
        if n_s <= 0 or geom is None or geom.is_empty:
            continue
        cells = np.flatnonzero(rasterio.features.rasterize(
            [(geom, 1)], out_shape=shape, transform=transform, fill=0, dtype=np.uint8,
        ).ravel())
        if len(cells):
            pools.append((cells, int(n_s)))

    rng = np.random.default_rng(seed)
    exceed = 0
    for start in range(0, n_iter, chunk):
        size = min(chunk, n_iter - start)
        null_hits = np.zeros(size, dtype=np.int64)
        for cells, n_s in pools:
            pos = rng.integers(0, len(cells), size=(size, n_s))
            null_hits += mask[cells[pos]].sum(axis=1, dtype=np.int64)
        exceed += int((null_hits >= n_hits).sum())
    return exceed / n_iter


def _mwu_u_and_var(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Mann–Whitney U for ``x`` (ties count ½) and its tie-corrected null variance."""
    from scipy.stats import rankdata

    n1, n2 = len(x), len(y)
    ranks = rankdata(np.concatenate([x, y]))
    u = float(ranks[:n1].sum() - n1 * (n1 + 1) / 2)
    n = n1 + n2
    _, t = np.unique(np.concatenate([x, y]), return_counts=True)
    tie = float((t ** 3 - t).sum()) / (n * (n - 1)) if n > 1 else 0.0
    var = n1 * n2 / 12.0 * ((n + 1) - tie)
    return u, var


def stratified_mannwhitney(
    deposit_scores: np.ndarray,
    deposit_strata: np.ndarray,
    background_scores: np.ndarray,
    background_strata: np.ndarray,
    min_background: int = 3,
) -> dict:
    """One-sided stratified rank-sum test: deposits vs background of the SAME stratum.

    Within each stratum s, U_s compares that stratum's deposit scores only to
    that stratum's background pixels, so a detector that simply scores
    bedrock above alluvial cover cannot generate signal. Strata combine with
    deposit-count weights w_s = n_s / n:

        AUC = Σ w_s · U_s / (n_s · m_s)       (0.5 = no discrimination)
        Var = Σ w_s² · Var0(U_s) / (n_s · m_s)²   (tie-corrected)
        z = (AUC − 0.5) / √Var,  p = P(Z ≥ z)

    AUC reads as P(a deposit pixel outscores a random background pixel from
    the deposit's own stratum). With a single stratum this is the standard
    asymptotic Mann–Whitney test without continuity correction. Deposits in
    strata with fewer than ``min_background`` background pixels are dropped.
    """
    deposit_scores = np.asarray(deposit_scores, dtype=float)
    background_scores = np.asarray(background_scores, dtype=float)
    deposit_strata = np.asarray(deposit_strata)
    background_strata = np.asarray(background_strata)

    parts = []
    for s in np.unique(deposit_strata):
        x = deposit_scores[deposit_strata == s]
        y = background_scores[background_strata == s]
        if len(x) == 0 or len(y) < min_background:
            continue
        u, var = _mwu_u_and_var(x, y)
        parts.append((len(x), len(y), u, var))

    n_used = sum(p[0] for p in parts)
    if n_used == 0:
        return {"p": float("nan"), "auc": float("nan"), "z": float("nan"), "n_used": 0}
    auc, var_auc = 0.0, 0.0
    for n_s, m_s, u, var in parts:
        w = n_s / n_used
        auc += w * u / (n_s * m_s)
        var_auc += w ** 2 * var / (n_s * m_s) ** 2
    if var_auc <= 0:
        return {"p": float("nan"), "auc": auc, "z": float("nan"), "n_used": n_used}
    from scipy.stats import norm

    z = (auc - 0.5) / np.sqrt(var_auc)
    return {"p": float(norm.sf(z)), "auc": float(auc), "z": float(z), "n_used": int(n_used)}


def significant_site_sets(repo_root) -> tuple[tuple[str, ...], frozenset[str]]:
    """(FDR-significant site ids, clustering-robust subset) from the current results.

    FDR = ``sig_fdr_binomial`` in ``results/phase4_fdr_corrected_significance.csv``;
    robust = FDR sites with ``sig_05`` at every declustering radius in
    ``results/phase5_clustering_sensitivity.csv``.
    """
    from pathlib import Path

    results = Path(repo_root) / "results"
    fdr = pd.read_csv(results / "phase4_fdr_corrected_significance.csv")
    fdr_sites = tuple(sorted(fdr.loc[fdr["sig_fdr_binomial"].astype(bool), "site_id"]))
    p5 = pd.read_csv(results / "phase5_clustering_sensitivity.csv")
    piv = p5.pivot(index="site_id", columns="radius_m", values="sig_05").fillna(False).astype(bool)
    robust = frozenset(piv.index[piv.all(axis=1)]) & frozenset(fdr_sites)
    return fdr_sites, robust


# ---------------------------------------------------------------------------
# Wilson score confidence interval
# ---------------------------------------------------------------------------


def wilson_ci(k: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Wilson score confidence interval for a proportion k/n.

    Returns (lower, upper) as fractions in [0, 1].

    Why Wilson over the naive Wald interval (p ± z·se)?
    The Wald interval can produce negative lower bounds or bounds > 1 when
    the observed rate is extreme or n is small — both common in this dataset.
    Wilson is algebraically constrained to [0, 1] and has substantially
    better coverage properties for proportions near 0 or 1.

    Reference: Wilson (1927), J. Amer. Statist. Assoc.
    """
    if n == 0:
        return 0.0, 1.0
    import math
    # z-score for two-sided interval (e.g. 1.96 for 95%)
    z = {0.90: 1.645, 0.95: 1.960, 0.99: 2.576}.get(confidence, 1.960)
    p_hat = k / n
    denom = 1.0 + z ** 2 / n
    center = (p_hat + z ** 2 / (2 * n)) / denom
    margin = z * math.sqrt(p_hat * (1 - p_hat) / n + z ** 2 / (4 * n ** 2)) / denom
    return max(0.0, center - margin), min(1.0, center + margin)


# ---------------------------------------------------------------------------
# Footprint-aware coverage fraction (exposed for pipeline use)
# ---------------------------------------------------------------------------


def coverage_fraction_footprint(
    zones: gpd.GeoDataFrame,
    footprint: gpd.GeoDataFrame | None,
    bbox: tuple,
) -> float:
    """Zone area / denominator area, where denominator is the TIR footprint polygon
    (if available) or the rectangular bbox.

    Preferred over ``coverage_fraction`` when the TIR scene has a diagonal
    boundary (standard ASTER swath geometry) — the footprint polygon is more
    accurate than the bounding rectangle.
    """
    union = unary_union(zones.geometry)
    if footprint is not None and len(footprint):
        fp = footprint if footprint.crs == zones.crs else footprint.to_crs(zones.crs)
        denom = float(fp.iloc[0].geometry.area)
    else:
        minx, miny, maxx, maxy = _bbox_in_crs(bbox, zones.crs)
        denom = (maxx - minx) * (maxy - miny)
    if denom == 0:
        return 0.0
    return min(float(union.area / denom), 1.0)


# ---------------------------------------------------------------------------
# Annotate a site-summary DataFrame with uncertainty columns
# ---------------------------------------------------------------------------


def add_uncertainty_columns(
    summary: pd.DataFrame,
    zones: gpd.GeoDataFrame,
    footprint: gpd.GeoDataFrame | None,
    bbox: tuple,
    confidence: float = 0.95,
) -> pd.DataFrame:
    """Add Wilson CIs, null hit rate, and binomial p-value to every summary row.

    Works on the DataFrame produced by ``metrics.compute_site_summary()``.
    All per-row calculations use the same null probability (zone coverage
    fraction of the site footprint), because the null model is: a deposit
    is placed uniformly at random within the TIR footprint, so its probability
    of landing in a zone equals the fraction of the footprint covered by zones.

    New columns added
    -----------------
    null_hit_rate_pct : float
        Expected hit rate under H₀ = p_cover × 100.
    hit_rate_ci_low : float
        Wilson lower bound on observed hit rate (%), at *confidence* level.
    hit_rate_ci_high : float
        Wilson upper bound on observed hit rate (%), at *confidence* level.
    p_binomial : float
        One-sided exact binomial p-value: P(hits ≥ observed | n, p_cover).
        1.0 when n_deposits_bbox == 0 or p_cover == 0.
    """
    out = summary.copy()

    p_cover = coverage_fraction_footprint(zones, footprint, bbox)
    null_pct = round(p_cover * 100, 2)

    ci_lows, ci_highs, p_vals = [], [], []

    for _, row in out.iterrows():
        n = int(row["n_deposits_bbox"])
        k = int(row["n_deposits_in_zones"])
        lo, hi = wilson_ci(k, n, confidence)
        ci_lows.append(round(lo * 100, 2))
        ci_highs.append(round(hi * 100, 2))
        p_val, _ = run_binomial(k, n, p_cover)
        p_vals.append(round(p_val, 4))

    out["null_hit_rate_pct"] = null_pct
    out["hit_rate_ci_low"] = ci_lows
    out["hit_rate_ci_high"] = ci_highs
    out["p_binomial"] = p_vals
    return out


# ---------------------------------------------------------------------------
# MRDS location uncertainty — buffer sensitivity
# ---------------------------------------------------------------------------


def buffer_sensitivity(
    deposits: gpd.GeoDataFrame,
    zones: gpd.GeoDataFrame,
    radii_m: tuple[int, ...] = (0, 250, 500, 1000),
) -> pd.DataFrame:
    """Test hit-rate sensitivity to MRDS coordinate uncertainty.

    MRDS deposit coordinates are point estimates georeferenced from historical
    reports; accuracy varies from tens of metres to >1 km.  This function
    expands a circular buffer around each deposit point at increasing radii
    and asks: "does the zone polygon come within *radius* metres of this
    deposit?"

    At radius = 0 the result equals the standard point-in-polygon hit rate.
    If the hit rate increases sharply at 250–500 m, zones are spatially offset
    from deposit centers — either because ASTER detects the alteration *halo*
    (which peaks away from the ore body) or because MRDS coordinates are
    systematically imprecise.

    Parameters
    ----------
    deposits :
        GeoDataFrame of deposit points in a *projected* CRS (metres).
    zones :
        GeoDataFrame of strong-anomaly zone polygons in the same CRS.
    radii_m :
        Buffer radii in metres.  0 uses point-in-polygon; >0 buffers
        the point and tests polygon intersection.

    Returns
    -------
    DataFrame with columns:
        radius_m, n_deposits, n_hits, hit_rate_pct,
        hit_rate_ci_low, hit_rate_ci_high, delta_pct
    where *delta_pct* = hit_rate_pct − hit_rate at radius = 0.
    """
    if deposits.crs is None or deposits.crs.is_geographic:
        raise ValueError("deposits must be in a projected CRS (metres) for buffer_sensitivity")
    if zones.crs != deposits.crs:
        zones = zones.to_crs(deposits.crs)

    records = []
    baseline_hr = None

    for r in radii_m:
        if r == 0:
            geoms = deposits.geometry
            predicate = "within"
        else:
            geoms = deposits.geometry.buffer(r)
            predicate = "intersects"

        buffered = deposits.copy()
        buffered["geometry"] = geoms
        joined = gpd.sjoin(buffered, zones, how="left", predicate=predicate)
        n = len(deposits)
        hits = int(joined["index_right"].notna().groupby(level=0).any().sum())
        hr = hits / n * 100 if n else 0.0
        lo, hi = wilson_ci(hits, n)

        if baseline_hr is None:
            baseline_hr = hr
        delta = round(hr - baseline_hr, 1)

        records.append({
            "radius_m": r,
            "n_deposits": n,
            "n_hits": hits,
            "hit_rate_pct": round(hr, 1),
            "hit_rate_ci_low": round(lo * 100, 1),
            "hit_rate_ci_high": round(hi * 100, 1),
            "delta_pct": delta,
        })

    return pd.DataFrame(records)
