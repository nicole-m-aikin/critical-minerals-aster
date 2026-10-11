"""
Lithology-conditioned null models for anomaly-zone / deposit overlap.

Addresses the reviewer risk that strong-anomaly zones and MRDS deposits both
simply track host lithology: a site-level geometric null (p0 = zone area /
footprint area) cannot distinguish "ASTER sees alteration where deposits are"
from "ASTER sees limestone, and deposits are in limestone".

Data source: USGS State Geologic Map Compilation (SGMC; Horton, San Juan &
Stoeser 2017, USGS DS 1052, doi:10.5066/F7WH2N65), polygon layer 3
``SGMC_Geology`` of the ScienceBase-hosted ArcGIS FeatureServer. Source map
scales range 1:50,000-1:1,000,000 (most western states ~1:500,000).

Functions
---------
fetch_sgmc_units          — paged, bbox-bounded FeatureServer query (retries/backoff)
load_sgmc_units           — cache-first loader → ``data/lithology/{site_id}_sgmc_units.geojson``
generalize_lithology      — SGMC GENERALIZED_LITH → coarse lithology class
build_strata              — disjoint per-class polygons clipped to the footprint
stratum_table             — per-stratum area, zone area, p0_L, n_L, hits_L
poisson_binomial_sf       — exact P(X ≥ k) for a sum of per-stratum binomials
stratified_test           — lithology-stratified expectation, enrichment, p
host_restricted_test      — binomial test restricted to deposit-hosting strata

Bedrock-vs-cover stratification (the PRIMARY null of the analysis chain)
-----------------------------------------------------------------------
cover_strata_from_units   — {bedrock, cover, unmapped} polygons within a footprint
cover_strata              — same, loading (and fetching if missing) the SGMC cache
assign_cover_stratum      — stratum label per deposit point
stratum_zone_fractions    — per-stratum area, zone area, p = zone ∩ stratum / stratum
stratified_null           — E = Σ p_i, effective p0 = E / n, Poisson-binomial tail
rasterize_cover_strata    — stratum codes on a raster grid (for pixel-level tests)
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.stats import binom
from shapely.geometry import Polygon, box
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union
from shapely.validation import make_valid

from critical_minerals_aster.significance import run_binomial

SGMC_GEOLOGY_URL = (
    "https://services.arcgis.com/v01gqwM5QqNysAAi/arcgis/rest/services"
    "/SB_5888bf4fe4b05ccb964bab9d_USGS_SGMC_feature/FeatureServer/3/query"
)
SGMC_FIELDS = (
    "STATE", "SGMC_LABEL", "UNIT_LINK", "UNIT_NAME", "AGE_MIN", "AGE_MAX",
    "MAJOR1", "MAJOR2", "MINOR1", "GENERALIZED_LITH",
)
SGMC_SOURCE = (
    "USGS SGMC (Horton et al. 2017, doi:10.5066/F7WH2N65), "
    "FeatureServer layer 3 SGMC_Geology"
)

UNMAPPED = "unmapped"

# SGMC GENERALIZED_LITH (33 values nationally) → 8 coarse classes + unmapped.
# Coarse enough that strata are not tiny at ~1:500k, fine enough to separate
# the carbonate hosts of skarn/CRD from clastic hosts of roll-front uranium.
LITH_CLASS_MAP: dict[str, str] = {
    "Igneous, intrusive": "intrusive",
    "Metamorphic, intrusive": "intrusive",
    "Igneous, volcanic": "volcanic",
    "Metamorphic, volcanic": "volcanic",
    "Igneous, undifferentiated": "igneous_mixed",
    "Igneous and Sedimentary, undifferentiated": "igneous_mixed",
    "Sedimentary, carbonate": "carbonate_chemical",
    "Metamorphic, carbonate": "carbonate_chemical",
    "Sedimentary, chemical": "carbonate_chemical",
    "Sedimentary, evaporite": "carbonate_chemical",
    "Sedimentary, iron formation, undifferentiated": "carbonate_chemical",
    "Sedimentary, clastic": "clastic",
    "Metamorphic, sedimentary clastic": "clastic",
    "Sedimentary, undifferentiated": "sedimentary_mixed",
    "Metamorphic, sedimentary": "sedimentary_mixed",
    "Metamorphic and Sedimentary, undifferentiated": "sedimentary_mixed",
    "Metamorphic, amphibolite": "metamorphic",
    "Metamorphic, gneiss": "metamorphic",
    "Metamorphic, granulite": "metamorphic",
    "Metamorphic, igneous": "metamorphic",
    "Metamorphic, other": "metamorphic",
    "Metamorphic, schist": "metamorphic",
    "Metamorphic, serpentinite": "metamorphic",
    "Metamorphic, undifferentiated": "metamorphic",
    "Igneous and Metamorphic, undifferentiated": "metamorphic",
    "Tectonite, undifferentiated": "metamorphic",
    "Melange": "metamorphic",
    "Unconsolidated, undifferentiated": "unconsolidated",
    "Unconsolidated and Sedimentary, undifferentiated": "unconsolidated",
    "Water": "other",
    "Ice": "other",
    "Dam": "other",
    "Unknown": "other",
}

def generalize_lithology(generalized_lith: str | None) -> str:
    """Map an SGMC GENERALIZED_LITH value to a coarse lithology class.

    ``sedimentary_mixed`` is deliberately not split by MAJOR1: in the 1:1M
    Arizona map, carbonate-dominated Paleozoic sections (e.g. Bisbee's
    Escabrosa/Naco hosts) carry MAJOR1 = "Sandstone" with limestone only as a
    minor component, so a MAJOR1 split would merge them into Cenozoic clastic
    basin fill.
    """
    return LITH_CLASS_MAP.get((generalized_lith or "").strip(), "other")


def add_lith_class(units: gpd.GeoDataFrame, col: str = "lith_class") -> gpd.GeoDataFrame:
    """Return a copy of SGMC units with a coarse ``lith_class`` column."""
    out = units.copy()
    out[col] = [generalize_lithology(g) for g in out["GENERALIZED_LITH"]]
    return out


# ---------------------------------------------------------------------------
# SGMC fetch / cache
# ---------------------------------------------------------------------------


def _http_get_json(url: str, timeout: int, max_retries: int) -> dict:
    last_exc: Exception | None = None
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "critical-minerals-aster/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            if "error" in data:
                raise RuntimeError(f"ArcGIS error: {data['error']}")
            return data
        except Exception as exc:  # network, JSON, or service error — retry
            last_exc = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"SGMC query failed after {max_retries} attempts: {last_exc!r}")


def fetch_sgmc_units(
    bbox_wgs84: tuple[float, float, float, float],
    *,
    page_size: int = 2000,
    timeout: int = 60,
    max_retries: int = 4,
    pause_s: float = 0.25,
) -> dict:
    """Fetch SGMC geologic-unit polygons intersecting a WGS84 bbox as GeoJSON.

    Pages through the FeatureServer (maxRecordCount = 2000) with
    ``resultOffset`` ordered by OBJECTID until ``exceededTransferLimit`` clears.
    """
    west, south, east, north = bbox_wgs84
    features: list[dict] = []
    offset = 0
    while True:
        params = {
            "where": "1=1",
            "geometry": f"{west},{south},{east},{north}",
            "geometryType": "esriGeometryEnvelope",
            "inSR": "4326",
            "outSR": "4326",
            "spatialRel": "esriSpatialRelIntersects",
            "outFields": ",".join(SGMC_FIELDS),
            "orderByFields": "OBJECTID",
            "resultOffset": offset,
            "resultRecordCount": page_size,
            "f": "geojson",
        }
        data = _http_get_json(SGMC_GEOLOGY_URL + "?" + urllib.parse.urlencode(params), timeout, max_retries)
        page = data.get("features", [])
        features.extend(page)
        exceeded = data.get("exceededTransferLimit") or (data.get("properties") or {}).get(
            "exceededTransferLimit"
        )
        if not page or not exceeded:
            break
        offset += len(page)
        time.sleep(pause_s)
    return {"type": "FeatureCollection", "features": features}


def load_sgmc_units(
    site_id: str,
    bbox_wgs84: tuple[float, float, float, float],
    repo_root: Path,
    *,
    refresh: bool = False,
) -> gpd.GeoDataFrame:
    """Cache-first SGMC unit loader (EPSG:4326), clipped to ``bbox_wgs84``.

    Polygons are clipped to the bbox before caching: SGMC units such as
    valley-fill alluvium can span hundreds of km, and the unclipped
    geometries would dominate the cache size.
    """
    cache = repo_root / "data" / "lithology" / f"{site_id}_sgmc_units.geojson"
    if cache.is_file() and not refresh:
        gdf = gpd.read_file(cache)
        return gdf.set_crs("EPSG:4326") if gdf.crs is None else gdf

    fc = fetch_sgmc_units(bbox_wgs84)
    if fc["features"]:
        gdf = gpd.GeoDataFrame.from_features(fc["features"], crs="EPSG:4326")
        gdf["geometry"] = gdf.geometry.apply(lambda g: make_valid(g) if g is not None else g)
        gdf = gpd.clip(gdf, box(*bbox_wgs84))
        gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()]
        gdf = gdf[gdf.geom_type.isin(["Polygon", "MultiPolygon", "GeometryCollection"])]
    else:
        gdf = gpd.GeoDataFrame(columns=[*SGMC_FIELDS, "geometry"], geometry="geometry", crs="EPSG:4326")
    cache.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(cache, driver="GeoJSON")
    return gdf


# ---------------------------------------------------------------------------
# Pure statistics on strata
# ---------------------------------------------------------------------------


def _polygonal(geom: BaseGeometry) -> BaseGeometry:
    """Drop point/line debris that make_valid / intersection can leave behind."""
    if geom.is_empty:
        return geom
    if geom.geom_type in ("Polygon", "MultiPolygon"):
        return geom
    if hasattr(geom, "geoms"):
        parts = [g for g in geom.geoms if g.geom_type in ("Polygon", "MultiPolygon")]
        return unary_union(parts) if parts else Polygon()
    return Polygon()


def build_strata(
    units: gpd.GeoDataFrame,
    footprint: BaseGeometry,
    class_col: str = "lith_class",
) -> gpd.GeoDataFrame:
    """Disjoint per-class polygons covering ``footprint`` exactly.

    SGMC is not reconciled across state lines, so units may overlap; classes
    are made disjoint by subtracting earlier classes (alphabetical order,
    deterministic). Footprint area not covered by any SGMC unit becomes the
    ``unmapped`` stratum so stratum areas always sum to the footprint area.
    ``units`` must already be in the same projected CRS as ``footprint``.
    """
    records: list[dict] = []
    covered = None
    if len(units):
        for cls, grp in sorted(units.groupby(class_col), key=lambda kv: str(kv[0])):
            parts = [_polygonal(make_valid(g)) for g in grp.geometry.values if g is not None]
            geom = _polygonal(unary_union(parts).intersection(footprint))
            if covered is not None and not geom.is_empty:
                geom = _polygonal(geom.difference(covered))
            if geom.is_empty or geom.area <= 0:
                continue
            records.append({class_col: str(cls), "geometry": geom})
            covered = geom if covered is None else covered.union(geom)
    rest = footprint if covered is None else _polygonal(footprint.difference(covered))
    if not rest.is_empty and rest.area > 0:
        records.append({class_col: UNMAPPED, "geometry": rest})
    return gpd.GeoDataFrame(records, geometry="geometry", crs=units.crs if len(units) else None)


def assign_strata(points: gpd.GeoDataFrame, strata: gpd.GeoDataFrame, class_col: str = "lith_class") -> pd.Series:
    """Stratum label per point; points on no stratum (numerical edge) → ``unmapped``."""
    labels = pd.Series(UNMAPPED, index=points.index, dtype=object)
    for cls, geom in zip(strata[class_col], strata.geometry):
        mask = points.geometry.intersects(geom) & (labels == UNMAPPED)
        labels[mask] = cls
    return labels


def stratum_table(
    strata: gpd.GeoDataFrame,
    zone_union: BaseGeometry,
    deposit_strata: pd.Series,
    deposit_hits: pd.Series,
    class_col: str = "lith_class",
) -> pd.DataFrame:
    """Per-stratum geometry and deposit counts.

    Columns: ``class_col``, area, zone_area, p0 (= zone_area / area),
    n (deposits), hits (deposits inside a zone). Areas are in CRS units².
    Strata with deposits but no area (points snapped to ``unmapped`` with no
    unmapped polygon) get area 0 and p0 0.
    """
    rows = []
    for cls, geom in zip(strata[class_col], strata.geometry):
        area = float(geom.area)
        zone_area = float(_polygonal(geom.intersection(zone_union)).area) if area > 0 else 0.0
        rows.append({class_col: cls, "area": area, "zone_area": zone_area})
    tab = pd.DataFrame(rows, columns=[class_col, "area", "zone_area"])
    counts = pd.DataFrame({"s": deposit_strata.values, "h": deposit_hits.astype(bool).values})
    n_by = counts.groupby("s").size()
    h_by = counts[counts["h"]].groupby("s").size()
    missing = sorted(set(n_by.index) - set(tab[class_col]))
    if missing:
        tab = pd.concat(
            [tab, pd.DataFrame({class_col: missing, "area": 0.0, "zone_area": 0.0})], ignore_index=True
        )
    tab["n"] = tab[class_col].map(n_by).fillna(0).astype(int)
    tab["hits"] = tab[class_col].map(h_by).fillna(0).astype(int)
    tab["p0"] = np.where(tab["area"] > 0, tab["zone_area"] / tab["area"].where(tab["area"] > 0, 1.0), 0.0)
    return tab


def poisson_binomial_sf(n_per: np.ndarray, p_per: np.ndarray, k: int) -> float:
    """Exact P(X ≥ k) for X = Σ_L Binomial(n_L, p_L), independent strata.

    Convolves the per-stratum binomial pmfs; the upper tail is summed
    directly (not 1 − cdf) so very small p-values keep their precision.
    """
    pmf = np.array([1.0])
    for n_l, p_l in zip(np.asarray(n_per, dtype=int), np.asarray(p_per, dtype=float)):
        if n_l <= 0:
            continue
        pmf = np.convolve(pmf, binom.pmf(np.arange(n_l + 1), n_l, min(max(p_l, 0.0), 1.0)))
    if k <= 0:
        return 1.0
    if k >= len(pmf):
        return 0.0
    return float(min(1.0, pmf[k:].sum()))


def stratified_test(tab: pd.DataFrame) -> dict:
    """Lithology-stratified null: expected = Σ n_L p0_L; exact Poisson-binomial p."""
    n = int(tab["n"].sum())
    hits = int(tab["hits"].sum())
    expected = float((tab["n"] * tab["p0"]).sum())
    p = poisson_binomial_sf(tab["n"].values, tab["p0"].values, hits) if n else 1.0
    return {
        "n": n,
        "hits": hits,
        "expected": expected,
        "enrichment": hits / expected if expected > 0 else float("nan"),
        "p": p,
    }


def host_restricted_test(tab: pd.DataFrame, class_col: str = "lith_class") -> dict:
    """Binomial test with footprint, zones and deposits restricted to host strata.

    Host strata = classes containing ≥ 1 deposit. Every deposit is in a host
    stratum by construction, so n and hits are unchanged; only p0 moves to
    zone area ∩ host / host area.
    """
    host = tab[tab["n"] > 0]
    n = int(host["n"].sum())
    hits = int(host["hits"].sum())
    area = float(host["area"].sum())
    p0 = float(host["zone_area"].sum() / area) if area > 0 else 0.0
    p, expected = run_binomial(hits, n, p0)
    return {
        "n": n,
        "hits": hits,
        "p0": p0,
        "expected": expected,
        "enrichment": (hits / n) / p0 if n and p0 > 0 else float("nan"),
        "p": p,
        "host_classes": sorted(host[class_col].astype(str).tolist()),
        "host_area_frac": area / float(tab["area"].sum()) if tab["area"].sum() > 0 else float("nan"),
    }


# ---------------------------------------------------------------------------
# Bedrock / cover / unmapped strata — primary null of the analysis chain
# ---------------------------------------------------------------------------
#
# Deposits sit in bedrock; ASTER TIR zones move between outcrop and alluvial
# valley fill with acquisition temperature. A geometric null (p0 = zone area /
# footprint area) therefore rewards any detector that preferentially flags
# outcrop. Conditioning each deposit's null probability on its own stratum
# removes that bias while keeping the null as coarse as possible.

BEDROCK = "bedrock"
COVER = "cover"
COVER_STRATA = (BEDROCK, COVER, UNMAPPED)
COVER_CLASSES = frozenset({"unconsolidated"})
STRATUM_CODES = {BEDROCK: 1, COVER: 2, UNMAPPED: 3}

_REPO_ROOT = Path(__file__).resolve().parents[2]


def cover_strata_from_units(units: gpd.GeoDataFrame, footprint: BaseGeometry) -> dict[str, BaseGeometry]:
    """Partition ``footprint`` into disjoint bedrock / cover / unmapped polygons.

    cover    = footprint ∩ SGMC units whose lith class is in ``COVER_CLASSES``
    unmapped = footprint not covered by any SGMC unit (e.g. Mexico, map gaps)
    bedrock  = the rest (every other SGMC class, incl. water/unknown "other")

    Where SGMC units overlap (state-line mismatches), cover wins. ``units``
    must be in the footprint's projected CRS; ``lith_class`` is derived from
    ``GENERALIZED_LITH`` when absent. Areas always sum to ``footprint.area``.
    """
    empty = Polygon()
    if len(units) == 0:
        return {BEDROCK: empty, COVER: empty, UNMAPPED: footprint}
    if "lith_class" not in units.columns:
        units = add_lith_class(units)
    geoms = [_polygonal(make_valid(g)) for g in units.geometry.values if g is not None]
    is_cover = [c in COVER_CLASSES for g, c in zip(units.geometry.values, units["lith_class"]) if g is not None]
    mapped = _polygonal(unary_union(geoms).intersection(footprint))
    cover_parts = [g for g, c in zip(geoms, is_cover) if c]
    cover = _polygonal(unary_union(cover_parts).intersection(footprint)) if cover_parts else empty
    bedrock = _polygonal(mapped.difference(cover)) if not cover.is_empty else mapped
    unmapped = _polygonal(footprint.difference(mapped)) if not mapped.is_empty else footprint
    return {BEDROCK: bedrock, COVER: cover, UNMAPPED: unmapped}


def cover_strata(
    site_id: str,
    footprint_geom: BaseGeometry,
    crs,
    *,
    bbox_wgs84: tuple[float, float, float, float] | None = None,
    repo_root: Path | None = None,
    refresh: bool = False,
) -> dict[str, BaseGeometry]:
    """Bedrock / cover / unmapped strata for a site's footprint (geometries in ``crs``).

    Loads the SGMC cache ``data/lithology/{site_id}_sgmc_units.geojson``,
    fetching it first if missing. ``bbox_wgs84`` (default: the footprint's
    WGS84 bounds) bounds the fetch only.
    """
    root = repo_root or _REPO_ROOT
    if bbox_wgs84 is None:
        bbox_wgs84 = tuple(gpd.GeoSeries([footprint_geom], crs=crs).to_crs("EPSG:4326").total_bounds)
    units = load_sgmc_units(site_id, bbox_wgs84, root, refresh=refresh)
    if len(units):
        units = units.to_crs(crs)
        units = units[units.geometry.intersects(footprint_geom)]
    return cover_strata_from_units(units, footprint_geom)


def assign_cover_stratum(points: gpd.GeoDataFrame, strata: dict[str, BaseGeometry]) -> pd.Series:
    """Stratum label per point.

    Points on no stratum polygon (boundary round-off, or points just outside
    the footprint) take the nearest non-empty stratum, so every deposit gets
    the null probability of the ground it actually sits on.
    """
    labels = pd.Series(None, index=points.index, dtype=object)
    for name in COVER_STRATA:
        geom = strata.get(name)
        if geom is None or geom.is_empty:
            continue
        mask = labels.isna() & points.geometry.intersects(geom)
        labels[mask] = name
    missing = labels.isna()
    if missing.any():
        live = [(n, strata[n]) for n in COVER_STRATA if strata.get(n) is not None and not strata[n].is_empty]
        if not live:
            labels[missing] = UNMAPPED
        else:
            for idx in labels.index[missing]:
                pt = points.geometry.loc[idx]
                labels.loc[idx] = min(live, key=lambda kv: kv[1].distance(pt))[0]
    return labels


def stratum_zone_fractions(strata: dict[str, BaseGeometry], zone_union: BaseGeometry) -> pd.DataFrame:
    """Per stratum: area, zone_area (zone ∩ stratum) and p = zone_area / area (0 if area 0)."""
    rows = []
    for name in COVER_STRATA:
        geom = strata.get(name)
        area = float(geom.area) if geom is not None and not geom.is_empty else 0.0
        zone_area = float(_polygonal(geom.intersection(zone_union)).area) if area > 0 else 0.0
        rows.append({"stratum": name, "area": area, "zone_area": zone_area,
                     "p": min(zone_area / area, 1.0) if area > 0 else 0.0})
    return pd.DataFrame(rows).set_index("stratum")


def poisson_binomial_sf_trials(p_i: np.ndarray, k: int) -> float:
    """Exact P(X ≥ k) for X = Σ_i Bernoulli(p_i), grouping trials with equal p."""
    p_i = np.asarray(p_i, dtype=float)
    if len(p_i) == 0:
        return 1.0 if k <= 0 else 0.0
    vals, counts = np.unique(p_i, return_counts=True)
    return poisson_binomial_sf(counts, vals, k)


def stratified_null(p_i: np.ndarray, hits: int) -> dict:
    """Stratified null for n deposits with per-deposit null probabilities ``p_i``.

    expected = Σ p_i; p0_eff = expected / n (the hit rate the null predicts,
    so hits / expected keeps the meaning of "enrichment"); p = exact
    Poisson-binomial upper tail P(X ≥ hits).
    """
    p_i = np.asarray(p_i, dtype=float)
    n = len(p_i)
    expected = float(p_i.sum())
    return {
        "n": n,
        "hits": int(hits),
        "expected": expected,
        "p0_eff": expected / n if n else 0.0,
        "enrichment": hits / expected if expected > 0 else float("nan"),
        "p": poisson_binomial_sf_trials(p_i, int(hits)) if n else 1.0,
    }


def rasterize_cover_strata(
    strata: dict[str, BaseGeometry],
    transform,
    shape: tuple[int, int],
) -> np.ndarray:
    """``uint8`` grid of ``STRATUM_CODES`` (0 = outside every stratum polygon)."""
    import rasterio.features

    shapes = [(strata[n], STRATUM_CODES[n]) for n in COVER_STRATA
              if strata.get(n) is not None and not strata[n].is_empty]
    if not shapes:
        return np.zeros(shape, dtype=np.uint8)
    return rasterio.features.rasterize(shapes, out_shape=shape, transform=transform,
                                       fill=0, dtype=np.uint8)
