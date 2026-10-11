"""USGS USMIN mine-feature data: fetch, cache, and verification geometry helpers.

Motivation
----------
CIM (2018) Mineral Exploration Best Practice Guidelines §2.1.1.1 asks that
legacy data be both *validated* (internal consistency -- ``mrds.validate_mrds``)
and *verified* against an independent source.  For MRDS the natural
independent reference is USMIN: "Prospect- and Mine-Related Features from
U.S. Geological Survey 7.5- and 15-Minute Topographic Quadrangle Maps of the
United States" (Horton & San Juan 2016, ver. 10.0, doi:10.5066/F78W3CHG).
USMIN was digitised from historical topographic map symbols, so it shares no
lineage with MRDS report-based coordinates.

The same data answer a second question: ASTER TIR strong-anomaly zones could
be mapping mine waste (tailings, dumps, pits) rather than natural alteration.
USMIN polygons (tailings, mine dumps, open pits, disturbed surface) give a
disturbed-ground mask that can be removed from both the zones and the TIR
footprint before re-running the hit test.

Sources
-------
* USMIN WFS at ``USMIN_WFS_URL`` (MapServer; GML only, native NAD27).  Layers
  ``points-high`` (shafts, adits, prospect pits, quarries ...) and
  ``polygons``.  ``points-low`` is the same records at a different render
  scale and is not fetched.  Queries use WFS 1.0.0 because MapServer's 1.1.0
  endpoint silently returns nothing for a lat/lon-ordered bbox.
* OpenStreetMap via Overpass, for mine-related polygons where USMIN has none
  (or as a supplementary layer).  Overpass mirrors are frequently overloaded,
  so a failed fetch is reported as ``"unavailable"``, never as "no features".

All raw responses are cached under ``data/usmin/`` keyed by bbox, so reruns are
offline.  NAD27 -> WGS84 uses pyproj's grid-based NADCON transform when the
PROJ network/grids are available (~50-100 m shift in the western US, which
matters at a 250 m verification radius).
"""

from __future__ import annotations

import json
import math
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from shapely.geometry import LineString, Point, Polygon, box
from shapely.geometry.base import BaseGeometry
from shapely.ops import polygonize, unary_union
from shapely.strtree import STRtree

USMIN_WFS_URL = "https://mrdata.usgs.gov/services/usmin"
USMIN_CRS = "EPSG:4267"
USMIN_POINT_LAYER = "points-high"
USMIN_POLYGON_LAYER = "polygons"

OVERPASS_URLS = (
    "https://overpass-api.de/api/interpreter",
    "https://overpass.private.coffee/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
)

USER_AGENT = "critical-minerals-aster/0.1 (research; MRDS verification)"

#: Bbox area guard (deg^2).  Study-site rasters are ~0.2-0.55 deg^2.
MAX_BBOX_DEG2 = 1.0

#: 7.5-minute quadrangle size in degrees; USMIN is compiled per quad.
QUAD_DEG = 0.125

_GML = "{http://www.opengis.net/gml}"
_MS = "{http://mapserver.gis.umn.edu/mapserver}"

BBox = tuple[float, float, float, float]


# ---------------------------------------------------------------------------
# HTTP + cache
# ---------------------------------------------------------------------------


def _bbox_key(bbox: BBox) -> str:
    return "_".join(f"{v:.4f}" for v in bbox)


def _check_bbox(bbox: BBox) -> None:
    lon0, lat0, lon1, lat1 = bbox
    if not (lon0 < lon1 and lat0 < lat1):
        raise ValueError(f"Degenerate bbox {bbox}")
    if (lon1 - lon0) * (lat1 - lat0) > MAX_BBOX_DEG2:
        raise ValueError(f"bbox {bbox} exceeds MAX_BBOX_DEG2={MAX_BBOX_DEG2}; tile it first")


def _request(
    method: str,
    url: str,
    *,
    params: dict | None = None,
    data: dict | None = None,
    retries: int = 4,
    backoff: float = 2.0,
    timeout: float = 180.0,
) -> bytes:
    """HTTP request with exponential backoff on network errors, 429 and 5xx."""
    last: Exception | None = None
    for attempt in range(retries):
        try:
            resp = requests.request(
                method, url, params=params, data=data, timeout=timeout,
                headers={"User-Agent": USER_AGENT},
            )
            if resp.status_code == 429 or resp.status_code >= 500:
                raise requests.HTTPError(f"HTTP {resp.status_code}", response=resp)
            resp.raise_for_status()
            return resp.content
        except (requests.RequestException, OSError) as exc:
            last = exc
            if attempt < retries - 1:
                time.sleep(backoff * (2 ** attempt))
    raise RuntimeError(f"{method} {url} failed after {retries} attempts: {last!r}")


# ---------------------------------------------------------------------------
# USMIN WFS
# ---------------------------------------------------------------------------


def _parse_coords(text: str) -> list[tuple[float, float]]:
    out = []
    for tok in text.split():
        x, y = tok.split(",")[:2]
        out.append((float(x), float(y)))
    return out


def _ring(el: ET.Element) -> list[tuple[float, float]]:
    coords = el.find(f".//{_GML}coordinates")
    return _parse_coords(coords.text) if coords is not None and coords.text else []


def _gml2_geometry(el: ET.Element) -> BaseGeometry | None:
    """Convert a GML 2 geometry element (as emitted by MapServer WFS 1.0.0)."""
    tag = el.tag.replace(_GML, "")
    if tag == "Point":
        pts = _ring(el)
        return Point(pts[0]) if pts else None
    if tag == "LineString":
        pts = _ring(el)
        return LineString(pts) if len(pts) >= 2 else None
    if tag == "Polygon":
        outer_el = el.find(f"{_GML}outerBoundaryIs")
        if outer_el is None:
            return None
        outer = _ring(outer_el)
        holes = [_ring(h) for h in el.findall(f"{_GML}innerBoundaryIs")]
        poly = Polygon(outer, [h for h in holes if len(h) >= 4])
        return poly if poly.is_valid else poly.buffer(0)
    if tag in ("MultiPolygon", "MultiPoint", "MultiLineString", "MultiGeometry"):
        parts = []
        for member in el:
            for child in member:
                g = _gml2_geometry(child)
                if g is not None:
                    parts.append(g)
        return unary_union(parts) if parts else None
    return None


def parse_wfs_gml(content: bytes | str, crs: str = USMIN_CRS) -> gpd.GeoDataFrame:
    """Parse a MapServer WFS 1.0.0 GML 2 FeatureCollection into a GeoDataFrame.

    Raises ``RuntimeError`` on a ``ServiceExceptionReport`` so that a server
    error is never mistaken for an empty (no-mines) result.
    """
    root = ET.fromstring(content.encode() if isinstance(content, str) else content)
    if root.tag.endswith("ServiceExceptionReport"):
        msg = " ".join(t.strip() for t in root.itertext() if t.strip())
        raise RuntimeError(f"WFS service exception: {msg}")

    records: list[dict] = []
    geoms: list[BaseGeometry] = []
    for member in root.iter(f"{_GML}featureMember"):
        feat = next(iter(member), None)
        if feat is None:
            continue
        props: dict = {"fid": feat.get("fid")}
        geom = None
        for child in feat:
            name = child.tag.replace(_MS, "")
            if name == "geometry":
                inner = next(iter(child), None)
                geom = _gml2_geometry(inner) if inner is not None else None
            elif child.tag.startswith(_MS):
                props[name] = (child.text or "").strip()
        if geom is None or geom.is_empty:
            continue
        records.append(props)
        geoms.append(geom)
    return gpd.GeoDataFrame(records, geometry=geoms, crs=crs)


def fetch_usmin_layer(
    layer: str,
    bbox: BBox,
    cache_dir: Path,
    refresh: bool = False,
) -> gpd.GeoDataFrame:
    """Fetch one USMIN WFS layer for a WGS84 bbox (cached raw GML), in EPSG:4326.

    The bbox is interpreted by the server in the layer's native NAD27; the
    datum difference (<100 m) is negligible against a ~50 km site bbox.
    """
    _check_bbox(bbox)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"usmin_{layer}_{_bbox_key(bbox)}.gml"
    if refresh or not path.is_file():
        params = {
            "service": "WFS",
            "version": "1.0.0",
            "request": "GetFeature",
            "typeName": layer,
            "bbox": ",".join(f"{v:.6f}" for v in bbox),
        }
        content = _request("GET", USMIN_WFS_URL, params=params)
        parse_wfs_gml(content)  # raise before caching a service exception
        path.write_bytes(content)
        time.sleep(1.0)
    gdf = parse_wfs_gml(path.read_bytes())
    if gdf.empty:
        return gpd.GeoDataFrame(gdf, geometry=gdf.geometry, crs="EPSG:4326")
    return gdf.to_crs("EPSG:4326")


def load_usmin(
    bbox: BBox, cache_dir: Path, refresh: bool = False
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Return (points, polygons) USMIN GeoDataFrames in EPSG:4326 for *bbox*."""
    pts = fetch_usmin_layer(USMIN_POINT_LAYER, bbox, cache_dir, refresh)
    polys = fetch_usmin_layer(USMIN_POLYGON_LAYER, bbox, cache_dir, refresh)
    if not polys.empty:
        polys = polys[polys.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    pts["source"] = "usmin"
    polys["source"] = "usmin"
    return pts, polys


# ---------------------------------------------------------------------------
# OpenStreetMap / Overpass fallback
# ---------------------------------------------------------------------------

#: Overpass tag filters for mine-related features.  landuse=industrial is
#: only included when it carries an explicit mining tag.
OSM_FILTERS = (
    '["landuse"="quarry"]',
    '["disused:landuse"="quarry"]',
    '["abandoned:landuse"="quarry"]',
    '["man_made"~"^(spoil_heap|mineshaft|adit)$"]',
    '["historic"~"^(mine|mine_shaft|mine_adit)$"]',
    '["industrial"~"^(mine|mining)$"]',
    '["landuse"="brownfield"]',
    '["content"="tailings"]',
    '["basin"~"tailings"]',
    '["name"~"[Tt]ailings"]',
)


def overpass_query(bbox: BBox, timeout_s: int = 120) -> str:
    lon0, lat0, lon1, lat1 = bbox
    b = f"({lat0:.5f},{lon0:.5f},{lat1:.5f},{lon1:.5f})"
    body = "".join(f"nwr{f}{b};" for f in OSM_FILTERS)
    return f"[out:json][timeout:{timeout_s}];({body});out geom;"


def _way_geometry(coords: list[tuple[float, float]]) -> BaseGeometry | None:
    if len(coords) >= 4 and coords[0] == coords[-1]:
        poly = Polygon(coords)
        return poly if poly.is_valid else poly.buffer(0)
    if len(coords) >= 2:
        return LineString(coords)
    return None


def osm_elements_to_gdf(elements: Iterable[dict]) -> gpd.GeoDataFrame:
    """Convert Overpass ``out geom`` JSON elements to a GeoDataFrame (EPSG:4326).

    Nodes -> Points; closed ways -> Polygons (open ways -> LineStrings);
    multipolygon relations -> outer rings polygonised minus inner rings.
    """
    rows: list[dict] = []
    geoms: list[BaseGeometry] = []
    for el in elements:
        etype = el.get("type")
        tags = el.get("tags", {}) or {}
        geom: BaseGeometry | None = None
        if etype == "node" and "lon" in el:
            geom = Point(el["lon"], el["lat"])
        elif etype == "way" and el.get("geometry"):
            geom = _way_geometry([(p["lon"], p["lat"]) for p in el["geometry"]])
        elif etype == "relation":
            outers, inners = [], []
            for m in el.get("members", []):
                if m.get("type") != "way" or not m.get("geometry"):
                    continue
                line = LineString([(p["lon"], p["lat"]) for p in m["geometry"]])
                (inners if m.get("role") == "inner" else outers).append(line)
            if outers:
                outer = unary_union(list(polygonize(unary_union(outers))))
                if inners:
                    outer = outer.difference(unary_union(list(polygonize(unary_union(inners)))))
                geom = outer if not outer.is_empty else None
        if geom is None:
            continue
        ftr = next(
            (f"{k}={tags[k]}" for k in ("landuse", "man_made", "historic", "industrial",
                                        "disused:landuse", "abandoned:landuse", "basin", "content")
             if k in tags),
            "name~tailings",
        )
        rows.append({"osm_id": f"{etype}/{el.get('id')}", "ftr_type": ftr,
                     "ftr_name": tags.get("name", "")})
        geoms.append(geom)
    return gpd.GeoDataFrame(rows, geometry=geoms, crs="EPSG:4326")


def fetch_osm_mining(
    bbox: BBox, cache_dir: Path, refresh: bool = False
) -> tuple[gpd.GeoDataFrame | None, gpd.GeoDataFrame | None]:
    """Return (points, polygons) of OSM mine-related features, or (None, None)
    when every Overpass mirror fails.  ``None`` means *unknown*, not *none*."""
    _check_bbox(bbox)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"osm_mining_{_bbox_key(bbox)}.json"
    if refresh or not path.is_file():
        query = overpass_query(bbox)
        payload = None
        for url in OVERPASS_URLS:
            try:
                raw = _request("POST", url, data={"data": query}, retries=2, backoff=5.0, timeout=130.0)
                parsed = json.loads(raw)
            except (RuntimeError, ValueError):
                continue
            if "remark" in parsed and "error" in str(parsed["remark"]).lower():
                continue
            payload = parsed
            break
        if payload is None:
            return None, None
        path.write_text(json.dumps(payload))
        time.sleep(2.0)
    gdf = osm_elements_to_gdf(json.loads(path.read_text()).get("elements", []))
    gdf["source"] = "osm"
    pts = gdf[gdf.geom_type == "Point"].copy()
    polys = gdf[gdf.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    return pts, polys


# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------


def quad_cells(bbox: BBox, size: float = QUAD_DEG) -> list[Polygon]:
    """7.5-minute quadrangle cells (aligned to whole-degree multiples) covering bbox."""
    lon0, lat0, lon1, lat1 = bbox
    cells = []
    x = math.floor(lon0 / size) * size
    while x < lon1:
        y = math.floor(lat0 / size) * size
        while y < lat1:
            cells.append(box(x, y, x + size, y + size))
            y += size
        x += size
    return cells


def quad_coverage(
    footprint_wgs84: BaseGeometry, features_wgs84: gpd.GeoSeries
) -> tuple[int, int]:
    """(n 7.5' quads touching the footprint, n of those containing >=1 feature).

    A quad without features may be a mapped quad with no mines, so this is a
    density/extent indicator, not proof of missing compilation.
    """
    cells = [c for c in quad_cells(footprint_wgs84.bounds) if c.intersects(footprint_wgs84)]
    if features_wgs84 is None or len(features_wgs84) == 0:
        return len(cells), 0
    tree = STRtree(list(features_wgs84.geometry))
    n_hit = sum(1 for c in cells if len(tree.query(c, predicate="intersects")))
    return len(cells), n_hit


# ---------------------------------------------------------------------------
# Geometry helpers (projected CRS, metres)
# ---------------------------------------------------------------------------


def _geoms(obj) -> list[BaseGeometry]:
    if obj is None:
        return []
    if isinstance(obj, BaseGeometry):
        return [] if obj.is_empty else [obj]
    return [g for g in obj if g is not None and not g.is_empty]


def disturbed_mask(
    polygons: gpd.GeoSeries | Iterable[BaseGeometry] | None,
    points: gpd.GeoSeries | Iterable[BaseGeometry] | None = None,
    buffer_m: float = 0.0,
) -> BaseGeometry:
    """Union of disturbed-ground polygons plus *buffer_m* discs around points.

    Inputs must already be in a projected CRS in metres.  Points are ignored
    at ``buffer_m == 0`` (a zero-area point cannot mask anything).
    """
    parts = [g for g in _geoms(polygons) if g.geom_type in ("Polygon", "MultiPolygon")]
    if buffer_m > 0:
        parts += [g.buffer(buffer_m, quad_segs=8) for g in _geoms(points)]
    if not parts:
        return Polygon()
    return unary_union(parts)


def overlap_fraction(geom: BaseGeometry, mask: BaseGeometry) -> float:
    """Fraction of *geom*'s area covered by *mask* (0 when geom has no area)."""
    if geom.is_empty or geom.area == 0 or mask.is_empty:
        return 0.0
    return float(geom.intersection(mask).area / geom.area)


def subtract_mask(geom: BaseGeometry, mask: BaseGeometry) -> BaseGeometry:
    return geom if mask.is_empty else geom.difference(mask)


def nearest_feature_distance(
    targets: gpd.GeoSeries | Iterable[BaseGeometry],
    features: gpd.GeoSeries | Iterable[BaseGeometry],
) -> np.ndarray:
    """Distance (CRS units) from each target to the nearest feature.

    Polygon features count as distance 0 for targets inside them.  Returns
    NaN for every target when *features* is empty (no reference available).
    """
    tgt = list(targets.geometry if isinstance(targets, gpd.GeoDataFrame) else targets)
    feats = _geoms(features.geometry if isinstance(features, gpd.GeoDataFrame) else features)
    if not tgt:
        return np.array([], dtype=float)
    if not feats:
        return np.full(len(tgt), np.nan)
    tree = STRtree(feats)
    idx = tree.query_nearest(tgt, return_distance=False, all_matches=False)
    out = np.full(len(tgt), np.nan)
    for ti, fi in zip(idx[0], idx[1]):
        out[ti] = tgt[ti].distance(feats[fi])
    return out


def points_in_mask(points: gpd.GeoSeries, mask: BaseGeometry) -> np.ndarray:
    """Boolean array: point intersects mask (boundary counts as inside)."""
    if mask.is_empty or len(points) == 0:
        return np.zeros(len(points), dtype=bool)
    return np.asarray(points.intersects(mask), dtype=bool)


def feature_type_counts(gdf: gpd.GeoDataFrame | None) -> str:
    """Compact 'type:n;type:n' summary for provenance columns."""
    if gdf is None or gdf.empty or "ftr_type" not in gdf.columns:
        return ""
    vc = gdf["ftr_type"].replace("", "unspecified").value_counts()
    return ";".join(f"{k}:{v}" for k, v in vc.items())


def concat_gdfs(frames: Iterable[gpd.GeoDataFrame | None], crs="EPSG:4326") -> gpd.GeoDataFrame:
    frames = [f for f in frames if f is not None and not f.empty]
    if not frames:
        return gpd.GeoDataFrame(geometry=gpd.GeoSeries([], crs=crs), crs=crs)
    return gpd.GeoDataFrame(pd.concat([f.to_crs(crs) for f in frames], ignore_index=True), crs=crs)


# ---------------------------------------------------------------------------
# Per-site verification inputs (cache-only for zones / MRDS)
# ---------------------------------------------------------------------------


@dataclass
class SiteInputs:
    """Everything the verification scripts need for one site, in ``crs``.

    ``deposits`` holds the critical-mineral MRDS records inside the TIR
    footprint with an ``inside_zone`` flag -- the same set and hit definition
    as ``n_crit`` / ``hits_crit`` in the Phase 2/4 significance tables, and
    ``p0`` is ``significance.coverage_fraction_footprint`` (zone union area /
    first footprint feature area).
    """

    site_id: str
    site_name: str
    crs: object
    raster_bbox: BBox
    zones_union: BaseGeometry
    footprint: BaseGeometry
    footprint_wgs84: BaseGeometry
    deposits: gpd.GeoDataFrame
    p0: float
    usmin_points: gpd.GeoDataFrame | None = None
    usmin_polygons: gpd.GeoDataFrame | None = None
    osm_points: gpd.GeoDataFrame | None = None
    osm_polygons: gpd.GeoDataFrame | None = None
    usmin_status: str = "not_fetched"
    osm_status: str = "not_fetched"


def load_site_inputs(
    site_id: str,
    repo_root: Path,
    cache_dir: Path,
    fetch_osm: bool = True,
    refresh: bool = False,
) -> SiteInputs:
    """Load cached zones / provenance / MRDS for a site and fetch USMIN (+OSM).

    Mirrors the loading in ``scripts/regenerate_site_summaries.py``.  USMIN
    and OSM layers are clipped to the TIR footprint; ``usmin_status`` /
    ``osm_status`` are ``"ok"``, ``"empty"`` (service answered, zero features
    in footprint) or ``"unavailable"`` (fetch failed -- coverage unknown).
    """
    from critical_minerals_aster.config import load_site_by_id
    from critical_minerals_aster.metrics import read_mrds_national
    from critical_minerals_aster.mrds import (
        filter_mrds_bbox,
        is_critical_mineral,
        mrds_to_points_gdf,
        reclassify_mrds_earth_mri,
        spatial_join_deposits_zones,
    )
    from critical_minerals_aster.paths import site_paths_for

    repo_root = Path(repo_root)
    site = load_site_by_id(site_id, repo_root / "sites")
    paths = site_paths_for(site, repo_root)
    zones = gpd.read_file(paths.strong_zones_geojson)
    with open(repo_root / "results" / f"{site_id}_provenance.json") as f:
        prov = json.load(f)
    raster_bbox = tuple(prov["raster_bbox_wgs84"]) if prov.get("raster_bbox_wgs84") else site.bbox_wgs84

    fp_gdf = None
    if prov.get("tir_footprint_wgs84") is not None:
        fp_gdf = gpd.GeoDataFrame.from_features(
            prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
        ).to_crs(zones.crs)
    if fp_gdf is None or fp_gdf.empty:
        fp_gdf = gpd.GeoDataFrame(geometry=[box(*raster_bbox)], crs="EPSG:4326").to_crs(zones.crs)
    footprint = fp_gdf.iloc[0].geometry
    footprint_wgs84 = gpd.GeoSeries([footprint], crs=zones.crs).to_crs("EPSG:4326").iloc[0]

    zones_union = unary_union(zones.geometry)
    p0 = min(zones_union.area / footprint.area, 1.0) if footprint.area else 0.0

    mrds = read_mrds_national(paths)
    deposits = mrds_to_points_gdf(filter_mrds_bbox(mrds, raster_bbox), zones.crs)
    deposits = gpd.clip(deposits, fp_gdf)
    deposits = reclassify_mrds_earth_mri(deposits)
    deposits = deposits[deposits["earth_mri_category"].apply(is_critical_mineral)].copy()
    _, hits, _ = spatial_join_deposits_zones(deposits, zones)
    deposits["inside_zone"] = deposits.index.isin(set(hits.index.unique()))

    inputs = SiteInputs(
        site_id=site_id, site_name=site.name, crs=zones.crs, raster_bbox=raster_bbox,
        zones_union=zones_union, footprint=footprint, footprint_wgs84=footprint_wgs84,
        deposits=deposits, p0=p0,
    )

    def _clip(gdf: gpd.GeoDataFrame | None) -> gpd.GeoDataFrame | None:
        if gdf is None:
            return None
        if gdf.empty:
            return gdf.set_crs(zones.crs, allow_override=True)
        proj = gdf.to_crs(zones.crs)
        return proj[proj.intersects(footprint)].copy()

    try:
        pts, polys = load_usmin(raster_bbox, cache_dir, refresh)
        inputs.usmin_points, inputs.usmin_polygons = _clip(pts), _clip(polys)
        n = len(inputs.usmin_points) + len(inputs.usmin_polygons)
        inputs.usmin_status = "ok" if n else "empty"
    except RuntimeError:
        inputs.usmin_status = "unavailable"

    if fetch_osm:
        opts, opolys = fetch_osm_mining(raster_bbox, cache_dir, refresh)
        if opts is None:
            inputs.osm_status = "unavailable"
        else:
            inputs.osm_points, inputs.osm_polygons = _clip(opts), _clip(opolys)
            n = len(inputs.osm_points) + len(inputs.osm_polygons)
            inputs.osm_status = "ok" if n else "empty"
    return inputs
