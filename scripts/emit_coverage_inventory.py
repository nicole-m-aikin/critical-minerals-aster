#!/usr/bin/env python3
"""EMIT (then EnMAP/PRISMA schema) coverage inventory over the 109 site bboxes.

Shares the gemx_coverage schema columns:
  site_id, collection, sensor, coverage_frac, n_granules, status

EMIT L2A Reflectance is queried via NASA CMR (no credentials for search).
EnMAP/PRISMA rows are emitted as ``status=not_queried`` placeholders until
portal credentials are configured.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
from shapely.geometry import box
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.config import list_site_ids, load_site_by_id
from critical_minerals_aster.hyperspectral import EQUAL_AREA_CRS

CMR_GRANULE_URL = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"

# EMIT L2A Reflectance — LP DAAC (V002 preferred; V001 retained for continuity)
EMIT_COLLECTIONS = (
    {
        "key": "emit_l2a_rfl",
        "concept_id": "C4079844428-LPCLOUD",
        "sensor": "EMIT",
        "gsd_m": 60,
        "note": "EMIT L2A Surface Reflectance V002 (ISS); arid-focused acquisition",
    },
    {
        "key": "emit_l2a_rfl_v001",
        "concept_id": "C2408750690-LPCLOUD",
        "sensor": "EMIT",
        "gsd_m": 60,
        "note": "EMIT L2A Estimated Reflectance V001 (legacy)",
    },
)

PLACEHOLDER_SENSORS = (
    {"key": "enmap_l2a", "sensor": "EnMAP", "gsd_m": 30, "status": "not_queried"},
    {"key": "prisma_l2", "sensor": "PRISMA", "gsd_m": 30, "status": "not_queried"},
)

OUT = ROOT / "results" / "emit_coverage.csv"
CACHE = ROOT / "data" / "gemx" / "cmr_cache"


def _cmr_page(concept_id: str, bbox: tuple[float, float, float, float], page: int) -> dict:
    west, south, east, north = bbox
    params = {
        "collection_concept_id": concept_id,
        "bounding_box": f"{west},{south},{east},{north}",
        "page_size": 50,
        "page_num": page,
    }
    url = CMR_GRANULE_URL + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _granule_polygon(umm: dict):
    """Best-effort footprint polygon from UMM SpatialExtent."""
    try:
        ext = umm["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"]
    except (KeyError, TypeError):
        return None
    if "GPolygons" in ext:
        rings = ext["GPolygons"][0]["Boundary"]["Points"]
        coords = [(p["Longitude"], p["Latitude"]) for p in rings]
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        from shapely.geometry import Polygon

        return Polygon(coords)
    if "BoundingRectangles" in ext:
        br = ext["BoundingRectangles"][0]
        return box(br["WestBoundingCoordinate"], br["SouthBoundingCoordinate"],
                   br["EastBoundingCoordinate"], br["NorthBoundingCoordinate"])
    return None


def coverage_for_site(
    site_id: str,
    bbox: tuple[float, float, float, float],
    collection: dict,
    *,
    use_cache: bool = True,
) -> dict:
    cache_path = CACHE / f"{site_id}__{collection['key']}.json"
    CACHE.mkdir(parents=True, exist_ok=True)
    if use_cache and cache_path.is_file():
        payload = json.loads(cache_path.read_text())
    else:
        items = []
        page = 1
        while page <= 20:
            try:
                data = _cmr_page(collection["concept_id"], bbox, page)
            except Exception as exc:
                return {
                    "site_id": site_id,
                    "collection": collection["key"],
                    "sensor": collection["sensor"],
                    "coverage_frac": 0.0,
                    "n_granules": 0,
                    "status": f"error:{exc.__class__.__name__}",
                    "gsd_m": collection["gsd_m"],
                }
            batch = data.get("items") or []
            if not batch:
                break
            items.extend(batch)
            hits = int(data.get("hits") or 0)
            if len(items) >= hits or len(batch) < 50:
                break
            page += 1
            time.sleep(0.2)
        payload = {"items": items, "n": len(items)}
        cache_path.write_text(json.dumps(payload))

    items = payload.get("items") or []
    site_poly = box(*bbox)
    geoms = []
    for it in items:
        umm = it.get("umm") or it
        g = _granule_polygon(umm)
        if g is not None and not g.is_empty:
            geoms.append(g)
    if not geoms:
        return {
            "site_id": site_id,
            "collection": collection["key"],
            "sensor": collection["sensor"],
            "coverage_frac": 0.0,
            "n_granules": 0,
            "status": "none",
            "gsd_m": collection["gsd_m"],
        }

    import geopandas as gpd

    union = unary_union(geoms).intersection(site_poly)
    site_g = gpd.GeoSeries([site_poly], crs="EPSG:4326").to_crs(EQUAL_AREA_CRS)
    hit_g = gpd.GeoSeries([union], crs="EPSG:4326").to_crs(EQUAL_AREA_CRS)
    area_site = float(site_g.area.iloc[0])
    area_hit = float(hit_g.area.iloc[0]) if not union.is_empty else 0.0
    frac = min(area_hit / area_site, 1.0) if area_site > 0 else 0.0
    status = "covered" if frac >= 0.05 else ("trace" if frac > 0 else "none")
    return {
        "site_id": site_id,
        "collection": collection["key"],
        "sensor": collection["sensor"],
        "coverage_frac": round(frac, 4),
        "n_granules": len(items),
        "status": status,
        "gsd_m": collection["gsd_m"],
    }


def build_inventory(*, site_ids: list[str] | None = None, use_cache: bool = True) -> pd.DataFrame:
    site_ids = site_ids or list_site_ids(ROOT / "sites")
    rows = []
    for i, sid in enumerate(site_ids):
        site = load_site_by_id(sid, ROOT / "sites")
        bbox = tuple(site.bbox_wgs84)
        print(f"[{i+1}/{len(site_ids)}] EMIT {sid}", flush=True)
        for coll in EMIT_COLLECTIONS:
            rows.append(coverage_for_site(sid, bbox, coll, use_cache=use_cache))
        for ph in PLACEHOLDER_SENSORS:
            rows.append({
                "site_id": sid,
                "collection": ph["key"],
                "sensor": ph["sensor"],
                "coverage_frac": 0.0,
                "n_granules": 0,
                "status": ph["status"],
                "gsd_m": ph["gsd_m"],
            })
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--site", action="append", dest="sites")
    p.add_argument("--refresh", action="store_true", help="Ignore CMR cache")
    p.add_argument("--limit", type=int, default=0)
    args = p.parse_args()
    sites = args.sites
    if sites is None and args.limit:
        sites = list_site_ids(ROOT / "sites")[: args.limit]
    df = build_inventory(site_ids=sites, use_cache=not args.refresh)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    emit = df[df["sensor"] == "EMIT"]
    print(f"Wrote {args.out.relative_to(ROOT)} ({len(df)} rows)")
    print("EMIT status:\n", emit["status"].value_counts().to_string())
    print(f"EMIT covered (frac>=0.05): {(emit['status']=='covered').sum()} / {len(emit)}")


if __name__ == "__main__":
    main()
