"""Fetch fault GeoJSONs for all 15 new expansion sites (rounds 1-3, reaching 30 total).

Strategy per site:
  1. Try USGS Quaternary Fault and Fold Database (layer 21) — fast, accurate for active western US.
  2. If zero features returned, auto-fall back to USGS SGMC FeatureServer (all geological ages).

Site YAMLs are already pre-wired with the correct data/structures/ paths; this script
only downloads the missing GeoJSON files.

Usage:
    python scripts/fetch_new_site_structures.py
    python scripts/fetch_new_site_structures.py --site bagdad      # single site
    python scripts/fetch_new_site_structures.py --skip-existing    # skip already-downloaded
"""

from __future__ import annotations

import argparse
import json
import urllib.request
import urllib.parse
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
STRUCTURES_DIR = REPO_ROOT / "data" / "structures"

# ---------------------------------------------------------------------------
# New sites: id -> [west, south, east, north]
# ---------------------------------------------------------------------------
NEW_SITES: dict[str, list[float]] = {
    # Porphyry Cu replication (AZ/NV/UT)
    "bagdad":           [-113.45, 34.35, -112.95, 34.75],
    "globe_miami":      [-111.05, 33.20, -110.55, 33.65],
    "ajo":              [-113.10, 32.15, -112.60, 32.60],
    "ely_robinson":     [-115.25, 39.00, -114.65, 39.45],
    "bingham_canyon":   [-112.35, 40.35, -111.90, 40.70],
    # Clay Li / carbonatite REE
    "thacker_pass":     [-118.20, 41.70, -117.85, 41.95],
    "elk_creek":        [-96.35,  40.10, -95.90,  40.50],
    # Negative controls
    "carlin_trend":     [-116.60, 40.90, -116.10, 41.50],
    "gas_hills":        [-107.90, 42.65, -107.40, 43.05],
    "green_river":      [-109.80, 41.30, -109.20, 41.80],
    "viburnum_trend":   [-91.40,  37.50, -90.90,  38.00],
    # Skarn
    "bisbee":           [-110.15, 31.25, -109.65, 31.70],
    "hanover_fierro":   [-108.25, 32.65, -107.80, 33.05],
    # IOCG REE / REE veins
    "pea_ridge":        [-90.90,  37.70, -90.50,  38.10],
    "lemhi_pass":       [-114.15, 44.75, -113.65, 45.20],
}

# ---------------------------------------------------------------------------
# API endpoints
# ---------------------------------------------------------------------------
QFAULTS_URL = (
    "https://earthquake.usgs.gov/arcgis/rest/services"
    "/haz/Qfaults/MapServer/21/query"
)

# Primary SGMC endpoint (FeatureServer with RuleID fault filtering)
SGMC_URL = (
    "https://services.arcgis.com/v01gqwM5QqNysAAi/arcgis/rest/services"
    "/SB_5888bf4fe4b05ccb964bab9d_USGS_SGMC_feature/FeatureServer/1/query"
)
# Fallback SGMC endpoint (mrdata.usgs.gov MapServer — broader, no RuleID filter)
SGMC_MRDATA_URL = "https://mrdata.usgs.gov/arcgis/rest/services/sgmc/MapServer/1/query"

# SGMC RuleID values that correspond to fault-type linework (from Horton et al. 2017)
SGMC_FAULT_RULE_IDS = (
    "11,12,13,21,22,23,24,29,30,31,33,34,35,36,"
    "42,43,44,45,46,47,48,49,50,51,52,53,54,62,63,64,65,66"
)
SGMC_PAGE_SIZE = 2000

HEADERS = {"User-Agent": "critical-minerals-aster/1.0"}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def fetch_json(url: str, timeout: int = 45) -> dict | None:
    try:
        req = urllib.request.Request(url, headers=HEADERS)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        if "error" in data:
            print(f"      API error: {data['error']}")
            return None
        return data
    except Exception as exc:
        print(f"      Request failed: {exc}")
        return None


def keep_line_geometries(features: list[dict]) -> list[dict]:
    return [
        f for f in features
        if f.get("geometry", {}).get("type", "") in ("LineString", "MultiLineString")
    ]


def tag_type(features: list[dict]) -> list[dict]:
    for f in features:
        if "properties" not in f or f["properties"] is None:
            f["properties"] = {}
        f["properties"].setdefault("type", "faults")
    return features


def save_geojson(features: list[dict], out_path: Path) -> None:
    fc = {"type": "FeatureCollection", "features": features}
    out_path.write_text(json.dumps(fc, indent=2))
    size_kb = out_path.stat().st_size / 1024
    print(f"      Saved {out_path.name} ({size_kb:.1f} KB, {len(features)} features)")


# ---------------------------------------------------------------------------
# USGS Quaternary Faults fetch
# ---------------------------------------------------------------------------

def fetch_quaternary(west: float, south: float, east: float, north: float) -> list[dict]:
    params = {
        "where": "1=1",
        "geometry": f"{west},{south},{east},{north}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "*",
        "f": "geojson",
    }
    url = QFAULTS_URL + "?" + urllib.parse.urlencode(params)
    data = fetch_json(url)
    if data is None:
        return []
    return data.get("features", [])


# ---------------------------------------------------------------------------
# USGS SGMC fetch (paginated)
# ---------------------------------------------------------------------------

def fetch_sgmc(west: float, south: float, east: float, north: float) -> list[dict]:
    """Paginated fetch from the primary SGMC FeatureServer with RuleID fault filter."""
    all_features: list[dict] = []
    offset = 0
    page = 0
    while True:
        page += 1
        params = {
            "where": f"RuleID IN ({SGMC_FAULT_RULE_IDS})",
            "geometry": f"{west},{south},{east},{north}",
            "geometryType": "esriGeometryEnvelope",
            "inSR": "4326",
            "outSR": "4326",
            "spatialRel": "esriSpatialRelIntersects",
            "outFields": "STATE,DESCRIPTION,RuleID",
            "resultOffset": offset,
            "resultRecordCount": SGMC_PAGE_SIZE,
            "f": "geojson",
        }
        url = SGMC_URL + "?" + urllib.parse.urlencode(params)
        data = fetch_json(url)
        if data is None:
            break
        batch = data.get("features", [])
        all_features.extend(batch)
        print(f"      SGMC page {page}: {len(batch)} features")
        if len(batch) < SGMC_PAGE_SIZE:
            break
        offset += SGMC_PAGE_SIZE
    return all_features


def fetch_sgmc_mrdata(west: float, south: float, east: float, north: float) -> list[dict]:
    """Fallback fetch from mrdata.usgs.gov MapServer — no RuleID filter, broader coverage."""
    params = {
        "where": "1=1",
        "geometry": f"{west},{south},{east},{north}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": "4326",
        "outSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "*",
        "f": "geojson",
    }
    url = SGMC_MRDATA_URL + "?" + urllib.parse.urlencode(params)
    data = fetch_json(url)
    if data is None:
        return []
    return data.get("features", [])


# ---------------------------------------------------------------------------
# Per-site orchestration
# ---------------------------------------------------------------------------

def process_site(site_id: str, bbox: list[float], skip_existing: bool) -> str:
    """Returns 'ok', 'skipped', or 'empty'."""
    west, south, east, north = bbox
    out_path = STRUCTURES_DIR / f"{site_id}_faults.geojson"

    if skip_existing and out_path.exists():
        print(f"  [{site_id}] already exists, skipping")
        return "skipped"

    print(f"\n  [{site_id}] bbox={bbox}")

    # --- 1. Try Quaternary ---
    print("      Trying USGS Quaternary Faults...")
    features = fetch_quaternary(west, south, east, north)
    features = keep_line_geometries(tag_type(features))
    source = "Quaternary"

    # --- 2. Fall back to SGMC FeatureServer ---
    if not features:
        print("      Zero Quaternary features — trying SGMC FeatureServer (all ages)...")
        features = fetch_sgmc(west, south, east, north)
        features = keep_line_geometries(tag_type(features))
        source = "SGMC-FeatureServer"

    # --- 3. Fall back to SGMC mrdata.usgs.gov ---
    if not features:
        print("      Zero SGMC FeatureServer features — trying SGMC mrdata fallback...")
        features = fetch_sgmc_mrdata(west, south, east, north)
        features = keep_line_geometries(tag_type(features))
        source = "SGMC-mrdata"

    if not features:
        print(f"      No fault features from any source for {site_id}")
        # Write an empty FeatureCollection so the pipeline doesn't crash on a missing file
        save_geojson([], out_path)
        return "empty"

    print(f"      Source: {source} — {len(features)} fault features")
    save_geojson(features, out_path)
    return "ok"


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site", help="Process a single site by id")
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip sites whose GeoJSON already exists in data/structures/",
    )
    args = parser.parse_args()

    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

    sites = NEW_SITES
    if args.site:
        if args.site not in NEW_SITES:
            print(f"Unknown site '{args.site}'. Available: {', '.join(NEW_SITES)}")
            raise SystemExit(1)
        sites = {args.site: NEW_SITES[args.site]}

    results: dict[str, list[str]] = {"ok": [], "skipped": [], "empty": []}
    for site_id, bbox in sites.items():
        status = process_site(site_id, bbox, skip_existing=args.skip_existing)
        results[status].append(site_id)

    print("\n" + "=" * 60)
    print(f"OK       ({len(results['ok'])}): {', '.join(results['ok']) or 'none'}")
    print(f"EMPTY    ({len(results['empty'])}): {', '.join(results['empty']) or 'none'}")
    print(f"SKIPPED  ({len(results['skipped'])}): {', '.join(results['skipped']) or 'none'}")
    if results["empty"]:
        print(
            "\nNOTE: Empty files written for sites with no fault coverage — pipeline will "
            "run without structure overlay for those sites. This is expected for some "
            "Midwest sites (Elk Creek, Pea Ridge, Viburnum Trend)."
        )


if __name__ == "__main__":
    main()
