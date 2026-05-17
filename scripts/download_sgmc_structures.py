"""Download USGS SGMC structural fault data for sites that lack Quaternary fault coverage.

Targets Bear Lodge (WY) and Jerome (AZ), which returned zero features from the
USGS Quaternary Faults API because both are geologically older terranes.

SGMC_Structure FeatureServer covers all geological ages and all 48 contiguous
states (Horton et al. 2017, USGS Data Series 1052; updated Dec 2025).

For Bear Lodge only, a second layer is downloaded from the Wyoming State
Geological Survey (WSGS) Precambrian Basement Map, which captures the
basement-involved fault fabric that directly controlled REE carbonatite
emplacement at the Black Hills margin.
"""

from __future__ import annotations

import json
import sys
import urllib.request
import urllib.parse
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).parent.parent
STRUCTURES_DIR = REPO_ROOT / "data" / "structures"
SITES_DIR = REPO_ROOT / "sites"

# ---------------------------------------------------------------------------
# SGMC FeatureServer — layer 1 = structural lines (faults, contacts, folds)
# ---------------------------------------------------------------------------
SGMC_URL = (
    "https://services.arcgis.com/v01gqwM5QqNysAAi/arcgis/rest/services"
    "/SB_5888bf4fe4b05ccb964bab9d_USGS_SGMC_feature/FeatureServer/1/query"
)

# RuleID values that correspond to fault-type features (excludes contacts,
# folds, dikes, glacial features). Derived from Horton et al. 2017 symbology.
FAULT_RULE_IDS = (
    "11,12,13,21,22,23,24,29,30,31,33,34,35,36,"
    "42,43,44,45,46,47,48,49,50,51,52,53,54,62,63,64,65,66"
)
SGMC_PAGE_SIZE = 2000

# ---------------------------------------------------------------------------
# WSGS Precambrian Basement Map — faults only (Bear Lodge supplement)
# ---------------------------------------------------------------------------
WSGS_BASEMENT_URL = (
    "https://portal.wsgs.wyo.gov/ags/rest/services"
    "/OilGas/PrecambrianBasement_WY/MapServer/1/query"
)

# ---------------------------------------------------------------------------
# Sites to process
# ---------------------------------------------------------------------------
SITES: dict[str, list[float]] = {
    "bear_lodge": [-104.52, 44.32, -103.90, 44.78],
    "jerome":     [-112.30, 34.52, -111.88, 35.00],
    # Stillwater Complex, MT: Precambrian layered intrusion — no USGS Quaternary faults, use SGMC.
    "stillwater": [-110.42, 45.18, -109.70, 45.72],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fetch_geojson(url: str, timeout: int = 30) -> dict | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "critical-minerals-aster/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        if "error" in data:
            print(f"    API error: {data['error']}")
            return None
        return data
    except Exception as exc:
        print(f"    Request failed: {exc}")
        return None


def ensure_type_property(geojson: dict) -> dict:
    """Add type='faults' property to any feature missing it."""
    for feat in geojson.get("features", []):
        if "properties" not in feat or feat["properties"] is None:
            feat["properties"] = {}
        feat["properties"].setdefault("type", "faults")
    return geojson


def keep_line_geometries(features: list[dict]) -> list[dict]:
    """Drop any non-line features that slipped through (Points, Polygons)."""
    return [
        f for f in features
        if f.get("geometry", {}).get("type", "") in ("LineString", "MultiLineString")
    ]


def save_geojson(features: list[dict], out_path: Path) -> None:
    fc = {"type": "FeatureCollection", "features": features}
    out_path.write_text(json.dumps(fc, indent=2))
    size_kb = out_path.stat().st_size / 1024
    print(f"    Saved -> {out_path.relative_to(REPO_ROOT)} ({size_kb:.1f} KB, {len(features)} features)")


# ---------------------------------------------------------------------------
# SGMC fetch (paginated)
# ---------------------------------------------------------------------------

def fetch_sgmc_faults(west: float, south: float, east: float, north: float) -> list[dict]:
    all_features: list[dict] = []
    offset = 0
    page = 0
    while True:
        page += 1
        params = {
            "where": f"RuleID IN ({FAULT_RULE_IDS})",
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
        data = fetch_geojson(url)
        if data is None:
            break
        batch = data.get("features", [])
        all_features.extend(batch)
        print(f"    Page {page}: {len(batch)} features (total so far: {len(all_features)})")
        if len(batch) < SGMC_PAGE_SIZE:
            break
        offset += SGMC_PAGE_SIZE
    return all_features


# ---------------------------------------------------------------------------
# WSGS Precambrian Basement fetch (Bear Lodge supplement)
# ---------------------------------------------------------------------------

_WSGS_FAULT_KEYWORDS = ("fault",)


def fetch_wsgs_basement_faults(west: float, south: float, east: float, north: float) -> list[dict]:
    params = {
        "where": "1=1",
        "geometry": f"{west},{south},{east},{north}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": "4326",
        "outSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "Type",
        "f": "geojson",
    }
    url = WSGS_BASEMENT_URL + "?" + urllib.parse.urlencode(params)
    data = fetch_geojson(url)
    if data is None:
        return []

    raw = data.get("features", [])
    # Keep only fault-type features (the Type field may be a code or a label)
    faults = [
        f for f in raw
        if any(kw in str(f.get("properties", {}).get("Type", "")).lower()
               for kw in _WSGS_FAULT_KEYWORDS)
    ]
    print(f"    WSGS raw features: {len(raw)} → fault-type kept: {len(faults)}")
    return faults


# ---------------------------------------------------------------------------
# YAML update — adds a new layer entry without clobbering existing ones
# ---------------------------------------------------------------------------

def update_yaml(site_id: str, geojson_path: Path, label: str) -> None:
    yaml_path = SITES_DIR / f"{site_id}.yaml"
    if not yaml_path.exists():
        print(f"    WARN: No YAML at {yaml_path}")
        return

    text = yaml_path.read_text()
    rel_path = str(geojson_path.relative_to(REPO_ROOT))

    if rel_path in text:
        print(f"    YAML already has entry for {rel_path}, skipping")
        return

    layer_entry = (
        f"  - path: {rel_path}\n"
        f"    type: faults\n"
        f"    label: \"{label}\"\n"
        f"    buffer_m: 500\n"
    )

    if "structure_layers: []" in text:
        block = "structure_layers:\n" + layer_entry.rstrip("\n")
        new_text = text.replace("structure_layers: []", block, 1)
    elif "structure_layers:" in text:
        new_text = text.rstrip() + "\n" + layer_entry
    else:
        new_text = text.rstrip() + "\nstructure_layers:\n" + layer_entry

    yaml_path.write_text(new_text)
    print(f"    Updated {yaml_path.relative_to(REPO_ROOT)}")


# ---------------------------------------------------------------------------
# Per-site orchestration
# ---------------------------------------------------------------------------

def process_bear_lodge() -> None:
    site_id = "bear_lodge"
    west, south, east, north = SITES[site_id]
    print(f"\n[{site_id}] Fetching SGMC structural faults...")

    features = fetch_sgmc_faults(west, south, east, north)
    features = keep_line_geometries(ensure_type_property({"features": features})["features"])
    if features:
        out = STRUCTURES_DIR / f"{site_id}_faults_sgmc.geojson"
        save_geojson(features, out)
        update_yaml(site_id, out, "SGMC Structural Faults (all ages)")
    else:
        print(f"    No SGMC fault features for {site_id}")

    print(f"\n[{site_id}] Fetching WSGS Precambrian Basement faults...")
    wsgs_features = fetch_wsgs_basement_faults(west, south, east, north)
    wsgs_features = keep_line_geometries(wsgs_features)
    for f in wsgs_features:
        if "properties" not in f or f["properties"] is None:
            f["properties"] = {}
        f["properties"].setdefault("type", "faults")

    if wsgs_features:
        out = STRUCTURES_DIR / f"{site_id}_faults_basement.geojson"
        save_geojson(wsgs_features, out)
        update_yaml(site_id, out, "WSGS Precambrian Basement Faults")
    else:
        print(f"    No WSGS basement fault features for {site_id} — skipping second layer")


def process_jerome() -> None:
    site_id = "jerome"
    west, south, east, north = SITES[site_id]
    print(f"\n[{site_id}] Fetching SGMC structural faults...")

    features = fetch_sgmc_faults(west, south, east, north)
    features = keep_line_geometries(ensure_type_property({"features": features})["features"])
    if features:
        out = STRUCTURES_DIR / f"{site_id}_faults_sgmc.geojson"
        save_geojson(features, out)
        update_yaml(site_id, out, "SGMC Structural Faults (all ages)")
    else:
        print(f"    No SGMC fault features for {site_id}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def process_sgmc_only(site_id: str) -> None:
    """Generic SGMC-only download (no secondary WSGS layer)."""
    west, south, east, north = SITES[site_id]
    print(f"\n[{site_id}] Fetching SGMC structural faults...")

    features = fetch_sgmc_faults(west, south, east, north)
    features = keep_line_geometries(ensure_type_property({"features": features})["features"])
    if features:
        out = STRUCTURES_DIR / f"{site_id}_faults.geojson"
        save_geojson(features, out)
        update_yaml(site_id, out, "SGMC Structural Faults (all ages)")
    else:
        print(f"    No SGMC fault features for {site_id}")


def main() -> None:
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

    process_bear_lodge()
    process_jerome()
    process_sgmc_only("stillwater")

    print("\nDone.")


if __name__ == "__main__":
    main()
