"""Download USGS Quaternary fault data for study sites and wire into site YAMLs."""

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

SITES = {
    "mcdermitt":    [-118.1, 41.8, -117.3, 42.4],
    "silver_peak":  [-118.1, 37.5, -117.3, 38.1],
    "goldfield":    [-117.4, 37.5, -116.8, 38.1],
    "tonopah":      [-117.4, 37.8, -116.8, 38.3],
    "yerington":    [-119.4, 38.7, -118.8, 39.1],
    "mountain_pass":[-115.7, 35.3, -115.3, 35.7],
    "climax":       [-106.5, 39.3, -106.1, 39.6],
    "jerome":       [-112.4, 34.5, -112.0, 34.9],
    "stillwater":   [-118.2, 39.4, -117.8, 39.7],
}

# USGS Quaternary Fault and Fold Database - layer 21 = National Database (polylines)
ENDPOINTS = [
    "https://earthquake.usgs.gov/arcgis/rest/services/haz/Qfaults/MapServer/21/query",
]


def build_url(base: str, west: float, south: float, east: float, north: float) -> str:
    params = {
        "where": "1=1",
        "geometry": f"{west},{south},{east},{north}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "*",
        "f": "geojson",
    }
    return base + "?" + urllib.parse.urlencode(params)


def fetch_geojson(url: str, timeout: int = 30) -> dict | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "critical-minerals-aster/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        if "error" in data:
            print(f"  API error: {data['error']}")
            return None
        return data
    except Exception as exc:
        print(f"  Request failed: {exc}")
        return None


def ensure_type_property(geojson: dict) -> dict:
    """Add type='faults' property to any feature missing it."""
    for feat in geojson.get("features", []):
        if "properties" not in feat or feat["properties"] is None:
            feat["properties"] = {}
        feat["properties"].setdefault("type", "faults")
    return geojson


def download_faults(site_id: str, bbox: list[float]) -> Path | None:
    west, south, east, north = bbox
    print(f"\n[{site_id}] bbox={bbox}")
    for endpoint in ENDPOINTS:
        url = build_url(endpoint, west, south, east, north)
        print(f"  Trying: {endpoint.split('/')[2]}...")
        data = fetch_geojson(url)
        if data is None:
            continue
        features = data.get("features", [])
        print(f"  Got {len(features)} features")
        if not features:
            continue
        data = ensure_type_property(data)
        out_path = STRUCTURES_DIR / f"{site_id}_faults.geojson"
        out_path.write_text(json.dumps(data, indent=2))
        size_kb = out_path.stat().st_size / 1024
        print(f"  Saved -> {out_path.relative_to(REPO_ROOT)} ({size_kb:.1f} KB)")
        return out_path
    print(f"  No fault features found for {site_id}")
    return None


def update_yaml(site_id: str, geojson_path: Path) -> None:
    """Patch structure_layers in-place using text replacement to preserve all other formatting."""
    yaml_path = SITES_DIR / f"{site_id}.yaml"
    if not yaml_path.exists():
        print(f"  WARN: No YAML found at {yaml_path}")
        return

    text = yaml_path.read_text()
    rel_path = str(geojson_path.relative_to(REPO_ROOT))

    # Skip if already wired up
    if rel_path in text:
        print(f"  YAML already has entry for {rel_path}, skipping")
        return

    layer_block = (
        f"structure_layers:\n"
        f"  - path: {rel_path}\n"
        f"    type: faults\n"
        f"    buffer_m: 500\n"
    )

    if "structure_layers: []" in text:
        new_text = text.replace("structure_layers: []", layer_block.rstrip("\n"), 1)
    elif "structure_layers:" in text:
        # Already has layers - append using yaml round-trip (should be safe here)
        data = yaml.safe_load(text)
        existing_layers = data.get("structure_layers") or []
        paths_already = [l.get("path") for l in existing_layers if isinstance(l, dict)]
        if rel_path in paths_already:
            print(f"  YAML already has entry for {rel_path}, skipping")
            return
        # Text-level: append to the existing block
        new_text = text.rstrip() + (
            f"\n  - path: {rel_path}\n"
            f"    type: faults\n"
            f"    buffer_m: 500\n"
        )
    else:
        new_text = text.rstrip() + "\n" + layer_block

    yaml_path.write_text(new_text)
    print(f"  Updated {yaml_path.relative_to(REPO_ROOT)}")


def main() -> None:
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

    successes = []
    failures = []

    for site_id, bbox in SITES.items():
        path = download_faults(site_id, bbox)
        if path:
            update_yaml(site_id, path)
            successes.append(site_id)
        else:
            failures.append(site_id)

    print("\n" + "=" * 60)
    print(f"SUCCESS ({len(successes)}): {', '.join(successes) or 'none'}")
    print(f"NO DATA ({len(failures)}): {', '.join(failures) or 'none'}")


if __name__ == "__main__":
    main()
