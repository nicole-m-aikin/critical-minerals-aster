#!/usr/bin/env python3
"""Download structural geology data for 7 remaining sites."""
import json
import urllib.request
import urllib.parse
import os
import sys

SITES = {
    "oatman":           {"bbox": [-114.7, 34.8, -114.16, 35.4],   "try_qfaults": True},
    "marysvale":        {"bbox": [-112.5, 38.2, -111.97, 38.8],   "try_qfaults": True},
    "jerritt_canyon":   {"bbox": [-116.05, 41.20, -115.40, 41.70],"try_qfaults": True},
    "darwin":           {"bbox": [-117.80, 36.10, -117.40, 36.55],"try_qfaults": True},
    "steamboat_springs":{"bbox": [-119.90, 39.25, -119.55, 39.58],"try_qfaults": True},
    "bear_lodge":       {"bbox": [-104.52, 44.32, -103.90, 44.78],"try_qfaults": True},
    "jerome":           {"bbox": [-112.30, 34.52, -111.88, 35.00],"try_qfaults": False},
}

OUT_DIR = "/Users/nicoleaikin/projects/critical-minerals-aster/data/structures"
SITES_DIR = "/Users/nicoleaikin/projects/critical-minerals-aster/sites"

def build_qfaults_url(bbox):
    w, s, e, n = bbox
    geom = f"{w},{s},{e},{n}"
    return (
        f"https://earthquake.usgs.gov/arcgis/rest/services/haz/Qfaults/MapServer/21/query"
        f"?where=1%3D1"
        f"&geometry={urllib.parse.quote(geom)}"
        f"&geometryType=esriGeometryEnvelope"
        f"&inSR=4326"
        f"&spatialRel=esriSpatialRelIntersects"
        f"&outFields=*"
        f"&f=geojson"
    )

def build_sgmc_url(bbox, layer=1):
    w, s, e, n = bbox
    geom = f"{w},{s},{e},{n}"
    return (
        f"https://mrdata.usgs.gov/arcgis/rest/services/sgmc/MapServer/{layer}/query"
        f"?where=1%3D1"
        f"&geometry={urllib.parse.quote(geom)}"
        f"&geometryType=esriGeometryEnvelope"
        f"&inSR=4326"
        f"&spatialRel=esriSpatialRelIntersects"
        f"&outFields=*"
        f"&f=geojson"
    )

def fetch_geojson(url, timeout=30):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        # ArcGIS error response
        if "error" in data:
            return None, f"API error: {data['error']}"
        features = data.get("features", [])
        return data, len(features)
    except Exception as e:
        return None, str(e)

def add_type_property(geojson):
    for feat in geojson.get("features", []):
        if "properties" not in feat or feat["properties"] is None:
            feat["properties"] = {}
        feat["properties"]["type"] = "faults"
    return geojson

def save_geojson(site_id, geojson):
    path = os.path.join(OUT_DIR, f"{site_id}_faults.geojson")
    with open(path, "w") as f:
        json.dump(geojson, f, indent=2)
    return path

def update_yaml(site_id):
    yaml_path = os.path.join(SITES_DIR, f"{site_id}.yaml")
    with open(yaml_path, "r") as f:
        content = f.read()
    
    new_block = (
        f"structure_layers:\n"
        f"  - path: data/structures/{site_id}_faults.geojson\n"
        f"    type: faults\n"
        f"    buffer_m: 500"
    )
    
    if "structure_layers: []" in content:
        content = content.replace("structure_layers: []", new_block)
    elif "structure_layers:" in content:
        # Replace existing block — find and replace until next top-level key or EOF
        lines = content.splitlines(keepends=True)
        out = []
        in_block = False
        replaced = False
        for line in lines:
            if line.startswith("structure_layers:") and not replaced:
                out.append(new_block + "\n")
                in_block = True
                replaced = True
            elif in_block:
                # Skip old block lines (indented or list items)
                if line.startswith(" ") or line.startswith("-"):
                    continue
                else:
                    in_block = False
                    out.append(line)
            else:
                out.append(line)
        content = "".join(out)
    else:
        content = content.rstrip() + "\n" + new_block + "\n"
    
    with open(yaml_path, "w") as f:
        f.write(content)

results = {}

for site_id, cfg in SITES.items():
    bbox = cfg["bbox"]
    print(f"\n{'='*60}")
    print(f"Site: {site_id}")
    source = None
    geojson = None
    n_features = 0

    # Step 1: Try Quaternary Faults
    if cfg["try_qfaults"]:
        url = build_qfaults_url(bbox)
        print(f"  Trying Qfaults...")
        data, result = fetch_geojson(url)
        if data is not None and isinstance(result, int) and result > 0:
            print(f"  Qfaults: {result} features")
            geojson = data
            n_features = result
            source = "USGS Qfaults"
        else:
            print(f"  Qfaults: 0 features or error ({result})")

    # Step 2: Try SGMC if no data yet
    if geojson is None:
        for layer in [1, 0, 2, 3]:
            url = build_sgmc_url(bbox, layer)
            print(f"  Trying SGMC layer {layer}...")
            data, result = fetch_geojson(url)
            if data is not None and isinstance(result, int) and result > 0:
                print(f"  SGMC layer {layer}: {result} features")
                geojson = data
                n_features = result
                source = f"USGS SGMC layer {layer}"
                break
            else:
                print(f"  SGMC layer {layer}: 0 features or error ({result})")

    if geojson is not None and n_features > 0:
        geojson = add_type_property(geojson)
        path = save_geojson(site_id, geojson)
        update_yaml(site_id)
        print(f"  Saved: {path}")
        print(f"  YAML updated.")
        results[site_id] = {"source": source, "features": n_features, "status": "success"}
    else:
        print(f"  No data found for {site_id}")
        results[site_id] = {"source": None, "features": 0, "status": "no_data"}

print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
for site_id, r in results.items():
    if r["status"] == "success":
        print(f"  {site_id:25s} {r['features']:4d} features  [{r['source']}]")
    else:
        print(f"  {site_id:25s}   NO DATA")
