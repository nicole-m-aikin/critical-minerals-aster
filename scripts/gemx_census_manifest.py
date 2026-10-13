#!/usr/bin/env python3
"""Build the 109-site TIR–SWIR census download/process manifest.

Writes ``results/gemx_census_manifest.csv`` with one row per site:
selected collection, coverage, existing TIFF/cube state, and free-space gate.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.aviris import COLLECTION_GSD_M, COLLECTION_SENSORS
from critical_minerals_aster.config import load_site_by_id
from gemx_site_select import (
    AVIRIS5_COLLECTION,
    census_collection_by_site,
)
from gemx_swir_pilot import swir_tifs_exist


OUT = ROOT / "results" / "gemx_census_manifest.csv"
MIN_FREE_GB = 40.0


def _cube_bytes(site_id: str) -> int:
    d = ROOT / "data" / "gemx" / "aviris" / site_id
    if not d.is_dir():
        return 0
    total = 0
    for p in d.iterdir():
        if p.suffix.lower() in {".nc", ".bin", ".img"} or p.name.endswith(".rfl.bin"):
            total += p.stat().st_size
    return total


def _provenance_collection(site_id: str) -> str | None:
    path = ROOT / "results" / f"gemx_{site_id}_swir_provenance.json"
    if not path.is_file():
        # Legacy AVIRIS-5 pilots have no provenance file.
        if swir_tifs_exist(site_id):
            return AVIRIS5_COLLECTION
        return None
    return json.loads(path.read_text()).get("collection")


def build_manifest(*, min_free_gb: float = MIN_FREE_GB) -> pd.DataFrame:
    from critical_minerals_aster.config import list_site_ids

    site_ids = list_site_ids(ROOT / "sites")

    picked = census_collection_by_site(ROOT)
    free_gb = shutil.disk_usage(ROOT).free / 1e9
    rows = []
    for site_id in site_ids:
        sel = picked.get(site_id)
        have = swir_tifs_exist(site_id)
        prov = _provenance_collection(site_id)
        if sel is None:
            status = "no_coverage"
            collection = ""
            sensor = ""
            coverage = 0.0
            tifs_match = False
            action = "skip"
        else:
            collection, coverage, sensor = sel
            tifs_match = bool(have and prov == collection)
            if tifs_match:
                status = "complete"
                action = "reuse_tif"
            elif have and prov != collection:
                status = "needs_reprocess"
                action = "download_selected"
            else:
                status = "pending"
                action = "download_selected"
        cache = ROOT / "data" / "gemx" / "cmr_cache" / f"{site_id}__{collection}.json" if collection else None
        rows.append({
            "site_id": site_id,
            "status": status,
            "action": action,
            "collection": collection,
            "sensor": sensor,
            "swir_coverage": coverage,
            "mosaic_resolution_m": COLLECTION_GSD_M.get(collection, ""),
            "have_tifs": have,
            "tifs_match_selection": tifs_match,
            "existing_collection": prov or "",
            "cube_bytes": _cube_bytes(site_id),
            "cache_exists": bool(cache and cache.is_file()),
            "free_gb_at_build": round(free_gb, 1),
            "space_ok": free_gb >= min_free_gb,
        })
    df = pd.DataFrame(rows)
    df.attrs["min_free_gb"] = min_free_gb
    return df


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--min-free-gb", type=float, default=MIN_FREE_GB)
    p.add_argument("--out", type=Path, default=OUT)
    args = p.parse_args()
    df = build_manifest(min_free_gb=args.min_free_gb)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    n = len(df)
    print(f"Wrote {args.out.relative_to(ROOT)} ({n} sites)")
    print(df["status"].value_counts().to_string())
    print(df.loc[df["status"] != "no_coverage", "collection"].value_counts().to_string())
    free = shutil.disk_usage(ROOT).free / 1e9
    print(f"free_gb={free:.1f} gate>={args.min_free_gb}")


if __name__ == "__main__":
    main()
