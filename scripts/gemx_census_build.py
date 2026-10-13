#!/usr/bin/env python3
"""Resume-friendly GEMx SWIR census builder for all covered sites.

Reads ``results/gemx_census_manifest.csv``, processes sites needing download
or reprocess, and stops when free disk falls below the gate.
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from gemx_census_manifest import MIN_FREE_GB, OUT as MANIFEST, build_manifest
from gemx_swir_pilot import run_site



def cleanup_cubes(site_id: str) -> int:
    """Remove reflectance cubes after successful TIFF mosaics to free disk."""
    from gemx_swir_pilot import swir_tifs_exist
    if not swir_tifs_exist(site_id):
        return 0
    d = ROOT / "data" / "gemx" / "aviris" / site_id
    freed = 0
    for pth in d.iterdir():
        name = pth.name.lower()
        # Cubes, incomplete downloads, and ENVI companions
        if (
            pth.suffix.lower() == ".nc"
            or name.endswith(".bin")
            or name.endswith(".bin.part")
            or name.endswith(".img")
            or name.endswith(".rfl.bin")
            or ".rfl.bin" in name
        ):
            freed += pth.stat().st_size
            pth.unlink()
        elif name.endswith(".hdr") and "rfl" in name:
            # Keep hdr only if a matching bin remains; otherwise remove orphans
            stem = pth.name[:-4]  # strip .hdr
            bin_cand = pth.parent / f"{stem}.bin"
            part_cand = pth.parent / f"{stem}.bin.part"
            if not bin_cand.is_file() and not part_cand.is_file():
                freed += pth.stat().st_size
                pth.unlink()
    return freed


def free_gb() -> float:
    return shutil.disk_usage(ROOT).free / 1e9


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--min-free-gb", type=float, default=MIN_FREE_GB)
    p.add_argument("--limit", type=int, default=0, help="Max sites to process (0=all)")
    p.add_argument("--from-tif-only", action="store_true",
                   help="Only rebuild figures for complete TIFF sites")
    p.add_argument("--actions", default="download_selected,reuse_tif",
                   help="Comma-separated manifest actions to run")
    p.add_argument("--cleanup-only", action="store_true",
                   help="Delete cubes for sites with complete TIFFs and exit")
    args = p.parse_args()

    if args.cleanup_only:
        man = build_manifest(min_free_gb=args.min_free_gb)
        total = 0
        for sid in man.loc[man["have_tifs"], "site_id"]:
            freed = cleanup_cubes(sid)
            if freed:
                print(f"{sid}: freed {freed/1e9:.1f} GB")
                total += freed
        print(f"total freed {total/1e9:.1f} GB; free now {free_gb():.1f} GB")
        return

    # Refresh manifest first.
    build_manifest(min_free_gb=args.min_free_gb).to_csv(MANIFEST, index=False)
    df = pd.read_csv(MANIFEST)
    actions = {a.strip() for a in args.actions.split(",") if a.strip()}
    todo = df[df["action"].isin(actions)].copy()
    # Prefer smaller pending AVIRIS-5 first, then Classic; complete reuse last.
    order = {"reuse_tif": 0, "download_selected": 1}
    todo["_ord"] = todo["action"].map(order).fillna(9)
    todo["_coll"] = (todo["collection"] != "aviris5_l2a_rfl").astype(int)
    todo = todo.sort_values(["_ord", "_coll", "site_id"])
    if args.limit:
        todo = todo.head(args.limit)

    print(f"Queue {len(todo)} sites; free={free_gb():.1f} GB; gate>={args.min_free_gb}")
    rows = []
    for i, row in enumerate(todo.itertuples(index=False), 1):
        if free_gb() < args.min_free_gb:
            print(f"STOP: free disk {free_gb():.1f} GB < gate {args.min_free_gb}", flush=True)
            break
        site_id = row.site_id
        print(f"=== [{i}/{len(todo)}] {site_id} ({row.collection}) free={free_gb():.1f}G ===", flush=True)
        t0 = time.time()
        from_tif = bool(args.from_tif_only or row.action == "reuse_tif")
        result = run_site(site_id, from_tif=from_tif, skip_if_no_coverage=True)
        result["elapsed_s"] = round(time.time() - t0, 1)
        rows.append(result)
        print(f"  -> {result.get('status')} ({result.get('elapsed_s')}s)", flush=True)
        if result.get("status") in {"ok", "ok_from_tif"}:
            freed = cleanup_cubes(site_id)
            if freed:
                print(f"  freed {freed/1e9:.1f} GB cubes", flush=True)

    out = ROOT / "results" / "gemx_census_build_log.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    build_manifest(min_free_gb=args.min_free_gb).to_csv(MANIFEST, index=False)
    print(f"Wrote {out.relative_to(ROOT)}")
    print(pd.read_csv(MANIFEST)["status"].value_counts().to_string())


if __name__ == "__main__":
    main()
