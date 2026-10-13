#!/usr/bin/env python3
"""Freeze the feature-ready status table from the GEMx census manifest.

Writes ``results/gemx_census_feature_ready.csv`` with columns shared by the
prospectivity feature store:

- site_id, collection, coverage_frac, status (census)
- tir_ok, hsi_ok, hsi_missing, feature_ready

A site is feature-ready when TIR mosaics exist. HSI is optional: missing SWIR
sets ``hsi_missing=1`` (TIR-only domain) and still allows feature-store rows.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.config import list_site_ids, load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.pipeline import resolve_granule_id
from gemx_census_manifest import build_manifest
from gemx_swir_pilot import swir_tifs_exist

OUT = ROOT / "results" / "gemx_census_feature_ready.csv"


def _tir_ok(site_id: str) -> bool:
    try:
        site = load_site_by_id(site_id, ROOT / "sites")
        paths = site_paths_for(site, ROOT)
        gid = resolve_granule_id(site, paths)
        return (paths.aster_dir / f"{gid}_ratio_silica.tif").is_file()
    except Exception:
        return False


def freeze(*, min_free_gb: float = 40.0, out: Path = OUT) -> pd.DataFrame:
    man = build_manifest(min_free_gb=min_free_gb)
    # Persist refreshed manifest alongside freeze
    man.to_csv(ROOT / "results" / "gemx_census_manifest.csv", index=False)

    rows = []
    for r in man.itertuples(index=False):
        tir = _tir_ok(r.site_id)
        hsi = bool(r.have_tifs) and swir_tifs_exist(r.site_id)
        # Census "complete" preferred; also accept have_tifs for pilots
        if r.status == "complete" or hsi:
            hsi_status = "ok"
        elif r.status == "no_coverage":
            hsi_status = "none"
        elif r.status in {"pending", "needs_reprocess"}:
            hsi_status = "pending"
        else:
            hsi_status = str(r.status)
        hsi_missing = int(not hsi)
        feature_ready = bool(tir)  # TIR required; HSI optional
        rows.append({
            "site_id": r.site_id,
            "collection": r.collection or "",
            "coverage_frac": float(r.swir_coverage or 0.0),
            "census_status": r.status,
            "tir_ok": int(tir),
            "hsi_ok": int(hsi),
            "hsi_missing": hsi_missing,
            "hsi_status": hsi_status,
            "feature_ready": int(feature_ready),
            "sensor": getattr(r, "sensor", "") or "",
            "mosaic_resolution_m": getattr(r, "mosaic_resolution_m", "") or "",
        })

    df = pd.DataFrame(rows)
    # Ensure all configured sites appear even if manifest drifted
    known = set(df["site_id"])
    for sid in list_site_ids(ROOT / "sites"):
        if sid not in known:
            tir = _tir_ok(sid)
            df = pd.concat([df, pd.DataFrame([{
                "site_id": sid,
                "collection": "",
                "coverage_frac": 0.0,
                "census_status": "unknown",
                "tir_ok": int(tir),
                "hsi_ok": 0,
                "hsi_missing": 1,
                "hsi_status": "unknown",
                "feature_ready": int(tir),
                "sensor": "",
                "mosaic_resolution_m": "",
            }])], ignore_index=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    df = df.sort_values("site_id").reset_index(drop=True)
    df.to_csv(out, index=False)
    return df


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--min-free-gb", type=float, default=40.0)
    args = p.parse_args()
    df = freeze(min_free_gb=args.min_free_gb, out=args.out)
    print(f"Wrote {args.out.relative_to(ROOT)} ({len(df)} sites)")
    print("feature_ready:", int(df["feature_ready"].sum()))
    print("hsi_ok:", int(df["hsi_ok"].sum()))
    print("tir_ok:", int(df["tir_ok"].sum()))
    print(df["hsi_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
