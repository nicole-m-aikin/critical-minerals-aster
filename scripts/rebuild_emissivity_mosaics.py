#!/usr/bin/env python3
"""Rebuild site index mosaics from AST_05 v004 surface emissivity.

Plan (CMR metadata only):
    python scripts/rebuild_emissivity_mosaics.py --all
    -> results/emissivity_mosaic_plan.csv

Prefetch granules only (EarthData via ~/.netrc):
    python scripts/rebuild_emissivity_mosaics.py --all --prefetch

Build mosaics + per-site pipeline:
    python scripts/rebuild_emissivity_mosaics.py --all --execute --skip-figures --workers 4

Night replication (separate directory, never feeds the primary score):
    python scripts/rebuild_emissivity_mosaics.py --all --execute --night

Day builds move the current radiance (AST_L1T) granules and mosaics to
``data/sites/{id}/aster_l1t_day/`` before writing the emissivity mosaic into
``aster/``. Night builds write ``data/sites/{id}/aster_night/`` and do not run
the pipeline (see scripts/day_night_replication.py).
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from critical_minerals_aster.config import list_site_ids, load_site_by_id  # noqa: E402
from critical_minerals_aster.paths import site_paths_for  # noqa: E402
from critical_minerals_aster.pipeline import (  # noqa: E402
    download_and_mosaic_ast05,
    plan_ast05_selection,
    run_site,
)
from critical_minerals_aster.spectral import extract_granule_id  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning, module="earthaccess")
warnings.filterwarnings("ignore", message=".*DataGranule.size.*")

RESULTS = REPO / "results"
L1T_ARCHIVE = "aster_l1t_day"
NIGHT_DIR = "aster_night"


def plan_path(night: bool) -> Path:
    return RESULTS / ("emissivity_mosaic_plan_night.csv" if night else "emissivity_mosaic_plan.csv")


def night_dir(site_id: str) -> Path:
    return REPO / "data" / "sites" / site_id / NIGHT_DIR


def plan_site(site_id: str, day_night: str) -> dict:
    site = load_site_by_id(site_id, REPO / "sites")
    sel = plan_ast05_selection(site, day_night=day_night)
    clouds = [r["cloud_pct"] for r in sel.records if r["cloud_pct"] is not None]
    return {
        "site_id": site_id,
        "day_night": day_night,
        "tier": sel.tier,
        "n_granules": len(sel.granules),
        "union_coverage": round(sel.union_coverage, 3),
        "max_cloud_pct": max(clouds) if clouds else None,
        "months": ";".join(str(r["month"]) for r in sel.records),
        "granules": ";".join(r["granule_id"] for r in sel.records),
    }


def write_plan(site_ids: list[str], day_night: str) -> pd.DataFrame:
    rows = []
    for i, sid in enumerate(site_ids, 1):
        try:
            r = plan_site(sid, day_night)
            print(f"[{i}/{len(site_ids)}] {sid}: {r['tier']} x{r['n_granules']} "
                  f"cov={r['union_coverage']:.2f} months={r['months']}", flush=True)
        except Exception as exc:
            print(f"[{i}/{len(site_ids)}] {sid}: plan failed ({exc})", file=sys.stderr)
            r = {"site_id": sid, "day_night": day_night, "tier": f"error: {exc}"}
        rows.append(r)
    plan = pd.DataFrame(rows)
    out = plan_path(day_night == "Night")
    plan.to_csv(out, index=False)
    print(f"Wrote {out.relative_to(REPO)}")
    return plan


def prefetch_site(site_id: str, day_night: str) -> str:
    import earthaccess

    earthaccess.login(strategy="netrc")
    site = load_site_by_id(site_id, REPO / "sites")
    paths = site_paths_for(site, REPO)
    dest = night_dir(site_id) if day_night == "Night" else paths.aster_dir
    sel = plan_ast05_selection(site, day_night=day_night)
    dest.mkdir(parents=True, exist_ok=True)
    missing = [
        g for g in sel.granules
        if not (dest / f"{extract_granule_id(g)}_SRE_TIR_QA_DataPlane.tif").is_file()
    ]
    if missing:
        earthaccess.download(missing, str(dest), show_progress=False)
    return f"{len(sel.granules)} granules ({len(missing)} downloaded)"


def _archive_l1t(aster_dir: Path) -> int:
    def _is_ast05_mosaic(p: Path) -> bool:
        sel = aster_dir / f"{p.name.split('_mosaic')[0]}_mosaic_granules.json"
        return sel.is_file() and json.loads(sel.read_text()).get("product") == "AST_05"

    files = [p for p in aster_dir.iterdir() if p.is_file()
             and ("AST_L1T" in p.name or ("_mosaic" in p.name and not _is_ast05_mosaic(p)))]
    if not files:
        return 0
    archive = aster_dir.parent / L1T_ARCHIVE
    archive.mkdir(exist_ok=True)
    for p in files:
        dest = archive / p.name
        if dest.exists():
            p.unlink()
        else:
            shutil.move(str(p), dest)
    return len(files)


def build_site(site_id: str, day_night: str, skip_figures: bool) -> str:
    import earthaccess

    earthaccess.login(strategy="netrc")
    site = load_site_by_id(site_id, REPO / "sites")
    paths = site_paths_for(site, REPO)
    if day_night == "Night":
        out = night_dir(site_id)
        download_and_mosaic_ast05(site, paths, interactive_login=False, day_night="Night",
                                  out_dir=out, mosaic_id=f"{site_id}_night_mosaic")
        return "night mosaic built"
    paths.aster_dir.mkdir(parents=True, exist_ok=True)
    moved = _archive_l1t(paths.aster_dir)
    gid = download_and_mosaic_ast05(site, paths, interactive_login=False)
    run_site(site, REPO, download=False, skip_figures=skip_figures)
    return f"{gid} ({moved} L1T files archived)"


def _run_parallel(fn, site_ids: list[str], workers: int, *args) -> list[str]:
    failed = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fn, sid, *args): sid for sid in site_ids}
        for i, fut in enumerate(as_completed(futs), 1):
            sid = futs[fut]
            try:
                print(f"[{i}/{len(site_ids)}] {sid}: OK {fut.result()}", flush=True)
            except Exception as exc:
                failed.append(sid)
                print(f"[{i}/{len(site_ids)}] {sid}: FAILED ({exc!r})", file=sys.stderr, flush=True)
    return failed


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--all", action="store_true")
    g.add_argument("--sites", nargs="+")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--prefetch", action="store_true", help="Download selected granules only")
    mode.add_argument("--execute", action="store_true", help="Build mosaics (+ pipeline for day)")
    ap.add_argument("--night", action="store_true", help="Night replication instead of day")
    ap.add_argument("--skip-figures", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args(argv)

    site_ids = list_site_ids(REPO / "sites") if args.all else args.sites
    day_night = "Night" if args.night else "Day"
    if args.prefetch:
        import earthaccess

        earthaccess.login(strategy="netrc")
        failed = _run_parallel(prefetch_site, site_ids, args.workers, day_night)
    elif args.execute:
        failed = _run_parallel(build_site, site_ids, args.workers, day_night, args.skip_figures)
    else:
        write_plan(site_ids, day_night)
        return 0
    if failed:
        print(f"FAILED: {' '.join(failed)}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
