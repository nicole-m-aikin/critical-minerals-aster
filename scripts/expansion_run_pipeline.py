"""Batched, error-tolerant pipeline runner for the +58 expansion sites.

Runs `critical_minerals_aster run --site <id> --mosaic` for each new site, one at a time,
catching and logging per-site failures (a bad granule search must not halt the batch).
Progress + outcomes are appended to scripts/expansion_run_pipeline.log so the run is
resumable and inspectable while it is going.

Usage:
    python scripts/expansion_run_pipeline.py                 # all 58, skipping done ones
    python scripts/expansion_run_pipeline.py --batch 0 8     # sites [0:8] only
    python scripts/expansion_run_pipeline.py --list          # print the site list + status

Designed to be launched with run_in_background; check the .log file for progress.
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import signal
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LOG = Path(__file__).with_suffix(".log")
HARD_TIMEOUT = 1200  # 20 min per site; sites were taking ~1-2 min, so this only trips on a hang

# The 58 expansion site_ids, in priority order (thin categories + arid-skarn replication first).
SITES = [
    # 1a Uranium/Energy
    "shirley_basin", "crooks_gap", "ambrosia_lake", "uravan", "henry_mountains",
    "lisbon_valley", "karnes_county", "arizona_strip",
    # 1b Mafic/Ultramafic
    "duluth_partridge_river", "duluth_south_kawishiwi", "webster_addie",
    "josephine_ophiolite", "new_idria", "state_line_chromite",
    # 1d Arid skarn/CRD replication (predicted POSITIVE — run early)
    "magdalena_kelly", "lake_valley", "organ_mountains", "victorio", "san_francisco_frisco",
    "shafter", "cherry_creek", "providence_mountains", "eagle_mountain", "contact_district",
    # 1c Sediment-hosted
    "getchell", "cortez", "alligator_ridge", "long_canyon", "tri_state",
    "upper_mississippi_valley", "mascot_jefferson_city", "austinville_ivanhoe", "black_butte",
    # 2a VMS
    "west_shasta", "east_shasta", "big_mike", "iron_king", "holden", "gossan_lead", "ore_knob",
    # 2b Humid skarn controls
    "cornwall_pa", "french_creek_pa", "dillsburg_pa", "willsboro_lewis_ny",
    "tilly_foster_ny", "snoqualmie_wa",
    # 3a Epithermal (exposed cap)
    "summitville", "round_mountain", "republic", "jarbidge",
    # 3b Alkaline/Carbonatite
    "wet_mountains", "round_top", "cornudas", "lemitar",
    # 3c Porphyry
    "safford_az", "san_manuel_az", "tyrone_nm", "questa_nm",
]


def _log(msg: str) -> None:
    line = f"{dt.datetime.now().isoformat(timespec='seconds')}  {msg}"
    print(line, flush=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def _is_done(site_id: str) -> bool:
    return (REPO / "results" / f"{site_id}_summary.csv").exists()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", nargs=2, type=int, metavar=("START", "END"))
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--force", action="store_true", help="re-run even if *_summary.csv exists")
    args = ap.parse_args()

    todo = SITES[args.batch[0]:args.batch[1]] if args.batch else SITES

    if args.list:
        for i, s in enumerate(SITES):
            print(f"{i:3d}  {'DONE' if _is_done(s) else '    '}  {s}")
        return 0

    _log(f"=== run start: {len(todo)} sites in scope "
         f"({'batch %d:%d' % tuple(args.batch) if args.batch else 'all'}) ===")
    ok, failed, skipped = [], [], []

    for i, site_id in enumerate(todo, 1):
        if _is_done(site_id) and not args.force:
            _log(f"[{i}/{len(todo)}] SKIP (summary exists): {site_id}")
            skipped.append(site_id)
            continue
        _log(f"[{i}/{len(todo)}] RUN: {site_id}")
        cmd = [sys.executable, "-m", "critical_minerals_aster", "run",
               "--site", site_id, "--mosaic"]
        # start_new_session=True so a hung download (which a plain subprocess.run
        # timeout does NOT reliably kill -- see the gossan_lead 9-hour hang) can be
        # taken down by its whole process group.
        proc = subprocess.Popen(cmd, cwd=REPO, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
        try:
            out, err = proc.communicate(timeout=HARD_TIMEOUT)
            if proc.returncode == 0 and _is_done(site_id):
                _log(f"[{i}/{len(todo)}] OK: {site_id}")
                ok.append(site_id)
            else:
                tail = (err or out or "").strip().splitlines()[-3:]
                _log(f"[{i}/{len(todo)}] FAIL rc={proc.returncode}: {site_id} :: {' | '.join(tail)}")
                failed.append(site_id)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            except ProcessLookupError:
                pass
            try:
                proc.communicate(timeout=30)
            except Exception:  # noqa: BLE001
                pass
            _log(f"[{i}/{len(todo)}] FAIL timeout({HARD_TIMEOUT}s, killed group): {site_id}")
            failed.append(site_id)
        except Exception as e:  # noqa: BLE001 - never let one site kill the batch
            _log(f"[{i}/{len(todo)}] FAIL exc: {site_id} :: {e!r}")
            failed.append(site_id)

    _log(f"=== run done: {len(ok)} ok, {len(failed)} failed, {len(skipped)} skipped ===")
    if failed:
        _log(f"FAILED: {failed}")
    _log(f"OK: {ok}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
