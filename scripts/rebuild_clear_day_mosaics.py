#!/usr/bin/env python3
"""Plan and (optionally) rebuild site mosaics from clear-sky, daytime ASTER granules.

Plan mode (default, CMR metadata only, no login):
    python scripts/rebuild_clear_day_mosaics.py --all
    -> results/clear_day_mosaic_plan.csv

Execute mode (EarthData download via ~/.netrc, then per-site pipeline run):
    python scripts/rebuild_clear_day_mosaics.py --sites bisbee tombstone --execute
    python scripts/rebuild_clear_day_mosaics.py --changed --execute --skip-figures

Selection is ``spectral.select_mosaic_granules``: clear (<= 30% scene cloud) day
scenes first, never mixing day and night, greedily covering the bbox with as few
granules as possible. Existing ASTER files are moved to
``{aster_dir}_archive_pre_clearday/`` (never deleted). Sites whose YAML pins a
granule are only rebuilt with ``--unpin``, which sets ``granule_id: null`` and
keeps the old id in a comment.

After executing, rerun the analysis chain (regenerate_site_summaries ->
site_specific_null_significance -> phase3..phase8 -> fig01..fig06).
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
import warnings
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from critical_minerals_aster.config import list_site_ids, load_site_by_id  # noqa: E402
from critical_minerals_aster.paths import site_paths_for  # noqa: E402
from critical_minerals_aster.pipeline import (  # noqa: E402
    ASTER_SEARCH_COUNT,
    download_and_mosaic_aster,
    plan_mosaic_selection,
    run_site,
)
from critical_minerals_aster.config import search_bbox  # noqa: E402
from critical_minerals_aster.spectral import (  # noqa: E402
    extract_granule_id,
    granule_cloud_cover,
    granule_day_night,
)

warnings.filterwarnings("ignore", category=FutureWarning, module="earthaccess")
warnings.filterwarnings("ignore", message=".*DataGranule.size.*")

RESULTS = REPO / "results"
PLAN_CSV = RESULTS / "clear_day_mosaic_plan.csv"
ARCHIVE_SUFFIX = "_archive_pre_clearday"


def _current_state() -> pd.DataFrame:
    audit_path = RESULTS / "granule_cloud_audit.csv"
    if not audit_path.is_file():
        return pd.DataFrame(columns=["site_id"])
    audit = pd.read_csv(audit_path)
    grp = audit.groupby("site_id")
    return pd.DataFrame({
        "current_granules": grp["granule_id"].apply(lambda s: ";".join(sorted(s.dropna()))),
        "current_n_granules": grp["granule_id"].nunique(),
        "current_blend_cloud_pct": grp["site_blend_cloud_pct"].first(),
        "current_night_fraction": grp["day_night"].apply(lambda s: float((s == "Night").mean())),
    }).reset_index()


def _significance() -> pd.DataFrame:
    path = RESULTS / "cloud_sensitivity_sites.csv"
    if not path.is_file():
        return pd.DataFrame(columns=["site_id"])
    return pd.read_csv(path)[["site_id", "sig_fdr_binomial", "robust"]]


def plan_site(site_id: str) -> dict:
    import earthaccess

    site = load_site_by_id(site_id, REPO / "sites")
    results = earthaccess.search_data(
        short_name="AST_L1T",
        bounding_box=search_bbox(site),
        temporal=(site.temporal_start, site.temporal_end),
        count=ASTER_SEARCH_COUNT,
    )
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        sel = plan_mosaic_selection(site, search_results=results)

    pinned_state = None
    if site.granule_id:
        for g in results:
            if site.granule_id in extract_granule_id(g):
                pinned_state = f"{granule_day_night(g)}/{granule_cloud_cover(g)}"
                break

    sizes = []
    for g in sel.granules:
        try:
            sizes.append(float(g.size()))
        except Exception:
            pass
    clouds = [r["cloud_pct"] for r in sel.records if r["cloud_pct"] is not None]
    n_day_clear = sum(
        1 for g in results
        if granule_day_night(g) == "Day"
        and (granule_cloud_cover(g) is None or granule_cloud_cover(g) <= 30)
    )
    return {
        "site_id": site_id,
        "pinned_granule": site.granule_id or "",
        "pinned_day_night_cloud": pinned_state or "",
        "n_search_results": len(results),
        "n_day_clear_available": n_day_clear,
        "new_tier": sel.tier,
        "new_n_granules": len(sel.granules),
        "new_union_coverage": round(sel.union_coverage, 3),
        "new_max_cloud_pct": max(clouds) if clouds else None,
        "new_mean_cloud_pct": round(sum(clouds) / len(clouds), 1) if clouds else None,
        "est_download_mb": round(sum(sizes), 1),
        "new_granules": ";".join(r["granule_id"] for r in sel.records),
    }


def write_plan(site_ids: list[str]) -> pd.DataFrame:
    rows = []
    for i, sid in enumerate(site_ids, 1):
        try:
            rows.append(plan_site(sid))
            r = rows[-1]
            print(f"[{i}/{len(site_ids)}] {sid}: {r['new_tier']} x{r['new_n_granules']} "
                  f"cov={r['new_union_coverage']:.2f} cloud_max={r['new_max_cloud_pct']}")
        except Exception as exc:
            print(f"[{i}/{len(site_ids)}] {sid}: plan failed ({exc})", file=sys.stderr)
            rows.append({"site_id": sid, "new_tier": f"error: {exc}"})
    plan = pd.DataFrame(rows)
    plan = plan.merge(_current_state(), on="site_id", how="left")
    plan = plan.merge(_significance(), on="site_id", how="left")
    plan["changed"] = [
        set(str(n).split(";")) != set(str(c).split(";"))
        for n, c in zip(plan.get("new_granules", ""), plan.get("current_granules", ""))
    ]
    RESULTS.mkdir(exist_ok=True)
    tmp = PLAN_CSV.with_suffix(".csv.tmp")
    plan.to_csv(tmp, index=False)
    tmp.replace(PLAN_CSV)
    print(f"Wrote {PLAN_CSV.relative_to(REPO)}")
    return plan


def _archive_aster_dir(aster_dir: Path) -> Path | None:
    files = [p for p in aster_dir.iterdir() if p.is_file() and p.suffix in (".tif", ".json")]
    files = [p for p in files if "AST_L1T" in p.name or "_mosaic" in p.name]
    if not files:
        return None
    archive = aster_dir.parent / f"{aster_dir.name}{ARCHIVE_SUFFIX}"
    archive.mkdir(exist_ok=True)
    for p in files:
        dest = archive / p.name
        if dest.exists():
            p.unlink()
        else:
            shutil.move(str(p), dest)
    return archive


def _unpin_yaml(site_id: str) -> None:
    yml = REPO / "sites" / f"{site_id}.yaml"
    text = yml.read_text()
    new = re.sub(
        r'^granule_id:\s*"?(AST_L1T_[0-9_]+)"?\s*$',
        r'granule_id: null  # was "\1"; unpinned for clear-sky day mosaic',
        text, count=1, flags=re.M,
    )
    if new != text:
        yml.write_text(new)


def execute(site_ids: list[str], *, unpin: bool, skip_figures: bool) -> None:
    import earthaccess

    earthaccess.login(strategy="netrc")
    for sid in site_ids:
        site = load_site_by_id(sid, REPO / "sites")
        if site.granule_id:
            if not unpin:
                print(f"{sid}: pinned granule {site.granule_id}; skipped (use --unpin)")
                continue
            _unpin_yaml(sid)
            site = load_site_by_id(sid, REPO / "sites")
        paths = site_paths_for(site, REPO)
        if paths.aster_dir.is_dir():
            archive = _archive_aster_dir(paths.aster_dir)
            if archive:
                print(f"{sid}: archived old ASTER files -> {archive.relative_to(REPO)}")
        print(f"{sid}: downloading clear-sky day granules ...", flush=True)
        try:
            gid = download_and_mosaic_aster(site, paths, interactive_login=False)
        except Exception as exc:
            print(f"{sid}: download/mosaic FAILED ({exc})", file=sys.stderr)
            continue
        print(f"{sid}: mosaic {gid}; running pipeline", flush=True)
        try:
            run_site(site, REPO, download=False, skip_figures=skip_figures)
        except Exception as exc:
            print(f"{sid}: run_site FAILED ({exc})", file=sys.stderr)
            continue
        print(f"{sid}: OK", flush=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--all", action="store_true", help="All sites in sites/index.yaml")
    g.add_argument("--sites", nargs="+")
    g.add_argument("--changed", action="store_true",
                   help="Sites whose planned granule set differs from the current one (needs a plan CSV)")
    ap.add_argument("--execute", action="store_true", help="Download + rebuild (default: plan only)")
    ap.add_argument("--unpin", action="store_true", help="Unpin YAML granule_id for pinned sites")
    ap.add_argument("--skip-figures", action="store_true")
    args = ap.parse_args(argv)

    if args.changed:
        if not PLAN_CSV.is_file():
            print("No plan CSV; run with --all first", file=sys.stderr)
            return 1
        plan = pd.read_csv(PLAN_CSV)
        site_ids = plan.loc[plan["changed"].astype(bool), "site_id"].tolist()
    elif args.all:
        site_ids = list_site_ids(REPO / "sites")
    else:
        site_ids = args.sites

    if not args.execute:
        if not args.changed:
            write_plan(site_ids)
        else:
            print("\n".join(site_ids))
        return 0

    execute(site_ids, unpin=args.unpin, skip_figures=args.skip_figures)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
