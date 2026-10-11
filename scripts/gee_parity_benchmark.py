#!/usr/bin/env python
"""Local rasterio vs Google Earth Engine parity benchmark.

Default (offline, no GEE credentials):
    For each site in the 8-site panel, run the local TIR workflow and
    four discrepancy probes that isolate the implementation choices a
    cloud-native port actually has to make:

        clip_order             — percentiles on the unclipped scene
        histogram_percentile   — GEE Reducer.percentile analogue
        radiance_ucc           — DN vs (DN-1)*UCC radiance ratios
        reproject_4326         — GEE's default reduceRegion CRS

    Writes:
        results/gee_parity_agreement.csv
        results/gee_parity_agreement.json
        results/gee_parity_recipes/{site}.js

Live Earth Engine (optional):
    python scripts/gee_parity_benchmark.py --gee --project YOUR_GCP_PROJECT

    Adds a ``gee_live`` row per site with reduceRegion percentiles and
    the selected ``ASTER/AST_L1T_003`` image id.  Pixel-array export is
    intentionally not the default (slow, size-capped); the probes already
    quantify the map-level effect of each GEE default.

Usage:
    conda run -n aster-minerals python scripts/gee_parity_benchmark.py
    conda run -n aster-minerals python scripts/gee_parity_benchmark.py \\
        --sites goldfield bisbee
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.gee_parity import (
    DEFAULT_PARITY_SITES,
    SITE_ROLES,
    gee_recipe_javascript,
    recipe_for_site,
    run_gee_workflow,
    run_site_probes,
    summarize_probe_effect,
)
from critical_minerals_aster.paths import site_paths_for

RESULTS_DIR = ROOT / "results"
RECIPES_DIR = RESULTS_DIR / "gee_parity_recipes"


def _write_recipes(site_ids: list[str]) -> None:
    RECIPES_DIR.mkdir(parents=True, exist_ok=True)
    for site_id in site_ids:
        site = load_site_by_id(site_id, ROOT / "sites")
        recipe = recipe_for_site(site)
        js_path = RECIPES_DIR / f"{site_id}.js"
        js_path.write_text(gee_recipe_javascript(recipe))
        print(f"  recipe  {js_path.relative_to(ROOT)}")


def _run_offline(site_ids: list[str]) -> pd.DataFrame:
    rows = []
    for site_id in site_ids:
        site = load_site_by_id(site_id, ROOT / "sites")
        paths = site_paths_for(site, ROOT)
        role = SITE_ROLES.get(site_id, "")
        if not paths.aster_dir.is_dir():
            print(f"  skip    {site_id} ({role}) — no ASTER dir", file=sys.stderr)
            continue
        try:
            local, site_rows = run_site_probes(site, paths)
        except FileNotFoundError as exc:
            print(f"  skip    {site_id} ({role}) — {exc}", file=sys.stderr)
            continue
        print(
            f"  local   {site_id:<16} {role:<22} "
            f"n={local.n_valid:,}  strong={local.strong_frac:.3f}  "
            f"{local.notes}"
        )
        for row in site_rows:
            print(
                f"          {row.probe:<32} "
                f"Jaccard={row.jaccard_strong:.3f}  "
                f"κ={row.kappa:.3f}  "
                f"Δstrong={row.strong_frac_delta:+.3f}"
            )
            rows.append(row.to_dict())
    return pd.DataFrame(rows)


def _run_gee(site_ids: list[str], project: str | None) -> list[dict]:
    live = []
    for site_id in site_ids:
        site = load_site_by_id(site_id, ROOT / "sites")
        recipe = recipe_for_site(site)
        # Prefer the local raster CRS so reduceRegion does not silently
        # resample into EPSG:4326 — that is a probe, not the live default.
        paths = site_paths_for(site, ROOT)
        try:
            from critical_minerals_aster.gee_parity import run_local_workflow

            local = run_local_workflow(site, paths)
            if local.crs:
                from dataclasses import replace

                recipe = replace(recipe, crs=local.crs)
        except Exception:
            pass
        print(f"  gee     {site_id} …")
        result = run_gee_workflow(site, recipe, project=project)
        live.append(result)
        image_id = (result.get("gee") or {}).get("image_id")
        print(f"          image_id={image_id}")
    return live


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sites",
        nargs="*",
        default=list(DEFAULT_PARITY_SITES),
        help="Site ids (default: 8-site representative panel)",
    )
    parser.add_argument(
        "--gee",
        action="store_true",
        help="Also run the live Earth Engine recipe (needs earthengine-api + auth)",
    )
    parser.add_argument(
        "--project",
        default=None,
        help="GCP project for ee.Initialize (required by current GEE)",
    )
    args = parser.parse_args(argv)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"GEE parity panel: {', '.join(args.sites)}")
    _write_recipes(args.sites)

    df = _run_offline(args.sites)
    if df.empty:
        print("No local ASTER rasters found for the requested sites.", file=sys.stderr)
        return 1

    csv_path = RESULTS_DIR / "gee_parity_agreement.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nWrote {csv_path.relative_to(ROOT)}  ({len(df)} rows)")

    from critical_minerals_aster.gee_parity import ParityRow

    rows = [
        ParityRow(**{k: rec[k] for k in ParityRow.__dataclass_fields__})
        for rec in df.to_dict("records")
    ]
    summary = summarize_probe_effect(rows)

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "panel": args.sites,
        "n_rows": int(len(df)),
        "summary_by_probe": summary,
        "gee_live": None,
        "notes": (
            "Offline probes isolate implementation choices. "
            "gee_live is populated only with --gee. "
            "GEE hosts ASTER/AST_L1T_003; local rasters are LP DAAC v004."
        ),
    }

    if args.gee:
        payload["gee_live"] = _run_gee(args.sites, args.project)

    json_path = RESULTS_DIR / "gee_parity_agreement.json"
    json_path.write_text(json.dumps(payload, indent=2, default=str))
    print(f"Wrote {json_path.relative_to(ROOT)}")

    print("\nProbe effect (mean across sites):")
    print(f"{'probe':<32} {'n':>4} {'Jaccard':>8} {'|Δ strong|':>10} {'κ':>8}")
    for probe, stats in sorted(summary.items()):
        print(
            f"{probe:<32} {int(stats['n_sites']):>4} "
            f"{stats['mean_jaccard_strong']:>8.3f} "
            f"{stats['mean_abs_strong_frac_delta']:>10.3f} "
            f"{stats['mean_kappa']:>8.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
