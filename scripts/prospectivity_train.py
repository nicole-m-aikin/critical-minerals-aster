#!/usr/bin/env python3
"""Train tabular prospectivity MVP under site- and era-blocked validation."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.prospectivity.adapters import load_aster_thermal_scene
from critical_minerals_aster.prospectivity.features import (
    build_site_feature_table,
    load_feature_store,
    model_matrix,
)
from critical_minerals_aster.prospectivity.models import (
    ablation_columns,
    fit_model,
    predict_proba,
    run_era_block,
    run_site_cv,
)


def _write_rank_map(site_id: str, scores_df: pd.DataFrame, out_path: Path) -> None:
    """Write a sparse ranked GeoTIFF from scored sample cells (nearest fill)."""
    site = load_site_by_id(site_id, ROOT / "sites")
    paths = site_paths_for(site, ROOT)
    try:
        thermal = load_aster_thermal_scene(site, paths)
    except FileNotFoundError:
        return
    h, w = thermal.qi.shape
    grid = np.full((h, w), np.nan, dtype=np.float32)
    for r, c, s in zip(scores_df["row"], scores_df["col"], scores_df["score_pred"]):
        if 0 <= int(r) < h and 0 <= int(c) < w:
            grid[int(r), int(c)] = float(s)
    # Fill remaining valid TIR cells with median score for a continuous map
    valid = thermal.qa_mask
    med = float(np.nanmedian(grid)) if np.isfinite(grid).any() else 0.0
    out = np.where(valid, np.where(np.isfinite(grid), grid, med), np.nan).astype(np.float32)
    profile = {
        "driver": "GTiff",
        "height": h,
        "width": w,
        "count": 1,
        "dtype": "float32",
        "crs": thermal.crs,
        "transform": thermal.transform,
        "compress": "deflate",
        "nodata": np.nan,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(out, 1)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--model", choices=("gbm", "logistic"), default="gbm")
    p.add_argument("--n-splits", type=int, default=5)
    p.add_argument("--ablations", default="all,tir,hsi")
    p.add_argument("--max-sites", type=int, default=None)
    p.add_argument("--write-maps", action="store_true")
    args = p.parse_args()

    df = load_feature_store(ROOT)
    if df.empty:
        print("Feature store empty — building a small cohort from feature-ready…")
        fr = ROOT / "results" / "gemx_census_feature_ready.csv"
        if fr.is_file():
            ready = pd.read_csv(fr)
            sids = ready.loc[ready["feature_ready"].astype(bool), "site_id"].tolist()
        else:
            sids = ["mcdermitt", "goldfield", "marysvale", "bodie", "morenci"]
        if args.max_sites:
            sids = sids[: args.max_sites]
        frames = []
        for sid in sids:
            tab = build_site_feature_table(sid, repo_root=ROOT)
            if tab is not None and len(tab):
                frames.append(tab)
        if not frames:
            raise SystemExit("No feature tables could be built")
        df = pd.concat(frames, ignore_index=True)

    if args.max_sites:
        keep = df["site_id"].drop_duplicates().tolist()[: args.max_sites]
        df = df[df["site_id"].isin(keep)].copy()

    print(f"Training on {df['site_id'].nunique()} sites, {len(df)} cells, "
          f"pos={int(df['label'].sum())}")

    results_dir = ROOT / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    site_cv_rows = []
    summaries = []

    for ab in [a.strip() for a in args.ablations.split(",") if a.strip()]:
        per_site, summary = run_site_cv(
            df, n_splits=args.n_splits, kind=args.model, ablation=ab
        )
        site_cv_rows.append(per_site)
        summaries.append(summary)
        print(f"site-CV [{ab}]:", json.dumps(summary, indent=2))

    site_cv = pd.concat(site_cv_rows, ignore_index=True) if site_cv_rows else pd.DataFrame()
    site_cv_path = results_dir / "prospectivity_site_cv.csv"
    site_cv.to_csv(site_cv_path, index=False)
    print(f"Wrote {site_cv_path.relative_to(ROOT)}")

    era = {}
    for ab in [a.strip() for a in args.ablations.split(",") if a.strip()]:
        era[ab] = run_era_block(df, kind=args.model, ablation=ab)
        print(f"era-block [{ab}] status={era[ab].get('status')} "
              f"n_train_pos={era[ab].get('n_train_pos')} n_test_pos={era[ab].get('n_test_pos')}")

    era_path = results_dir / "prospectivity_era_block.json"
    era_path.write_text(json.dumps(era, indent=2, default=str))
    print(f"Wrote {era_path.relative_to(ROOT)}")

    # Fit final all-data model (for maps / importance only — not a claim)
    cols = ablation_columns("all")
    X, y, _ = model_matrix(df, cols=cols)
    if len(np.unique(y)) >= 2:
        fit = fit_model(X, y, kind=args.model, feature_names=cols)
        df = df.copy()
        df["score_pred"] = predict_proba(fit, X)
        if args.write_maps:
            map_dir = ROOT / "data" / "prospectivity" / "maps"
            for sid, sub in df.groupby("site_id"):
                _write_rank_map(str(sid), sub, map_dir / f"{sid}_prospectivity.tif")
            print(f"Wrote ranked maps under data/prospectivity/maps/")

    summary_path = results_dir / "prospectivity_site_cv_summary.json"
    summary_path.write_text(json.dumps(summaries, indent=2, default=str))
    print(f"Wrote {summary_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
