#!/usr/bin/env python3
"""Build per-site prospectivity feature/label parquet store."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.prospectivity.features import build_feature_store


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--site", action="append", dest="sites", help="Limit to site id(s)")
    p.add_argument(
        "--feature-ready",
        type=Path,
        default=ROOT / "results" / "gemx_census_feature_ready.csv",
    )
    p.add_argument("--max-sites", type=int, default=None)
    p.add_argument("--all-sites", action="store_true", help="Ignore feature-ready filter")
    args = p.parse_args()

    fr = None if args.all_sites or args.sites else args.feature_ready
    df = build_feature_store(
        ROOT,
        site_ids=args.sites,
        feature_ready_csv=fr,
        max_sites=args.max_sites,
    )
    print(f"rows={len(df)} sites={df['site_id'].nunique() if len(df) else 0}")
    if len(df):
        print(df.groupby("site_id")["label"].agg(["count", "sum"]).head(20).to_string())


if __name__ == "__main__":
    main()
