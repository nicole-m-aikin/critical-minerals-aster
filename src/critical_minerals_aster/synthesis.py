"""Aggregate per-site summaries into a national comparison table."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_site_summaries(results_dir: Path) -> pd.DataFrame:
    """Load all *_summary.csv files; return site-level rows only."""
    results_dir = Path(results_dir)
    frames: list[pd.DataFrame] = []
    for path in sorted(p for p in results_dir.glob("*_summary.csv") 
                   if "national" not in p.name):
        df = pd.read_csv(path)
        site_rows = df[df["row_type"] == "site"] if "row_type" in df.columns else df
        frames.append(site_rows)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def write_national_summary(results_dir: Path) -> Path:
    """Write national_summary.csv and national_summary.parquet under results/."""
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    national = load_site_summaries(results_dir)
    csv_path = results_dir / "national_summary.csv"
    national.to_csv(csv_path, index=False)
    parquet_path = results_dir / "national_summary.parquet"
    try:
        national.to_parquet(parquet_path, index=False)
    except ImportError:
        parquet_path = None
    return csv_path
