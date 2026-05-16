"""Aggregate per-site summaries into a national comparison table."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_site_summaries(
    results_dir: Path,
    row_types: list[str] | None = None,
) -> pd.DataFrame:
    """Load all *_summary.csv files, optionally filtered by row_type.

    Parameters
    ----------
    row_types:
        If given, keep only rows whose ``row_type`` is in this list.
        Defaults to ``["site"]`` to preserve the original behaviour for
        callers that expect one row per site.  Pass ``None`` to return all
        rows (site + commodity + earth_mri).
    """
    _ALL = object()  # sentinel: include every row_type
    _filter = _ALL if row_types is None else row_types
    results_dir = Path(results_dir)
    frames: list[pd.DataFrame] = []
    for path in sorted(
        p for p in results_dir.glob("*_summary.csv") if "national" not in p.name
    ):
        df = pd.read_csv(path)
        if "row_type" in df.columns and _filter is not _ALL:
            df = df[df["row_type"].isin(_filter)]  # type: ignore[arg-type]
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def write_national_summary(results_dir: Path) -> Path:
    """Write national_summary.csv and national_summary.parquet under results/."""
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    national = load_site_summaries(results_dir, row_types=None)
    csv_path = results_dir / "national_summary.csv"
    national.to_csv(csv_path, index=False)
    parquet_path = results_dir / "national_summary.parquet"
    try:
        national.to_parquet(parquet_path, index=False)
    except ImportError:
        parquet_path = None
    return csv_path
