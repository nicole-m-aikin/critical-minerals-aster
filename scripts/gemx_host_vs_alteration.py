"""Batch TIR strong-zone vs SWIR band-depth correspondence (109-site capable).

Default site list = phase-4 primary FDR sites (19).  Sites without AVIRIS-5
coverage are skipped with a ``no_coverage`` status row.  Reuses download /
mosaic / band-depth logic from ``gemx_swir_pilot`` and overlap metrics from
``gemx_swir_correspondence``.

Does not download cubes unless ``--run-pilot`` is set and GeoTIFFs are absent.
Prefer ``--from-tif`` / existing mosaics for resume-friendly runs.

Usage:
    conda run -n aster-minerals python scripts/gemx_host_vs_alteration.py
    conda run -n aster-minerals python scripts/gemx_host_vs_alteration.py --primary-fdr
    conda run -n aster-minerals python scripts/gemx_host_vs_alteration.py --all-covered
    conda run -n aster-minerals python scripts/gemx_host_vs_alteration.py --sites goldfield marysvale bodie
    conda run -n aster-minerals python scripts/gemx_host_vs_alteration.py --site eureka --from-tif

Verdict heuristic (Al-OH primary; see ``gemx_site_select.correspondence_verdict``):

- no_coverage — no AVIRIS-5 coverage
- insufficient — missing TIFFs/TIR, tiny overlap, or mixed band
  (Jaccard < 0.15 and 1.1 ≤ lift < 1.3)
- overlapping — Jaccard ≥ 0.15 or lift ≥ 1.3 (TIR coincides with clay)
- complementary — Jaccard < 0.15 and lift < 1.1 (TIR not SWIR-clay)
- error — exception during compute
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from gemx_site_select import (  # noqa: E402
    add_site_select_args,
    sites_from_namespace,
)
from gemx_swir_correspondence import run_site_correspondence  # noqa: E402
from gemx_swir_pilot import run_site as run_pilot_site  # noqa: E402

CSV_FIELDS = [
    "site_id",
    "status",
    "swir_coverage",
    "jaccard_aloh",
    "lift_aloh",
    "jaccard_kaolinite",
    "jaccard_alunite",
    "n_tir_strong",
    "n_swir_high",
    "verdict",
]

OUT_CSV = ROOT / "results" / "gemx_host_vs_alteration.csv"


def write_csv(rows: list[dict], path: Path = OUT_CSV) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return path


def run_batch(
    site_ids: list[str],
    *,
    run_pilot: bool = False,
    from_tif: bool = True,
    write_figure: bool = False,
) -> list[dict]:
    """Correspondence (and optional pilot) for each site; return CSV rows."""
    rows: list[dict] = []
    for sid in site_ids:
        if run_pilot:
            pilot_row = run_pilot_site(sid, from_tif=from_tif)
            if pilot_row["status"] == "no_coverage":
                rows.append(
                    {
                        "site_id": sid,
                        "status": "no_coverage",
                        "swir_coverage": pilot_row.get("swir_coverage", 0.0),
                        "jaccard_aloh": None,
                        "lift_aloh": None,
                        "jaccard_kaolinite": None,
                        "jaccard_alunite": None,
                        "n_tir_strong": None,
                        "n_swir_high": None,
                        "verdict": "no_coverage",
                    }
                )
                continue
            if pilot_row["status"] == "error":
                rows.append(
                    {
                        "site_id": sid,
                        "status": "error",
                        "swir_coverage": pilot_row.get("swir_coverage"),
                        "jaccard_aloh": None,
                        "lift_aloh": None,
                        "jaccard_kaolinite": None,
                        "jaccard_alunite": None,
                        "n_tir_strong": None,
                        "n_swir_high": None,
                        "verdict": "error",
                    }
                )
                continue
        row = run_site_correspondence(
            sid,
            write_figure=write_figure,
            write_json=True,
        )
        rows.append(row)
    return rows


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Batch TIR vs SWIR host-vs-alteration correspondence; "
            "default = primary FDR sites"
        ),
    )
    add_site_select_args(p)
    p.add_argument(
        "--run-pilot",
        action="store_true",
        help="Run band-depth mosaic first (downloads cubes if TIFFs missing)",
    )
    p.add_argument(
        "--from-tif",
        action="store_true",
        default=True,
        help="Skip cube download when GeoTIFFs exist (default: on)",
    )
    p.add_argument(
        "--download",
        action="store_true",
        help="Allow cube download when TIFFs are missing (implies --run-pilot)",
    )
    p.add_argument(
        "--figures",
        action="store_true",
        help="Write per-site correspondence PNGs",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=OUT_CSV,
        help=f"Output CSV (default: {OUT_CSV.relative_to(ROOT)})",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    site_ids = sites_from_namespace(args, ROOT)
    run_pilot = bool(args.run_pilot or args.download)
    from_tif = not args.download
    print(f"Sites ({len(site_ids)}): {', '.join(site_ids)}", flush=True)
    rows = run_batch(
        site_ids,
        run_pilot=run_pilot,
        from_tif=from_tif,
        write_figure=args.figures,
    )
    path = write_csv(rows, args.out)
    n_ok = sum(1 for r in rows if r.get("status") == "ok")
    n_skip = sum(1 for r in rows if r.get("verdict") == "no_coverage")
    print(
        f"Wrote {path.relative_to(ROOT)}  "
        f"({n_ok} ok, {n_skip} no_coverage, {len(rows)} total)",
        flush=True,
    )


if __name__ == "__main__":
    main()
