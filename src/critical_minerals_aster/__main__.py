"""CLI: python -m critical_minerals_aster run --site mcdermitt"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from critical_minerals_aster.config import list_site_ids, load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.pipeline import (
    download_and_mosaic_aster,
    run_batch,
    run_batch_parallel,
    run_site,
)
from critical_minerals_aster.synthesis import write_national_summary


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _should_skip(site_id: str, repo_root: Path) -> bool:
    """Return True if site outputs already exist and can be skipped."""
    from critical_minerals_aster.config import load_site_by_id

    try:
        site = load_site_by_id(site_id, repo_root / "sites")
        paths = site_paths_for(site, repo_root)
        overlay = paths.figures_dir / "03_deposit_overlay.png"
        prov = paths.site_provenance_json
        return overlay.exists() and prov.exists()
    except Exception:
        return False


def cmd_run(args: argparse.Namespace) -> int:
    repo = Path(args.repo_root) if args.repo_root else _repo_root()
    skip_existing = getattr(args, "skip_existing", False)

    if skip_existing and _should_skip(args.site, repo):
        print(
            f"Skipping {args.site} (outputs exist, use --force to regenerate)",
            file=sys.stderr,
        )
        return 0

    site = load_site_by_id(args.site, repo / "sites")
    if getattr(args, "mosaic", False):
        paths = site_paths_for(site, repo)
        download_and_mosaic_aster(site, paths)
        # Re-run without download so run_site picks up the freshly built mosaic.
        run_site(
            site, repo,
            download=False,
            skip_figures=args.skip_figures,
            skip_existing=False,
        )
    else:
        run_site(
            site,
            repo,
            download=args.download,
            skip_figures=args.skip_figures,
            skip_existing=skip_existing,
        )
    print(f"Finished site {args.site}; outputs under {repo / 'results'}")
    return 0


def cmd_run_batch(args: argparse.Namespace) -> int:
    repo = Path(args.repo_root) if args.repo_root else _repo_root()
    sites_dir = repo / "sites"
    site_ids = list_site_ids(sites_dir) if args.all_sites else args.sites
    if not site_ids:
        print("No sites specified.", file=sys.stderr)
        return 1

    skip_existing = getattr(args, "skip_existing", False)
    workers = getattr(args, "workers", 1)

    if workers > 1 and len(site_ids) > 1:
        run_batch_parallel(
            site_ids,
            repo,
            workers=workers,
            download=args.download,
            skip_figures=args.skip_figures,
            skip_existing=skip_existing,
        )
    else:
        run_batch(
            site_ids,
            repo,
            download=args.download,
            skip_figures=args.skip_figures,
            skip_existing=skip_existing,
        )
    write_national_summary(repo / "results")
    print(f"Batch complete; national summary in {repo / 'results'}")
    return 0


def cmd_synthesize(args: argparse.Namespace) -> int:
    repo = Path(args.repo_root) if args.repo_root else _repo_root()
    path = write_national_summary(repo / "results")
    print(f"Wrote {path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="critical_minerals_aster",
        description="ASTER TIR alteration pipeline (multi-site)",
    )
    parser.add_argument(
        "--repo-root",
        help="Repository root (default: parent of src/)",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="Process one study site")
    p_run.add_argument("--site", required=True, help="Site id (sites/{id}.yaml)")
    p_run.add_argument(
        "--download",
        action="store_true",
        help="Download ASTER from EarthData before processing",
    )
    p_run.add_argument(
        "--mosaic",
        action="store_true",
        help="Download ALL covering ASTER granules, merge per-band, then process",
    )
    p_run.add_argument("--skip-figures", action="store_true")
    p_run.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip site if 03_deposit_overlay.png and provenance JSON already exist",
    )
    p_run.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help="Number of parallel workers (default: 1; single-site run ignores this)",
    )
    p_run.set_defaults(func=cmd_run)

    p_batch = sub.add_parser("run-batch", help="Process multiple sites")
    p_batch.add_argument(
        "--sites",
        nargs="*",
        default=[],
        help="Site ids (default: all in sites/index.yaml with --all-sites)",
    )
    p_batch.add_argument(
        "--all-sites",
        action="store_true",
        help="Run every site listed in sites/index.yaml",
    )
    p_batch.add_argument("--download", action="store_true")
    p_batch.add_argument("--skip-figures", action="store_true")
    p_batch.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip sites whose outputs already exist",
    )
    p_batch.add_argument(
        "--workers",
        type=int,
        default=1,
        metavar="N",
        help="Number of parallel worker processes (default: 1)",
    )
    p_batch.set_defaults(func=cmd_run_batch)

    p_syn = sub.add_parser("synthesize", help="Aggregate results/*_summary.csv")
    p_syn.set_defaults(func=cmd_synthesize)

    args = parser.parse_args(argv)
    return args.func(args)


def cli() -> None:
    raise SystemExit(main())


if __name__ == "__main__":
    cli()
