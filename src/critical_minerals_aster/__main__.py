"""CLI: python -m critical_minerals_aster run --site mcdermitt"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from critical_minerals_aster.config import list_site_ids, load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.pipeline import (
    download_and_mosaic_aster,
    fig03_outputs_current,
    run_batch,
    run_batch_parallel,
    run_site,
)
from critical_minerals_aster.synthesis import write_national_summary


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _should_skip(site_id: str, repo_root: Path) -> bool:
    """Return True if dual-panel fig 03 and matching provenance exist."""
    from critical_minerals_aster.config import load_site_by_id

    try:
        site = load_site_by_id(site_id, repo_root / "sites")
        paths = site_paths_for(site, repo_root)
        return fig03_outputs_current(repo_root, paths)
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
    use_mosaic = getattr(args, "mosaic", False)

    if use_mosaic:
        for site_id in site_ids:
            site = load_site_by_id(site_id, sites_dir)
            paths = site_paths_for(site, repo)
            print(f"  [mosaic] {site_id}: downloading all covering granules …")
            try:
                download_and_mosaic_aster(site, paths)
            except Exception as exc:
                print(f"  [mosaic] {site_id}: failed ({exc}), skipping", file=sys.stderr)
        # Data is on disk; run without download flag so run_batch picks up mosaics.
        download = False
    else:
        download = args.download

    if workers > 1 and len(site_ids) > 1:
        run_batch_parallel(
            site_ids,
            repo,
            workers=workers,
            download=download,
            skip_figures=args.skip_figures,
            skip_existing=skip_existing,
        )
    else:
        run_batch(
            site_ids,
            repo,
            download=download,
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


def cmd_demo(args: argparse.Namespace) -> int:
    """Print a formatted summary from pre-built results — no API keys or rasters needed."""
    import textwrap

    repo = Path(args.repo_root) if args.repo_root else _repo_root()
    db_path = repo / "results" / "results.duckdb"
    sig_path = repo / "results" / "phase4_fdr_corrected_significance.csv"

    print("\n" + "=" * 70)
    print("  ASTER TIR Critical Minerals — Demo (pre-built results)")
    print("=" * 70)

    # --- Significance summary from CSV (no duckdb dep needed for demo) ---
    import pandas as pd

    if sig_path.exists():
        sig = pd.read_csv(sig_path)
        n_total = len(sig)
        binomial = sig["sig_fdr_binomial"].fillna(False).astype(bool)
        mannwhitney = sig["sig_fdr_mannwhitney"].fillna(False).astype(bool)
        sig_sites = sig[binomial].sort_values("enrichment", ascending=False)
        print(
            f"\n{int(binomial.sum())} of {n_total} sites significant on the primary test"
            " (BH-FDR corrected exact binomial, α = 0.05)"
        )
        print(
            f"{int(mannwhitney.sum())} of {n_total} on the threshold-free continuous-score test"
            f"  ·  union {int((binomial | mannwhitney).sum())}\n"
        )
        print(
            f"{'Site':<28} {'n crit':>8} {'Hits':>6} {'Hit rate':>10}"
            f" {'Null':>7} {'Enrich':>8} {'q (binom)':>10}"
        )
        print("-" * 82)
        for _, row in sig_sites.iterrows():
            name = str(row["site_name"])
            if len(name) > 26:
                name = name[:25] + "\u2026"
            print(
                f"  {name:<26} {int(row['n_crit']):>8}"
                f" {int(row['hits_crit']):>6} {row['hit_rate_crit_pct']:>9.1f}%"
                f" {row['p0_null'] * 100:>6.1f}% {row['enrichment']:>7.2f}x"
                f" {row['q_binomial']:>10.4f}"
            )
        print(
            "\nEnrichment is the observed hit rate divided by each site's own geometric"
            "\nnull (zone area / TIR footprint area). Of these, 13 also survive DBSCAN"
            "\nspatial declustering at every radius — see docs/results.md."
        )
    else:
        print(
            "\n  (phase4_fdr_corrected_significance.csv not found —"
            " run scripts/phase4_fdr_correction.py first)"
        )

    # --- DuckDB query if available ---
    if db_path.exists():
        try:
            import duckdb
            con = duckdb.connect(str(db_path), read_only=True)
            totals = con.execute(
                "SELECT COUNT(DISTINCT site_id) AS n_sites,"
                " SUM(n_deposits_bbox) AS total_deposits,"
                " SUM(n_deposits_in_zones) AS total_hits"
                " FROM site_summaries WHERE row_type = 'site'"
            ).fetchone()
            if totals:
                n_sites, total_dep, total_hit = totals
                print(f"\nNational totals ({n_sites} sites):")
                print(f"  Deposits evaluated : {int(total_dep):,}")
                print(f"  Deposits in zones  : {int(total_hit):,}")
                print(f"  Overall hit rate   : {total_hit/total_dep*100:.1f}%")
            con.close()
        except Exception as exc:
            print(f"\n  (DuckDB query skipped: {exc})")

    # --- Discovery bias summary (unconditional framing; see docs/results.md) ---
    bias_path = repo / "results" / "phase7_discovery_bias_pooled.csv"
    if bias_path.exists():
        bias = pd.read_csv(bias_path)
        primary = bias[bias["framing"].str.contains("PRIMARY", na=False)]
        if not primary.empty:
            n_sites = int(primary["n_sites"].iloc[0])
            print(
                f"\nDiscovery-bias stratification (dated deposits only, {n_sites} sites):"
            )
            for _, row in primary.iterrows():
                print(
                    f"  {row['cohort']:<10}: {int(row['hits'])}/{int(row['n'])} hits"
                    f" ({row['hr_pct']:.1f}% vs {row['null_hr_pct']:.1f}% null)"
                    f" — {row['enrichment']:.2f}x, p = {row['pooled_p']:.2f}"
                )
            print(
                "  No enrichment in either cohort. This test is uninformative rather than"
                "\n  exculpatory; the earlier circularity claim is retracted."
            )

    # --- Where to find figures ---
    gallery = repo / "figures" / "index.html"
    print(f"\nFigures:")
    if gallery.exists():
        print(f"  open {gallery}")
    else:
        print(f"  figures/index.html not found locally.")
    print(
        textwrap.dedent("""
        To run the full pipeline (requires EarthData account + ASTER rasters):
          python -m critical_minerals_aster run --site mcdermitt --mosaic
          python -m critical_minerals_aster synthesize

        To reproduce the corrected significance chain:
          conda run -n aster-minerals python scripts/site_specific_null_significance.py
          conda run -n aster-minerals python scripts/phase3_monte_carlo_and_continuous_score.py
          conda run -n aster-minerals python scripts/phase4_fdr_correction.py
          conda run -n aster-minerals python scripts/phase5_clustering_sensitivity.py
        """)
    )
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
    p_batch.add_argument(
        "--mosaic",
        action="store_true",
        help="Download ALL covering ASTER granules per site, merge per-band, then process",
    )
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

    p_demo = sub.add_parser(
        "demo",
        help="Print pre-built results summary (no API keys or rasters required)",
    )
    p_demo.set_defaults(func=cmd_demo)

    args = parser.parse_args(argv)
    return args.func(args)


def cli() -> None:
    raise SystemExit(main())


if __name__ == "__main__":
    cli()
