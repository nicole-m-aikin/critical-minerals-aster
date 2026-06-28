"""
Generate pixel-level prospectivity probability maps for significant sites.

Applies a threshold-perturbation ensemble (3 scenarios: loose/nominal/tight)
to the pre-computed ratio mosaics and writes:

    data/sites/{id}/rasters/prospectivity_mean.tif  — P(strong anomaly)
    data/sites/{id}/rasters/prospectivity_std.tif   — classification uncertainty

See src/critical_minerals_aster/uncertainty.py for the full methodology.

Usage
-----
  conda run -n aster-minerals python scripts/generate_prospectivity_maps.py

  # Specific sites only:
  conda run -n aster-minerals python scripts/generate_prospectivity_maps.py \
      --sites bisbee eureka mcdermitt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.spectral import load_ratio_mosaic
from critical_minerals_aster.uncertainty import (
    DEFAULT_SCENARIOS,
    save_prospectivity_rasters,
    threshold_ensemble,
)


def _ratio_mosaic_id(site_id: str) -> str:
    """Resolve the granule_id prefix used for a site's ratio mosaic files.

    Prefers the canonical {site_id}_mosaic naming even when provenance records
    an individual granule ID — the pipeline always writes ratios as {site_id}_mosaic_ratio_*.tif.
    """
    aster_dir = ROOT / "data" / "sites" / site_id / "aster"
    mosaic_id = f"{site_id}_mosaic"
    if (aster_dir / f"{mosaic_id}_ratio_silica.tif").exists():
        return mosaic_id
    # Fall back to provenance granule_id for sites that pre-date the mosaic naming
    import json
    prov = ROOT / "results" / f"{site_id}_provenance.json"
    if prov.exists():
        with open(prov) as f:
            d = json.load(f)
        gid = d.get("granule_id", "")
        if gid:
            return gid
    return mosaic_id


def process_site(site_id: str) -> dict:
    aster_dir = ROOT / "data" / "sites" / site_id / "aster"
    out_dir = ROOT / "data" / "sites" / site_id / "rasters"

    mosaic_id = _ratio_mosaic_id(site_id)
    silica_path = aster_dir / f"{mosaic_id}_ratio_silica.tif"
    if not silica_path.exists():
        return {"site_id": site_id, "status": "skipped — ratio mosaic not found"}

    silica, carbonate, mafic, _b10, transform, crs = load_ratio_mosaic(aster_dir, mosaic_id)

    mean, std = threshold_ensemble(silica, carbonate, mafic, scenarios=DEFAULT_SCENARIOS)

    n_valid = int(np.isfinite(mean).sum())
    n_high = int((mean >= 1.0)[np.isfinite(mean)].sum())   # strong in all scenarios
    n_moderate = int(((mean >= 0.5) & (mean < 1.0))[np.isfinite(mean)].sum())
    n_marginal = int(((mean > 0.0) & (mean < 0.5))[np.isfinite(mean)].sum())
    mean_std = float(np.nanmean(std))

    mean_path, std_path = save_prospectivity_rasters(mean, std, transform, crs, out_dir)

    return {
        "site_id": site_id,
        "status": "ok",
        "n_valid_pixels": n_valid,
        "n_high_confidence": n_high,       # P=1.0
        "n_moderate_confidence": n_moderate,  # P=0.5–1.0
        "n_marginal": n_marginal,           # P=0.0–0.5
        "pct_high": round(n_high / n_valid * 100, 1) if n_valid else 0.0,
        "mean_uncertainty": round(mean_std, 4),
        "mean_path": str(mean_path.relative_to(ROOT)),
        "std_path": str(std_path.relative_to(ROOT)),
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sites", nargs="*",
        help="Site IDs to process (default: all significant sites)",
    )
    args = parser.parse_args(argv)

    sig_path = ROOT / "results" / "significance_critical_only.csv"
    sig_df = pd.read_csv(sig_path)
    default_sites = sig_df[sig_df["sig_crit"]]["site_id"].tolist()
    site_ids = args.sites if args.sites else default_sites

    print(f"Generating prospectivity maps for {len(site_ids)} sites")
    print(f"Scenarios: {DEFAULT_SCENARIOS}\n")

    records = []
    for site_id in sorted(site_ids):
        print(f"  {site_id} ...", end=" ", flush=True)
        result = process_site(site_id)
        records.append(result)
        if result["status"] == "ok":
            print(
                f"ok  high={result['pct_high']:.1f}%  "
                f"mean_σ={result['mean_uncertainty']:.3f}"
            )
        else:
            print(result["status"])

    # Summary table
    ok = [r for r in records if r["status"] == "ok"]
    print(f"\n{'Site':<28} {'High P=1.0':>11} {'Mod P≥0.5':>10} {'Marginal':>9} {'Mean σ':>8}")
    print("-" * 70)
    for r in ok:
        n = r["n_valid_pixels"]
        print(
            f"  {r['site_id']:<26}"
            f" {r['n_high_confidence']:>7} ({r['pct_high']:>4.1f}%)"
            f" {r['n_moderate_confidence']:>7} ({r['n_moderate_confidence']/n*100:>4.1f}%)"
            f" {r['n_marginal']:>7} ({r['n_marginal']/n*100:>4.1f}%)"
            f" {r['mean_uncertainty']:>7.3f}"
        )

    print(
        "\nInterpretation:"
        "\n  P=1.0 (High)     → strong anomaly in ALL threshold scenarios; high-confidence target"
        "\n  P=0.67 (Moderate) → strong in 2/3 scenarios; worth investigating"
        "\n  P=0.33 (Marginal) → strong in 1/3 scenarios only; threshold-sensitive, low confidence"
        "\n  P=0.0  (None)    → never strong; background or non-anomalous terrain"
        "\n  Mean σ           → average pixel uncertainty; sites with higher σ have more"
        "\n                     pixels near classification boundaries"
        "\n\nOutputs: data/sites/{id}/rasters/prospectivity_mean.tif + prospectivity_std.tif"
    )


if __name__ == "__main__":
    main()
