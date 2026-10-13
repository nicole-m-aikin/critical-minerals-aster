#!/usr/bin/env python3
"""Before/after destriping diagnostic for one site's AST_05 QI mosaic.

Rebuilds the emissivity index mosaic twice (destripe off / on) from cached
granules, measures detector-mean amplitude on B11 and QI, and writes a
two-panel figure.

Usage:
    conda run -n aster-minerals python scripts/destripe_diagnostic.py --site pioche
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.spectral import (
    alteration_indices,
    destripe_amplitude,
    destripe_emissivity_bands,
    load_ast05_granule,
)

REPO = Path(__file__).resolve().parents[1]


def _qi_from_granule(src_dir: Path, gid: str, *, destripe: bool) -> tuple[np.ndarray, dict]:
    bands, _, _, _, _ = load_ast05_granule(src_dir, gid)
    stats = {"destripe": False}
    if destripe:
        bands, stats = destripe_emissivity_bands(bands)
    qi, _, _ = alteration_indices(*(bands[b] for b in (10, 11, 12, 13, 14)))
    return qi, stats


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--site", default="pioche")
    args = ap.parse_args(argv)

    site = load_site_by_id(args.site, REPO / "sites")
    paths = site_paths_for(site, REPO)
    sel_path = paths.aster_dir / f"{site.id}_mosaic_granules.json"
    if not sel_path.is_file():
        raise SystemExit(f"Missing {sel_path}; run an emissivity rebuild first")
    sel = json.loads(sel_path.read_text())
    gids = [r["granule_id"] for r in sel["granules"]]
    if not gids:
        raise SystemExit(f"No granules in {sel_path}")

    ref_gid = gids[0]
    qi_off, _ = _qi_from_granule(paths.aster_dir, ref_gid, destripe=False)
    qi_on, dstats = _qi_from_granule(paths.aster_dir, ref_gid, destripe=True)
    amp_off = destripe_amplitude(qi_off)
    amp_on = destripe_amplitude(qi_on)

    # Plot native-geometry QI (where detector rows still align). Mosaic amplitude
    # after warp is not a valid striping metric.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    finite_both = qi_off[np.isfinite(qi_off) & np.isfinite(qi_on)]
    vmin, vmax = np.percentile(finite_both, [2, 98]) if finite_both.size else (0.9, 1.1)
    for ax, arr, title, amp in (
        (axes[0], qi_off, f"{site.id} QI — no destripe", amp_off),
        (axes[1], qi_on, f"{site.id} QI — moment-match destripe", amp_on),
    ):
        im = ax.imshow(arr, cmap="viridis", vmin=vmin, vmax=vmax, origin="upper")
        ax.set_title(f"{title}\ndetector-mean amplitude = {amp:.5f}", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(
        f"ASTER TIR QI destriping — native granule {ref_gid}\n"
        f"QI amplitude: {amp_off:.5f} → {amp_on:.5f}; B11 axis={dstats.get('axis')}",
        fontsize=11,
    )
    fig.tight_layout()
    out = REPO / "figures" / f"verification_destripe_{site.id}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)

    report = {
        "site_id": site.id,
        "reference_granule": ref_gid,
        "qi_amplitude_granule_off": round(amp_off, 6),
        "qi_amplitude_granule_on": round(amp_on, 6),
        "reduction_factor": round(amp_off / amp_on, 1) if amp_on > 0 else None,
        "destripe_stats": dstats,
        "figure": str(out.relative_to(REPO)),
    }
    report_path = REPO / "results" / f"destripe_diagnostic_{site.id}.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"Saved {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
