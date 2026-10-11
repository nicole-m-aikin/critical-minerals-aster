#!/usr/bin/env python
"""Figure: local ASTER TIR vs cloud-native implementation choices.

Reads results/gee_parity_agreement.csv from scripts/gee_parity_benchmark.py.

Four panels:

    A  Goldfield strong zones, local (clip-then-percentile)
    B  Goldfield after the clip-order probe (percentiles on the full scene)
    C  Strong-zone Jaccard by site × probe  — do the *same pixels* fire?
    D  Mean |Δ strong-zone fraction| by probe — does *p_null* move?

Usage:
    conda run -n aster-minerals python scripts/fig_gee_parity.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.gee_parity import (
    DEFAULT_PARITY_SITES,
    probe_clip_order,
    run_local_workflow,
)
from critical_minerals_aster.paths import site_paths_for

OUT = ROOT / "figures" / "gee_parity_agreement.png"
CSV = ROOT / "results" / "gee_parity_agreement.csv"

PROBE_LABELS = {
    "probe:clip_order": "clip order",
    "probe:histogram_percentile": "histogram pctl",
    "probe:radiance_ucc": "radiance UCC",
    "probe:reproject_4326": "EPSG:4326",
    "gee_live": "GEE live",
}
PROBE_COLORS = {
    "probe:clip_order": "#D55E00",
    "probe:histogram_percentile": "#0072B2",
    "probe:radiance_ucc": "#009E73",
    "probe:reproject_4326": "#CC79A7",
    "gee_live": "#000000",
}

# nodata / not-strong / strong — binary on purpose so 10% orange reads
STRONG_CMAP = ListedColormap(["#FFFFFF", "#C5C0B8", "#D55E00"])

SITE_SHORT = {
    "goldfield": "goldfield",
    "bisbee": "bisbee",
    "christmas": "christmas",
    "tombstone": "tombstone",
    "mcdermitt": "mcdermitt",
    "eureka": "eureka",
    "mountain_pass": "mt. pass",
    "marysvale": "marysvale",
}


def _load_table() -> pd.DataFrame:
    if not CSV.is_file():
        raise FileNotFoundError(
            f"{CSV} is missing. Run scripts/gee_parity_benchmark.py first."
        )
    return pd.read_csv(CSV)


def _valid_window(valid: np.ndarray) -> tuple[slice, slice]:
    rows = np.where(np.any(valid, axis=1))[0]
    cols = np.where(np.any(valid, axis=0))[0]
    return slice(int(rows[0]), int(rows[-1]) + 1), slice(int(cols[0]), int(cols[-1]) + 1)


def _goldfield_panels(ax_local, ax_probe) -> None:
    site = load_site_by_id("goldfield", ROOT / "sites")
    paths = site_paths_for(site, ROOT)
    local = run_local_workflow(site, paths)
    probed = probe_clip_order(local, site)
    if probed is None:
        ax_local.set_axis_off()
        ax_probe.set_axis_off()
        ax_local.set_title("Goldfield local (rasters missing)")
        return

    rsl, csl = _valid_window(local.valid)

    def _show(ax, score, title, frac):
        window_valid = local.valid[rsl, csl]
        window_score = score[rsl, csl]
        show = np.zeros(window_score.shape, dtype=np.uint8)
        show[window_valid & (window_score < 3)] = 1
        show[window_valid & (window_score >= 3)] = 2
        ax.imshow(show, cmap=STRONG_CMAP, vmin=0, vmax=2, interpolation="nearest")
        ax.set_title(f"{title}\nstrong fraction = {frac:.3f}", fontsize=10, pad=6)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)

    local_frac = float(np.mean(local.combined[local.valid] >= 3))
    probe_frac = float(np.mean(probed.combined[local.valid] >= 3))
    _show(ax_local, local.combined, "A   Goldfield · local  (clip, then percentile)", local_frac)
    _show(
        ax_probe,
        probed.combined,
        "B   Goldfield · clip-order probe  (percentile on full scene)",
        probe_frac,
    )
    ax_local.legend(
        handles=[
            Patch(facecolor="#FFFFFF", edgecolor="0.6", label="nodata"),
            Patch(facecolor="#C5C0B8", edgecolor="0.6", label="score < 3"),
            Patch(facecolor="#D55E00", edgecolor="0.6", label="strong (score ≥ 3)"),
        ],
        loc="lower left",
        fontsize=7.5,
        framealpha=0.92,
        borderpad=0.4,
    )


def _jaccard_panel(ax, df: pd.DataFrame) -> None:
    sites = [s for s in DEFAULT_PARITY_SITES if s in set(df["site_id"])]
    probes = [p for p in PROBE_LABELS if p in set(df["probe"])]
    x = np.arange(len(sites))
    width = 0.78 / max(1, len(probes))
    for i, probe in enumerate(probes):
        sub = df[df["probe"] == probe].set_index("site_id")
        vals = [
            float(sub.loc[s, "jaccard_strong"]) if s in sub.index else np.nan
            for s in sites
        ]
        ax.bar(
            x + (i - (len(probes) - 1) / 2) * width,
            vals,
            width=width * 0.90,
            color=PROBE_COLORS[probe],
            label=PROBE_LABELS[probe],
            edgecolor="white",
            linewidth=0.3,
        )
    ax.set_xticks(x)
    ax.set_xticklabels([SITE_SHORT.get(s, s) for s in sites], fontsize=8, rotation=20, ha="right")
    ax.set_ylabel("Strong-zone Jaccard vs local")
    ax.set_ylim(0, 1.08)
    ax.axhline(1.0, color="0.55", lw=0.7, ls=":")
    ax.set_title("C   Do the same pixels fire?", fontsize=10)
    ax.legend(fontsize=7.5, loc="lower left", ncol=2, framealpha=0.92)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="0.88", lw=0.5)


def _attribution_panel(ax, df: pd.DataFrame) -> None:
    probes = [p for p in PROBE_LABELS if p in set(df["probe"])]
    means = [
        float(np.nanmean(np.abs(df.loc[df["probe"] == p, "strong_frac_delta"])))
        for p in probes
    ]
    # Largest effect at the top.
    order = np.argsort(means)
    probes = [probes[i] for i in order]
    means = [means[i] for i in order]
    y = np.arange(len(probes))
    ax.barh(y, means, color=[PROBE_COLORS[p] for p in probes], edgecolor="white", height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels([PROBE_LABELS[p] for p in probes], fontsize=9)
    ax.set_xlabel("|Δ strong-zone fraction|   (mean of 8 sites)")
    ax.set_title("D   Does the area behind p_null move?", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    xmax = max(means) if means else 0.01
    ax.set_xlim(0, xmax * 1.22 if xmax > 0 else 0.03)
    for yi, val in zip(y, means):
        ax.text(val + xmax * 0.03, yi, f"{val:.3f}", va="center", fontsize=8)


def main() -> int:
    df = _load_table()
    fig = plt.figure(figsize=(12.2, 8.6))
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.12, 1.0], hspace=0.42, wspace=0.30,
        left=0.07, right=0.98, top=0.90, bottom=0.10,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    _goldfield_panels(ax_a, ax_b)
    _jaccard_panel(ax_c, df)
    _attribution_panel(ax_d, df)

    fig.suptitle(
        "ASTER TIR implementation parity   ·   local rasterio vs cloud-native choices",
        fontsize=13,
        y=0.97,
    )
    fig.text(
        0.50,
        0.025,
        "8-site panel.  GEE collection ASTER/AST_L1T_003; local granules are LP DAAC v004.  "
        "Clip-then-percentile is the local contract.  "
        "Live Earth Engine:  python -m critical_minerals_aster gee-parity --gee",
        ha="center",
        fontsize=8,
        color="0.35",
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=170, facecolor="white")
    plt.close(fig)
    print(f"Wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
