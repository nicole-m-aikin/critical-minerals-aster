#!/usr/bin/env python3
"""Combined TIR host screening vs GEMx SWIR clay check (primary FDR sites).

Reads results/gemx_host_vs_alteration.csv and phase4 stratified enrichment.
Writes figures/verification_tir_swir_combined.png.

Usage:
    conda run -n aster-minerals python scripts/fig_tir_swir_combined.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from critical_minerals_aster.significance import significant_site_sets

REPO = Path(__file__).resolve().parents[1]
CSV = REPO / "results" / "gemx_host_vs_alteration.csv"
PHASE4 = REPO / "results" / "phase4_fdr_corrected_significance.csv"
OUT = REPO / "figures" / "verification_tir_swir_combined.png"

VERDICT_STYLE = {
    "complementary": ("#2c7fb8", "o"),
    "overlapping": ("#c0392b", "s"),
    "insufficient": ("#e67e22", "D"),
    "no_coverage": ("0.65", "x"),
    "error": ("0.4", "+"),
}


def main() -> None:
    gemx = pd.read_csv(CSV)
    p4 = pd.read_csv(PHASE4)[["site_id", "enrichment", "sig_fdr_binomial"]]
    fdr, robust = significant_site_sets(REPO)
    df = gemx.merge(p4, on="site_id", how="left")
    df["acquisition_replicate"] = df["site_id"].isin(
        {"courtland_gleeson", "pioche", "providence_mountains"}
    )
    df["robust"] = df["site_id"].isin(robust)

    ok = df[df["status"] == "ok"].copy()
    if len(ok) < 3:
        raise SystemExit(f"Need ≥3 ok primary sites; got {len(ok)}")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.5, 5.8),
                                     gridspec_kw={"width_ratios": [1.35, 1]})

    # A — enrichment vs Al-OH Jaccard
    for verdict, (colour, marker) in VERDICT_STYLE.items():
        sub = ok[ok["verdict"] == verdict]
        if sub.empty:
            continue
        n_all = int((df["verdict"] == verdict).sum())
        ax_a.scatter(
            sub["enrichment"], sub["jaccard_aloh"],
            c=colour, marker=marker, s=55, edgecolors="k", linewidths=0.4,
            label=f"{verdict} (n={n_all})", zorder=3,
        )
    for _, r in ok.iterrows():
        weight = "bold" if r["acquisition_replicate"] else "normal"
        ax_a.annotate(
            r["site_id"], (r["enrichment"], r["jaccard_aloh"]),
            xytext=(4, 3), textcoords="offset points", fontsize=7,
            fontweight=weight, color="0.2",
        )
    ax_a.axhline(0.15, color="0.5", ls="--", lw=0.9, label="Jaccard = 0.15")
    ax_a.set_xlabel("TIR enrichment (bedrock/cover stratified null)")
    ax_a.set_ylabel("GEMx Al-OH Jaccard with TIR strong zones")
    ax_a.set_title("A. TIR host signal vs SWIR clay overlap", loc="left", fontsize=11)
    ax_a.legend(fontsize=8, loc="upper right", framealpha=0.92)
    ax_a.grid(alpha=0.25)
    ax_a.set_xlim(0, max(4.5, float(ok["enrichment"].max()) * 1.05))
    ax_a.set_ylim(0, max(0.22, float(ok["jaccard_aloh"].max()) * 1.15))

    # B — verdict counts among the 19 primary FDR sites
    order = ["complementary", "overlapping", "insufficient", "no_coverage", "error"]
    counts = df["verdict"].value_counts()
    labels = [v for v in order if v in counts.index]
    vals = [int(counts[v]) for v in labels]
    colours = [VERDICT_STYLE.get(v, ("0.5", "o"))[0] for v in labels]
    y = np.arange(len(labels))
    ax_b.barh(y, vals, color=colours, edgecolor="k", linewidth=0.4)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(labels)
    ax_b.invert_yaxis()
    for yi, v in zip(y, vals):
        ax_b.text(v + 0.1, yi, str(v), va="center", fontsize=10)
    ax_b.set_xlabel("Primary FDR sites (n = 19)")
    ax_b.set_xlim(0, max(vals) + 2.5)
    ax_b.set_title("B. Host-vs-alteration verdicts", loc="left", fontsize=11)
    n_comp = int(counts.get("complementary", 0))
    ax_b.text(
        0.98, 0.02,
        f"{n_comp}/{len(df)} complementary\n"
        f"(TIR ≠ SWIR clay)\n"
        f"{len(fdr)} FDR · {len(robust)} clustering-robust",
        transform=ax_b.transAxes, ha="right", va="bottom", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.7"),
    )

    fig.suptitle(
        "Combined method: ASTER TIR host screening + GEMx SWIR clay check",
        fontsize=12, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT}")
    print(counts.to_string())


if __name__ == "__main__":
    main()
