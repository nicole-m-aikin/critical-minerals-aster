#!/usr/bin/env python
"""
Figure 3 -- representative success/failure comparison.

Scientific question: what physically distinguishes a site where the method
works from one where it doesn't, holding other variables as close to fixed
as possible?

Five sites (post-109-site expansion), one replication and three contrasts:
  - Bisbee (success, AZ) vs Magdalena-Kelly (success, NM) -- a REPLICATION
    across an independent region. Magdalena-Kelly is a Zn-Pb skarn/carbonate-
    replacement district in a different state and structural province that,
    like Bisbee, survives every robustness check including spatial-clustering
    correction. The arid skarn/CRD signal is no longer confined to one
    Cochise County belt.
  - Bisbee (success) vs Cornwall, PA (failure) -- the SAME deposit genesis
    (Cornwall-type Fe skarn: calc-silicate replacement of carbonate against
    an intrusion) but humid, forested Piedmont instead of arid Sonoran
    Desert. Cornwall is raw-enriched (1.7x) but does not survive FDR. This
    is the systematic version of the old single Bisbee-vs-Ducktown pair:
    genesis held fixed, terrain varied -> VEGETATION/EXPOSURE isolated.
  - Bisbee vs Ducktown, TN (failure) -- humid-terrain VMS; climate + genesis
    both unfavourable.
  - Bisbee vs Climax, CO (failure) -- arid, but a deep porphyry-Mo stockwork
    with no surface alteration expression. DEPTH/EXPOSURE isolated.

Form: grouped bar chart, hit rate vs. site-specific null side by side per
site, plus enrichment as a direct-labeled value -- avoids a dual-axis chart
by keeping both bars on the same 0-100% scale and reporting enrichment as
text rather than a second y-axis.

Usage:
    conda run -n aster-minerals python scripts/fig03_success_failure.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent

SITES = ["bisbee", "magdalena_kelly", "cornwall_pa", "ducktown", "climax"]
LABELS = {
    "bisbee": "Bisbee, AZ\nSkarn/carbonate Cu-Ag\nArid, exposed\nSUCCESS",
    "magdalena_kelly": "Magdalena–Kelly, NM\nZn-Pb skarn/replacement\nArid, exposed\nSUCCESS (independent region)",
    "cornwall_pa": "Cornwall, PA\nCornwall-type Fe skarn\nHumid, forested\nraw 1.7×, NOT FDR-sig.",
    "ducktown": "Ducktown, TN\nVMS Cu-Zn\nHumid, forested\nFAILURE (climate)",
    "climax": "Climax, CO\nPorphyry Mo\nArid, but deep/buried\nFAILURE (depth)",
}
OUTCOME_COLOR = {
    "bisbee": "#2E7D32", "magdalena_kelly": "#2E7D32",
    "cornwall_pa": "#B71C1C", "ducktown": "#B71C1C", "climax": "#B71C1C",
}


def main() -> None:
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv").set_index("site_id")
    fdr = pd.read_csv(REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv").set_index("site_id")

    fig, ax = plt.subplots(figsize=(11, 5.8))
    x = np.arange(len(SITES))
    width = 0.32

    observed = [ss.loc[s, "hit_rate_crit_pct"] for s in SITES]
    null = [ss.loc[s, "p0_null"] * 100 for s in SITES]
    enrichment = [ss.loc[s, "enrichment"] for s in SITES]
    q = [fdr.loc[s, "q_binomial"] for s in SITES]

    b1 = ax.bar(x - width / 2, observed, width, label="Observed critical-mineral hit rate",
                color="#D55E00", zorder=3)
    b2 = ax.bar(x + width / 2, null, width, label="Site-specific null (chance)",
                color="#999999", zorder=3)

    for i, (o, n, e, qv) in enumerate(zip(observed, null, enrichment, q)):
        ax.text(i, max(o, n) + 2.5, f"{e:.2f}× enrichment\nq={qv:.3f}",
                ha="center", fontsize=8.5, color="#222222")

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in SITES], fontsize=9, linespacing=1.7)
    ax.set_ylabel("Hit rate (%)")
    ax.set_ylim(0, max(observed) + 12)
    ax.set_title(
        "Figure 3 — matched success/failure comparison\n"
        "Arid skarn signal replicates in an independent region (Magdalena–Kelly); the same genesis "
        "fails in humid terrain (Cornwall), as does VMS in humid terrain (Ducktown) and buried porphyry (Climax)",
        fontsize=9.5,
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.7, zorder=0)

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.24)
    out_path = REPO_ROOT / "figures" / "fig03_success_failure.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
