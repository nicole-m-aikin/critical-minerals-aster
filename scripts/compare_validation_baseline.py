#!/usr/bin/env python
"""
Before/after comparison for the MRDS validation step.

Compares the headline significance sets computed from raw MRDS (snapshot in
results/_baseline_pre_validation/) with the current results/ tables computed
from ``validate_mrds()`` output, so any change in a paper number is traceable
to a named validation rule.

Outputs
-------
    results/validation_baseline_diff.csv   one row per site
    stdout                                  headline counts + membership changes

Usage
-----
    conda run -n aster-minerals python scripts/compare_validation_baseline.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
BASE = REPO_ROOT / "results" / "_baseline_pre_validation"
CUR = REPO_ROOT / "results"


def headline(root: Path) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    fdr = pd.read_csv(root / "phase4_fdr_corrected_significance.csv").set_index("site_id")
    p5 = pd.read_csv(root / "phase5_clustering_sensitivity.csv")
    piv = p5.pivot(index="site_id", columns="radius_m", values="sig_05")
    robust_all = set(piv.index[(piv == True).all(axis=1)])  # noqa: E712
    sets = {
        "fdr_binomial": set(fdr.index[fdr["sig_fdr_binomial"]]),
        "fdr_mannwhitney": set(fdr.index[fdr["sig_fdr_mannwhitney"]]),
    }
    sets["fdr_union"] = sets["fdr_binomial"] | sets["fdr_mannwhitney"]
    sets["robust"] = sets["fdr_binomial"] & robust_all
    return fdr, sets


def main() -> None:
    b_fdr, b_sets = headline(BASE)
    c_fdr, c_sets = headline(CUR)

    print(f"{'set':18s} {'baseline':>9s} {'validated':>10s}")
    for key in b_sets:
        print(f"{key:18s} {len(b_sets[key]):9d} {len(c_sets[key]):10d}")
        gained = sorted(c_sets[key] - b_sets[key])
        lost = sorted(b_sets[key] - c_sets[key])
        if gained or lost:
            print(f"{'':18s} gained={gained} lost={lost}")

    cols = ["n_crit", "hits_crit", "enrichment", "p_binomial", "q_binomial", "sig_fdr_binomial"]
    diff = b_fdr[cols].join(c_fdr[cols], lsuffix="_raw", rsuffix="_validated", how="outer")
    diff["n_removed"] = diff["n_crit_raw"] - diff["n_crit_validated"]
    for key in c_sets:
        diff[f"{key}_raw"] = diff.index.isin(b_sets[key])
        diff[f"{key}_validated"] = diff.index.isin(c_sets[key])
    out = CUR / "validation_baseline_diff.csv"
    diff.sort_values("n_removed", ascending=False).to_csv(out)
    print(f"\nTotal critical-mineral records: {int(diff['n_crit_raw'].sum())} -> "
          f"{int(diff['n_crit_validated'].sum())}")
    print(f"Saved {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
