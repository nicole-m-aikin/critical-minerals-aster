#!/usr/bin/env python
"""
Phase 6 -- deposit-system category analysis.

Groups the 45 sites into 9 coarse, geologically-motivated deposit-system
categories (results/deposit_system_categories.csv -- hand-curated from
standard economic-geology identity of each named USGS district, cross-
referenced against docs/results.md's existing literature labels where
available; see that CSV's "confidence"/"basis" columns for provenance of
each assignment). Categories intentionally do NOT force every site into a
clean bucket -- "Mixed/Other" holds sites that are genuinely polygenetic
(Battle Mountain), non-metallic (Green River), or sui generis lithium-clay/
caldera systems with no dedicated slot in the 24-system OFR taxonomy
(McDermitt, Thacker Pass) -- per the explicit instruction not to force
ambiguous sites into categories.

For each category, computes:
  - n_sites
  - median enrichment (site_specific_null_significance.csv's `enrichment`
    column: observed hit rate / site-specific null hit rate)
  - a bootstrap 95% CI on that median
  - fraction of sites surviving BH-FDR on EITHER the primary binomial test
    or the secondary continuous-score test (phase4_fdr_corrected_significance.csv)

Bootstrap methodology
----------------------
scipy.stats.bootstrap(data, np.median, n_resamples=10_000, method="percentile",
confidence_level=0.95, random_state=42).

Method note: an initial run used method="basic" (the "reverse percentile"
correction) and produced an impossible negative lower bound (-1.03) for
enrichment -- a ratio that cannot be negative -- for the Skarn/Carbonate
Replacement category, whose 7 values are heavily right-skewed (Bisbee and
Eureka near 3.5-3.7x, the rest near or below 1x). "basic" computes
2*point_estimate - percentile, which can legitimately fall outside the
observed data range for small, skewed samples. Switched to method=
"percentile" (using the resample distribution's own 2.5th/97.5th
percentiles directly as the bounds), which is guaranteed to stay within the
range of values the bootstrap actually produced -- always non-negative
here, at some cost in accuracy relative to a bias-corrected method for
strongly skewed distributions, but far safer to report and interpret for a
small, non-negative ratio in a short paper's supplementary table.

What resampling with replacement means: from a category's n enrichment
values, draw a new sample of size n BY SAMPLING WITH REPLACEMENT (so the
same site's value can appear zero, one, or several times in a given
resample) -- do this 10,000 times, compute the median of each resample, and
use the resulting distribution of 10,000 medians as an empirical estimate
of how much the true median could plausibly vary given only n observations.

What the CI means: a 95% bootstrap CI is the range that would contain the
true population median in ~95% of hypothetical repeated samples of this
size from the same underlying population -- NOT "95% probability the true
median is in this range" (that is a Bayesian-sounding misreading of a
frequentist interval).

Why bootstrap is appropriate here: category enrichment values are a small,
non-normal, irregularly-spaced sample (n_sites per category ranges 1-12)
-- there is no clean parametric formula for the sampling distribution of a
MEDIAN of arbitrary shape at this n, and the bootstrap needs no
distributional assumption beyond "this sample is representative of the
population it was drawn from."

Limitations with small n: with n=1 (Uranium/Energy: gas_hills only) or n=2
(Volcanogenic/VMS: ducktown, jerome), every bootstrap resample is
necessarily built from the same 1-2 values, so the "CI" collapses to a
degenerate point or a handful of repeated values -- it conveys no real
uncertainty information and is reported as "insufficient n" rather than a
misleadingly narrow interval. Categories need at least ~5 sites before a
bootstrap CI is worth interpreting as anything beyond a rough sketch, and
even at n=6-7 the interval should be read as suggestive, not precise --
this project reports it anyway, labelled, because CLAUDE.md's operating
rule is to always report uncertainty rather than hide it, not because n=6
enrichment values support a strong quantitative claim.

Outputs
-------
    results/phase6_deposit_system_categories.csv

Usage:
    conda run -n aster-minerals python scripts/phase6_deposit_system_categories.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import bootstrap

REPO_ROOT = Path(__file__).parent.parent
SEED = 42
N_RESAMPLES = 10_000
MIN_N_FOR_BOOTSTRAP = 3


def main() -> None:
    cats = pd.read_csv(REPO_ROOT / "results" / "deposit_system_categories.csv")
    ss = pd.read_csv(REPO_ROOT / "results" / "site_specific_null_significance.csv")
    fdr = pd.read_csv(REPO_ROOT / "results" / "phase4_fdr_corrected_significance.csv")

    df = cats.merge(ss[["site_id", "n_crit", "hit_rate_crit_pct", "p0_null", "enrichment"]],
                     on="site_id", how="left")
    df = df.merge(fdr[["site_id", "sig_fdr_binomial", "sig_fdr_mannwhitney"]],
                   on="site_id", how="left")
    df["sig_fdr_either"] = df["sig_fdr_binomial"].fillna(False) | df["sig_fdr_mannwhitney"].fillna(False)

    rng = np.random.default_rng(SEED)
    rows = []

    for category, grp in df.groupby("category"):
        enrichments = grp["enrichment"].dropna().values
        n_sites = len(grp)
        n_valid = len(enrichments)

        median_enr = float(np.median(enrichments)) if n_valid else float("nan")

        ci_lo, ci_hi, ci_note = float("nan"), float("nan"), "insufficient n"
        if n_valid >= MIN_N_FOR_BOOTSTRAP:
            res = bootstrap(
                (enrichments,), np.median,
                n_resamples=N_RESAMPLES, method="percentile",
                confidence_level=0.95,
                random_state=np.random.default_rng(SEED),
            )
            ci_lo = float(res.confidence_interval.low)
            ci_hi = float(res.confidence_interval.high)
            ci_note = "ok"

        n_fdr = int(grp["sig_fdr_either"].sum())
        n_high_conf = int((grp["confidence"] == "high").sum())

        rows.append({
            "category": category,
            "n_sites": n_sites,
            "n_high_confidence": n_high_conf,
            "median_enrichment": round(median_enr, 2) if n_valid else None,
            "bootstrap_ci_low": round(ci_lo, 2) if ci_note == "ok" else None,
            "bootstrap_ci_high": round(ci_hi, 2) if ci_note == "ok" else None,
            "ci_note": ci_note,
            "n_fdr_significant": n_fdr,
            "frac_fdr_significant": round(n_fdr / n_sites, 2),
            "sites": ", ".join(sorted(grp["site_id"])),
        })

    out = pd.DataFrame(rows).sort_values("median_enrichment", ascending=False)
    out_path = REPO_ROOT / "results" / "phase6_deposit_system_categories.csv"
    out.to_csv(out_path, index=False)

    print(f"Saved {out_path}\n")
    disp_cols = ["category", "n_sites", "n_high_confidence", "median_enrichment",
                 "bootstrap_ci_low", "bootstrap_ci_high", "ci_note",
                 "n_fdr_significant", "frac_fdr_significant"]
    print(out[disp_cols].to_string(index=False))

    print(f"\nTotal sites classified: {df['site_id'].nunique()} / {ss['site_id'].nunique()}")
    low_conf = df[df["confidence"] == "low"]
    print(f"Low-confidence / genuinely mixed assignments ({len(low_conf)}): "
          f"{', '.join(low_conf['site_id'])}")


if __name__ == "__main__":
    main()
