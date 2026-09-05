# Results

Systematic ASTER thermal-infrared (TIR) alteration screening across **109 U.S. critical-mineral
districts**, validated against USGS MRDS deposit records under a corrected statistical framework.
All numbers below are regenerable from the scripts in [`scripts/`](../scripts/) against the
per-site outputs in `results/`; figures are in [`figures/`](../figures/).

---

## Summary

| Stage | Sites (of 109) |
|---|---|
| Site-specific geometric null, raw *p* < 0.05 | 27 |
| + Benjamini–Hochberg FDR, primary binomial test | **20** |
| + BH-FDR, threshold-free continuous-score (Mann–Whitney U) test | 31 (union 33) |
| + robust at every spatial-clustering declustering radius (0–1000 m) | **13** |

**The 13 fully-robust sites are dominated by one deposit system in one climate regime:** nine are
arid, well-exposed **skarn / carbonate-replacement** districts spanning Arizona, New Mexico, Texas
and California. Skarn/Carbonate Replacement is the **only** deposit-system category whose 95%
bootstrap confidence interval on median enrichment excludes 1.0.

The method behaves as a **deposit-type selector, not a universal detector**.

---

## Method (brief)

1. **Band ratios** target silica (B13/B14), carbonate (B13/B12) and mafic (B12/B13) mineralogy.
   Each is split into background / moderate / strong classes by per-scene 70th/90th percentile
   thresholds; an additive 0–6 combined score ≥ 3 defines a "strong anomaly zone."
2. **Site-specific null:** for each site, chance = (strong-anomaly-zone area) ÷ (valid TIR
   footprint area) — not a pooled survey-wide rate. MRDS deposits are clipped to the true diagonal
   ASTER footprint before counting.
3. **Two test families, FDR-corrected independently:** a one-sided exact binomial test and a
   Mann–Whitney U test on the continuous score (reported with its AUC effect size).
4. **Robustness:** footprint-aware Monte Carlo permutation (10,000 iters, seed 42); DBSCAN
   spatial-clustering declustering at 0/250/500/1000 m in each site's local UTM projection;
   classification-threshold sweep (65/85–75/95 percentiles, score cutoffs 2–4).
5. **Deposit-system classification** by each site's identity as a named USGS district
   (`results/deposit_system_categories.csv`), not MRDS's unreliable automated field.

Half the sample (58 sites) was added in a pre-registered expansion: each carries a written
NULL / POSITIVE / POSITIVE-ISH / NULL-TO-WEAK detectability expectation in
`results/deposit_system_categories.csv` (`[EXPANSION +58]` rows), recorded before its imagery was
processed.

---

## Site-level results

### Robust to every check (13 sites)

FDR-significant on the primary binomial test **and** significant at every declustering radius.

| Site | State | Deposit system | Enrichment | Note |
|---|---|---|---|---|
| Bisbee | AZ | Skarn / carbonate-replacement | 3.4× | original |
| Courtland-Gleeson | AZ | Skarn / carbonate-replacement | 2.7× | original |
| Tombstone | AZ | Carbonate-replacement | 2.2× | original |
| Christmas | AZ | Cu skarn | 1.9× | original (threshold-sensitive at tightest cutoff) |
| Magdalena–Kelly | NM | Zn-Pb skarn / replacement | 2.0× | replication, independent province |
| Organ Mountains | NM | Cu-Pb-Zn-Ag skarn | 2.3× | replication |
| Lake Valley | NM | Ag carbonate-replacement | 2.6× | replication |
| Shafter | TX | Ag carbonate-replacement | 2.3× | replication, Trans-Pecos |
| Eagle Mountain | CA | Fe skarn | 2.4× | replication, cleanest independent Fe-skarn test |
| Providence Mountains | CA | Ag-Pb replacement | 1.7× | replication (binomial test) |
| crooks_gap | WY | Roll-front uranium | 3.3× | **unexpected** — see Limitations |
| karnes_county | TX | Roll-front uranium | 3.1× | **unexpected** — see Limitations |
| holden | WA | Volcanogenic massive sulfide | 2.4× | **unexpected**, likely mine-tailings signal — see Limitations |

Two sites significant in an earlier 51-site version (Eureka, NV; Lordsburg, NM) are **not** in
this set: both lose significance once spatially-clustered MRDS records are deduplicated at the
1000 m radius (Eureka *p* = 0.76; Lordsburg *p* = 0.088).

### Filtered out by the robustness cascade

`tri_state`, `uravan`, `lisbon_valley`, `mascot_jefferson_city` are FDR-significant on **both**
primary tests but collapse under declustering (e.g. Tri-State: *p* ≈ 0 at 0–250 m → *p* ≈ 1 by
500 m) — MRDS point clusters, not alteration.

---

## Deposit-system results

Category median enrichment (observed ÷ site-specific null hit rate), 95% bootstrap CI, and FDR
survival rate (`results/phase6_deposit_system_categories.csv`):

| Category | n | Median enrichment | 95% CI | FDR-significant |
|---|---|---|---|---|
| **Skarn / Carbonate Replacement** | 28 | **1.64×** | **1.02 – 2.25** | 13 / 28 |
| Volcanogenic / VMS | 9 | 1.19× | 0.27 – 1.41 | 4 / 9 |
| Mafic / Ultramafic | 7 | 1.16× | 0.15 – 1.35 | 0 / 7 |
| Uranium / Energy | 9 | 1.10× | 0.19 – 3.13 | 5 / 9 |
| Sediment-hosted | 12 | 0.86× | 0.46 – 1.73 | 5 / 12 |
| Porphyry Cu-Mo-Au | 17 | 0.78× | 0.32 – 1.36 | 4 / 17 |
| Alkaline / Carbonatite | 10 | 0.70× | 0.00 – 1.34 | 0 / 10 |
| Epithermal | 11 | 0.44× | 0.20 – 1.05 | 0 / 11 |

**Skarn / carbonate-replacement is the only category whose CI excludes 1.0.** Its robust set spans
four states and three physiographic provinces — the signal is a regime defined by deposit genesis
and surface exposure, not a single mining belt.

**Pre-registered expectations that held:** mafic/ultramafic (0 / 6 new sites significant),
alkaline/carbonatite (0 / 4), Carlin-type sediment-hosted (0 / 4), and a systematic humid-terrain
skarn control cohort — six districts of the *same* genesis as the robust arid sites but in
vegetated Appalachian/Cascade terrain — **0 / 6 FDR-significant** (three raw-enriched, none
surviving correction). Genesis held fixed, terrain varied, signal disappeared.

**Pre-registered expectations that failed:** four epithermal districts with well-exposed,
unburied alteration caps were added expecting a modest signal; **0 / 4 are FDR-significant** and
the category median stays 0.44× — the lowest of any category. Exposed alteration caps do not make
epithermal systems TIR-detectable here. Four additional porphyry Cu-Mo-Au districts: 0 / 4.

The **Uranium / Energy** category is bimodal: two roll-front districts (crooks_gap,
karnes_county) are strongly and robustly enriched while the other seven — including the reference
site Gas Hills — are indistinguishable from chance. The category median (1.10×) and its very wide
CI (0.19–3.13) reflect that split, not a uniform effect.

---

## Discovery-bias stratification

To test whether the spatial association could reflect historical ascertainment bias (deposits
found by the same visible alteration ASTER detects), critical-mineral deposits with a usable
discovery / first-production year were split at 1950 and each cohort tested against its site's own
null, unconditionally across all 99 sites with dated records.

Result: **no enrichment in either cohort** (pre-1950 1.01×, *p* = 0.47; post-1950 0.92×,
*p* = 0.73). This analysis is uninformative, not exculpatory — it does not resolve the
ascertainment-bias question in either direction, given ~14% dated-record coverage and a
demographically skewed dated subsample. It is reported as strictly secondary.

---

## Limitations

1. **The method is not demonstrated to be alteration-specific.** Three sites — `crooks_gap` and
   `karnes_county` (roll-front sandstone uranium) and `holden` (VMS) — are FDR-significant on both
   tests and robust to every filter, despite deposit genesis that predicts no intrusion-related
   surface alteration halo. `holden` is a former mine with a large exposed sulphide-oxidation
   tailings field and its signal is plausibly mine waste rather than a natural halo; the two
   uranium districts are well-exposed and may carry a genuine limonite/hematite/calcite alteration
   front at the redox boundary. These are only 2 / 9 uranium and 1 / 9 VMS sites — the effect is
   not category-wide — but until they are explained the method can only be said to predict deposit
   *location* in some settings, not to detect *alteration*.

2. **SWIR bands are unavailable** for these areas in LP DAAC v004, so argillic / phyllic alteration
   (relevant to Carlin-type gold and some porphyry systems) is invisible to this method regardless
   of any statistical correction. The below-chance sediment-hosted and epithermal medians may
   partly reflect this detection gap.

3. **Arid / semi-arid applicability only.** The humid-terrain control cohort shows the signal does
   not survive vegetation cover even when deposit genesis is favourable.

4. **Small per-category sample sizes** (7–28 sites) mean the deposit-system comparison supports a
   directional ranking, not precise category-level estimates; every CI except skarn/CRD overlaps
   1.0.

5. **Category assignment is not fully blind to outcome.** The 58 expansion-site categories and
   expectations were fixed before results; the original 51 were not.

6. **BH-FDR's dependence assumption** was not stress-tested against the more conservative
   Benjamini–Yekutieli variant.
