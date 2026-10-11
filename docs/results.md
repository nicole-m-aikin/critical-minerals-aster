# Results

Systematic ASTER thermal-infrared (TIR) alteration screening across **109 U.S. critical-mineral
districts**, validated against USGS MRDS deposit records under a corrected statistical framework.
All numbers below are regenerable from the scripts in [`scripts/`](../scripts/) against the
per-site outputs in `results/`; figures are in [`figures/`](../figures/).

---

## Summary

Input since 2026-10-09: **AST_05 v004 daytime surface emissivity**, Ninomiya quartz/carbonate/mafic
indices, per-pixel QA and AG100 water/vegetation masks. Primary null: **bedrock-vs-cover stratified**.
See [`docs/data_verification.md`](data_verification.md) §2.4–2.6 for why: the earlier radiance
results depended on acquisition temperature. Earlier numbers are in `results/_baseline_pre_clearday/`
(mixed day/night radiance) and `results/_baseline_pre_emissivity/` (clear-day radiance).

| Stage | Sites (of 108 testable) |
|---|---|
| Bedrock/cover stratified null, raw *p* < 0.05 | 23 |
| + Benjamini–Hochberg FDR, primary binomial test | **19** |
| + BH-FDR, threshold-free continuous-score (Mann–Whitney U) test | 32 |
| + robust at every spatial-clustering declustering radius (0–1000 m) | **10** |
| Geometric null, BH-FDR binomial (for comparison) | 23 |

**Independent replication of the 19.** Night AST_05 mosaics and the AG100 multi-year mean emissivity
(both outside the score) check whether a site's result survives a different acquisition:

| Checks passed (clustering-robust, night *p* < 0.05, AG100 FDR) | Sites |
|---|---|
| 3 | courtland_gleeson, pioche, providence_mountains |
| 2 | eureka, lordsburg, shafter, sierrita, tombstone |
| 1 | contact_district, iron_king, magdalena_kelly, san_francisco_frisco, silver_peak, tintic, tri_state, victorio |
| 0 | darwin, mcdermitt, mountain_pass |

Under the finer 8-class SGMC lithology null, 13 of the 19 stay significant. Of the three-check
sites, only Providence Mountains survives; Courtland-Gleeson and Pioche collapse to ≤ 1×, i.e.
their zones track the carbonate host package. That makes Providence Mountains the only site to
pass every replication and verification test.

Most primary sites are still arid **skarn / carbonate-replacement** and polymetallic replacement
districts (14 of 28 skarn/CRD sites are FDR-significant). But the category CI now touches 1.0
(1.44×, 1.00–1.63). Only Mafic/Ultramafic (1.59×, 1.26–1.96, 1 of 7 FDR) and Mixed/Other have CIs
above 1.

The method behaves as a **deposit-type selector, not a universal detector**, and the per-site
results are less stable than the earlier headline implied: 11 of the original 20 FDR sites do not
reproduce on emissivity.

---

## Method (brief)

1. **Emissivity indices.** AST_05 v004 TES emissivity from clear, daytime, May–September scenes
   (never mixed with night). Bad, cloudy and cloud-adjacent pixels are masked per pixel, and AG100
   water and vegetation (NDVI > 0.3) are removed. Ninomiya quartz QI = B11²/(B10·B12), carbonate
   CI = B13/B14 and mafic MI = B12/B13 are computed per granule and priority-stacked, not blended.
   Each index is split into background / moderate / strong classes by per-scene 70th/90th
   percentile thresholds; an additive 0–6 combined score ≥ 3 defines a "strong anomaly zone."
2. **Bedrock/cover stratified null:** the TIR footprint is split into bedrock, cover (SGMC
   "Unconsolidated") and unmapped. Each deposit's chance of a hit is the zone coverage of its own
   stratum; the test is an exact Poisson-binomial tail. The geometric null (zone area ÷ footprint
   area) is kept for comparison. MRDS deposits are clipped to the valid footprint before counting.
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

### FDR-significant sites (19), with replication checks

Enrichment = hits ÷ stratified expected hits. Clustering = significant at every declustering radius.
Night = night AST_05 enrichment (✓ *p* < 0.05). AG100 = mean-emissivity enrichment (✓ FDR-significant).

| Site | Deposit system | n | Enrichment | q | Clustering | Night | AG100 |
|---|---|---|---|---|---|---|---|
| Courtland-Gleeson, AZ | Skarn/CRD | 255 | 1.52× | 0.001 | ✓ | 1.50× ✓ | 1.62× ✓ |
| Pioche, NV | Skarn/CRD | 274 | 1.44× | 0.003 | ✓ | 1.54× ✓ | 2.11× ✓ |
| Providence Mountains, CA | Skarn/CRD | 300 | 1.32× | 0.003 | ✓ | 1.19× ✓ | 1.31× ✓ |
| Eureka, NV | Skarn/CRD | 187 | 1.84× | < 0.001 | ✓ | 1.25× | 2.16× ✓ |
| Shafter, TX | Skarn/CRD | 49 | 3.04× | < 0.001 | ✓ | 1.44× | 2.53× ✓ |
| Sierrita, AZ | Porphyry | 265 | 1.99× | < 0.001 | ✓ | 0.86× | 1.97× ✓ |
| Lordsburg, NM | Porphyry | 275 | 1.72× | < 0.001 | – | 1.59× ✓ | 1.44× ✓ |
| Tombstone, AZ | Skarn/CRD | 201 | 1.54× | 0.013 | – | 1.52× ✓ | 2.30× ✓ |
| San Francisco (Frisco), UT | Skarn/CRD | 245 | 1.94× | < 0.001 | – | 1.02× | 1.94× ✓ |
| Tri-State, MO-KS-OK | Sediment-hosted (MVT) | 120 | 2.55× | < 0.001 | – | 1.31× | 2.55× ✓ |
| Victorio Mountains, NM | Skarn/CRD | 66 | 1.99× | 0.011 | – | 1.42× | 2.69× ✓ |
| Iron King, AZ | VMS | 212 | 1.73× | 0.002 | – | 1.51× ✓ | 0.98× |
| Tintic, UT | Skarn/CRD | 284 | 1.85× | < 0.001 | ✓ | 0.92× | 1.08× |
| Contact District, NV | Porphyry | 164 | 1.92× | 0.001 | ✓ | 0.88× | 0.57× |
| Silver Peak, NV | Mixed/Other | 146 | 1.43× | 0.031 | ✓ | 1.00× | 1.09× |
| Magdalena (Kelly), NM | Skarn/CRD | 69 | 1.64× | 0.041 | ✓ | 1.25× | 0.56× |
| Darwin, CA | Skarn/CRD | 206 | 1.50× | 0.004 | – | 1.08× | 1.17× |
| Mountain Pass, CA | Alkaline/Carbonatite | 343 | 1.21× | 0.026 | – | 1.14× | 1.08× |
| McDermitt Caldera, NV/OR | Mixed/Other | 36 | 4.49× | < 0.001 | – | 0.38× | 1.26× |

**Reading the table.** The three sites that pass all three replication checks (Courtland-Gleeson,
Pioche, Providence Mountains) are arid carbonate-replacement districts.

The finer 8-class SGMC lithology null ([`data_verification.md`](data_verification.md) §3.5) splits
them:

- **Providence Mountains survives (1.33×).** It is the only site that passes every replication *and*
  verification test (§3.5a).
- **Courtland-Gleeson and Pioche fall to ≤ 1×.** Their zones are repeatable, but they map the
  carbonate host package rather than anything beyond it.

McDermitt has the highest single-scene enrichment but does not replicate at night, and AG100 is not
significant there. Treat it as acquisition-specific.

### Changes from the radiance results

- **Still significant:** contact_district, courtland_gleeson, eureka, lordsburg,
  magdalena_kelly, providence_mountains, shafter, tombstone and tri_state.
- **Lost:** bisbee, christmas, crooks_gap, eagle_mountain, holden, karnes_county, lake_valley,
  lisbon_valley, mascot_jefferson_city, organ_mountains and uravan. The three "unexpected" robust
  sites of the radiance run (crooks_gap, karnes_county and holden: roll-front uranium and VMS) are
  gone. So is Bisbee, which the stratified null removes; it is still AG100-significant.
- **New:** darwin, iron_king, mcdermitt, mountain_pass, pioche, san_francisco_frisco, sierrita,
  silver_peak, tintic and victorio.

---

## Deposit-system results

Category median enrichment (observed ÷ stratified expected), 95% bootstrap CI, and FDR survival
rate (`results/phase6_deposit_system_categories.csv`):

| Category | n | Median enrichment | 95% CI | FDR-significant |
|---|---|---|---|---|
| Mafic / Ultramafic | 7 | 1.59× | 1.26 – 1.96 | 1 / 7 |
| Mixed / Other | 6 | 1.56× | 1.01 – 3.53 | 3 / 6 |
| **Skarn / Carbonate Replacement** | 28 | **1.44×** | **1.00 – 1.63** | **14 / 28** |
| Porphyry Cu-Mo-Au | 17 | 1.22× | 0.64 – 1.28 | 6 / 17 |
| Volcanogenic / VMS | 9 | 1.04× | 0.00 – 1.73 | 2 / 9 |
| Epithermal | 11 | 0.90× | 0.61 – 1.03 | 2 / 11 |
| Sediment-hosted | 12 | 0.88× | 0.00 – 1.16 | 2 / 12 |
| Uranium / Energy | 9 | 0.86× | 0.70 – 0.98 | 0 / 9 |
| Alkaline / Carbonatite | 10 | 0.64× | 0.00 – 1.21 | 2 / 10 |

(FDR counts in this table include the continuous-score test.)

**Skarn / carbonate-replacement still holds half of all significant sites**, but its CI now
reaches 1.0. The mafic/ultramafic median of 1.59× rests on few testable sites: vegetation masking
leaves only 2 of the 6 expansion mafic sites with ≥ 5 deposits.

**Pre-registered expectations:**

- **Held:** epithermal with exposed caps, 0 / 4; Carlin-type sediment-hosted, 0 / 6
  (0.76–1.16×); alkaline/carbonatite, 0 / 4; uranium, 0 / 9. The radiance-era roll-front
  "exceptions" (crooks_gap, karnes_county) are gone.
- **No longer testable:** the humid-terrain skarn control (cornwall_pa, french_creek_pa,
  dillsburg_pa, willsboro_lewis_ny, tilly_foster_ny, snoqualmie_wa). AG100 classes 77–98% of each
  footprint as vegetation or water. After masking, 5 of the 6 have ≤ 1 critical deposit, and
  Snoqualmie has 14 deposits with 2 hits. The radiance-era "0 / 6" came from classifying
  vegetated pixels; TIR has no mineral information there. 19 sites in total have < 5 deposits
  after masking, mostly eastern.
- **Mixed:** the 5 expansion porphyry districts give 1 / 5 (contact_district).

---

## Discovery-bias stratification

To test whether the spatial association could reflect historical ascertainment bias (deposits
found by the same visible alteration ASTER detects), critical-mineral deposits with a usable
discovery / first-production year were split at 1950 and each cohort tested against its site's own
null, unconditionally across the 87 sites with dated records.

Result: **both cohorts are enriched by similar amounts** (pre-1950 1.39×, *p* < 0.001; post-1950
1.32×, *p* = 0.005). A pure ascertainment effect, where deposits were found *because* of visible
alteration, would predict a stronger pre-1950 signal; the difference here is small. This is still
not exculpatory. Only ~14% of records have a usable year, the dated subsample is skewed towards
historic producers, and post-1950 exploration still built on earlier surface mapping. It is
reported as strictly secondary.

---

## Limitations

1. **The method is not demonstrated to be alteration-specific.** The primary null removes the
   bedrock/cover contrast, but not host-rock preference within bedrock. Carbonate and quartz
   indices respond to the carbonate and quartzite packages that host skarn/CRD deposits, so the
   zones can mark the host units without detecting their alteration (see `data_verification.md`
   §3.5). The defensible claim is host-package targeting in arid terrain.

2. **Single-scene results are acquisition-sensitive.** Night and day emissivity agree on site
   enrichment only moderately (ρ 0.35), pixel-level QI and CI agreement is near zero, and 11 of the
   radiance-era FDR sites did not reproduce on emissivity. Only the 3 sites passing every
   replication check should be cited as robust.

3. **ASTER SWIR bands are unavailable** in LP DAAC v004, so argillic / phyllic alteration
   (relevant to Carlin-type gold, high-sulfidation epithermal, and some porphyry systems) is
   invisible to the TIR method regardless of any statistical correction. That detection-gap
   hypothesis is now *testable*, not closed: see [GEMx / AVIRIS coverage](#gemx--aviris-coverage).
   It has not yet been tested. The below-chance epithermal, porphyry, and sediment-hosted medians
   may still reflect the gap — or they may not.

4. **Arid / semi-arid applicability only.** Vegetated ground is masked (AG100 NDVI > 0.3), so humid
   sites are mostly untestable rather than null.

5. **Small per-category sample sizes** (7–28 sites) mean the deposit-system comparison supports a
   directional ranking, not precise category-level estimates.

6. **Category assignment is not fully blind to outcome.** The 58 expansion-site categories and
   expectations were fixed before results; the original 51 were not.

7. **BH-FDR's dependence assumption** was not stress-tested against the more conservative
   Benjamini–Yekutieli variant.

---

## GEMx / AVIRIS coverage

This is a **coverage inventory**, not a mineral-identification result. USGS Earth MRI's GEMx
campaign (Geological Earth Mapping Experiment, PIs Kokaly and Hoefen, 2023–) flies NASA AVIRIS-C
and AVIRIS-5 over the arid Southwest to map critical-mineral alteration, including materials in
historic mine waste. Granule footprints from NASA's CMR were intersected with each of the 109
site bboxes (`scripts/gemx_coverage_report.py`). AVIRIS collections are facility-wide archives,
so only granules from the GEMx epoch (2023-04-25 onward) count. Coverage fraction is the union of
those footprints over the bbox in an equal-area projection, not a granule count.

Figures: [`figures/gemx_coverage_map.png`](../figures/gemx_coverage_map.png) (geography) and
[`figures/gemx_coverage_by_system.png`](../figures/gemx_coverage_by_system.png) (deposit system).

| | Sites |
|---|---|
| Any GEMx-era SWIR (AVIRIS-C or AVIRIS-5) | **71 / 109** |
| ≥75% of the bbox covered | 66 |
| Partial (25–75%) | 5 |
| None | 38 |

**What the map shows.** Coverage is a Southwest phenomenon. Green sites cluster in California,
Nevada, Arizona, and New Mexico — the GEMx flight boxes. Grey sites are the humid-terrain and
eastern control districts. That geographic split is the same split that already structures the
TIR result.

**The SWIR-gap test is now possible.** Limitation 2 said the below-chance epithermal and
porphyry medians "may partly reflect" missing SWIR. Those two categories are the ones GEMx
covers best:

| Category | Covered | Notes |
|---|---|---|
| Porphyry Cu-Mo-Au | 15 / 17 | Uncovered: `climax` (CO), `bingham_canyon` (UT) |
| Epithermal | 9 / 11 | Uncovered: `summitville` (CO), `republic` (WA) |
| **SWIR-gap cohort** | **24 / 28** | The sites that decide whether missing SWIR explains the TIR nulls |

**Sediment-hosted splits on the same geography.** All six Nevada Carlin-type sites
(`carlin_trend`, `cortez`, `getchell`, `jerritt_canyon`, `alligator_ridge`, `long_canyon`) have
~100% coverage. The six with none are humid-terrain MVT and related districts (Tri-State,
Viburnum, Mascot–Jefferson City, Austinville, Upper Mississippi Valley, Black Butte). A SWIR
test can speak to Carlin-type argillic alteration in Nevada. It cannot speak to the Midcontinent
MVT nulls.

**The emissivity primary sites are well covered.** 17 of the 19 FDR-significant sites have
≥ 75% GEMx-era SWIR coverage, including all three that pass every replication check. Only
`tintic` (UT) and `tri_state` (MO-KS-OK) have none. So GEMx can test whether the skarn/CRD
detections are alteration-specific or host-package-only (Limitation 1). (The radiance-era
"unexplained" sites holden, crooks_gap and karnes_county had zero coverage; they are no longer
significant.)

**Goldfield + Marysvale + Bodie SWIR pilots (band-depth indices on the
ASTER study box, not Tetracorder).** `scripts/gemx_swir_pilot.py --site {id}`
mosaics AVIRIS-5 L2A scenes onto `site.bbox_wgs84` (Goldfield and Bodie
EPSG:32611; Marysvale 32612). Admission is a greedy bbox set-cover (ASTER
`--mosaic` analogue:
paint the box). A single AVIRIS strip cannot pass ASTER's 30% gate, so lines
are added until the next one would add <2% new area; town and the science
target are seeded first.

| Site | Scenes | SWIR % of box | TIR % of box (valid pixels) | Overlap |
|---|---|---|---|---|
| Goldfield | 10 | 69% | 82% | 59% |
| Marysvale | 13 | 72% | (fig03 prints 61% of the unclipped scene) | 171k 90 m cells |
| Bodie | 13 | 79% | (cls from existing ratio mosaic) | 292k 90 m cells |

The 62% TIR label on Goldfield fig03 is the valid fraction of the
**unclipped ASTER scene raster**, which continues south of the study box.

**Correspondence** (`scripts/gemx_swir_correspondence.py`): on cells where
*both* sensors are valid, TIR strong (combined score ≥ 3) vs SWIR high
(depth ≥ overlap p90). Complementary hypothesis — the maps do not occupy
the same pixels — holds at n=2 HS plus one LS control.

| Site | Alunite Jaccard | Lift | Cliff’s δ | Pins |
|---|---|---|---|---|
| Goldfield | 0.075 | 1.83× | 0.15 | GNIS town is **not** SWIR-high on the live 90 m grid (depth 0.009) |
| Marysvale | 0.010 | 0.59× | −0.10 | Ridge SWIR-high, not TIR-strong; town TIR-strong, no alunite |
| Bodie (LS control) | 0.020 | 0.58× | −0.18 | Town 0.004 / Bluff 0.009; neither high, neither TIR-strong |

Goldfield in-box depths: Al-OH p50/p90/% >0.10 = 0.039 / 0.097 / 9.3%;
kaolinite-like 0.009 / 0.042 / 1.7%; alunite-like 0.014 / 0.051 / 2.3%.
Marysvale generic Al-OH **failed** (p50 −0.10; vegetation/continuum).
Doublets still fire but are rare (~0.8 M of 23 M pixels).

**Cutoff sensitivity** (`scripts/gemx_swir_sensitivity.py`), pre-registered
before the sweep: TIR score ≥ {2, 3, 4} × SWIR overlap percentile {80, 90,
95}. Jaccard stayed below 0.20 at all 18 cells (Goldfield max 0.135 at
score ≥ 2 / p80; Marysvale max 0.060). No pin was both TIR-strong and
SWIR-high. The published Goldfield-town pin story does **not** hold on
the live 90 m grid: GNIS town depth is 0.009, below even p80 (0.027).
The 0.137 value in an older `goldfield_correspondence.json` is stale.
Marysvale lift **changes sign** at the looser TIR cut (1.60–1.78× at
score ≥ 2 vs 0.45–0.63× at score ≥ 3). Complementary (low Jaccard) is
not an artifact of the published cuts. The lift direction at Marysvale
is. Figure: [`gemx_swir_sensitivity.png`](../figures/gemx_swir_sensitivity.png).

**Bodie low-sulfidation control** (prediction written before download:
alunite-like stays rare at the LS district; if it lights up at
Goldfield-like rates, the doublet is not HS-specific). GNIS town and
Bodie Bluff pins stayed quiet. Whole-box deep alunite (>0.10) is 0.49%
of the mosaic vs 1.26% at Goldfield; a positive doublet exists on 11%
of Bodie pixels vs 54% at Goldfield. The index is not globally silent —
scattered speckles remain — and it is not a Goldfield-style district
blob. Generic Al–OH median is negative (−0.035), Marysvale-like.
Complementary still holds (Jaccard 0.020). This does **not** make the
doublet HS-specific. n=2 HS + 1 LS control is not a 24-site test.

Figures: [`gemx_goldfield_aloh.png`](../figures/gemx_goldfield_aloh.png),
[`gemx_goldfield_correspondence.png`](../figures/gemx_goldfield_correspondence.png),
[`gemx_marysvale_aloh.png`](../figures/gemx_marysvale_aloh.png),
[`gemx_marysvale_correspondence.png`](../figures/gemx_marysvale_correspondence.png),
[`gemx_bodie_aloh.png`](../figures/gemx_bodie_aloh.png),
[`gemx_bodie_correspondence.png`](../figures/gemx_bodie_correspondence.png),
[`gemx_swir_sensitivity.png`](../figures/gemx_swir_sensitivity.png).
Stats: `results/gemx_*_pilot.json`, `results/*_correspondence.json`,
`results/gemx_swir_sensitivity.json`.

**What this is not.** USGS has not released derived mineral maps. Tetracorder /
MICA has not been run. Complementary means the two sensors see different ground,
not that SWIR “explains” the epithermal TIR nulls.

---

## Google Earth Engine parity (8-site panel)

Not a 109-site rewrite. A reproducibility module
([`docs/gee_parity.md`](gee_parity.md), `gee_parity.py`) ports the same six
TIR steps — scene selection, B10–B14, ratios, bbox clip, 70/90 percentiles,
score ≥ 3 — to `ASTER/AST_L1T_003` and measures agreement against the local
rasterio path.

Default run is **offline probes** that isolate the choices a cloud-native
port actually has to make (clip order, histogram vs exact percentile,
DN vs radiance UCC, GEE's default EPSG:4326 `reduceRegion` CRS). Live
`ee` is optional (`gee-parity --gee`). Figure:
[`gee_parity_agreement.png`](../figures/gee_parity_agreement.png).
