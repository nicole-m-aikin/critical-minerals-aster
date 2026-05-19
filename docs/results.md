# Results

ASTER TIR alteration mapping across 30 US critical mineral sites. Full figures: [`figures/index.html`](../figures/index.html).

---

## Site-level hit rates

**3 of 30 sites** show critical-mineral hit rates significantly above chance (binomial + spatial permutation, p < 0.05 on both tests):

| Site | Zone coverage¹ | Critical hit rate | Binomial p | Permutation p | Deposit type |
|---|---|---|---|---|---|
| Bisbee, AZ | 11.6% | 43.2% (41/95) | < 0.001 | < 0.001 | Skarn / VMS Cu-Zn |
| Jerome, AZ | 10.2% | 25.7% (19/74) | < 0.001 | < 0.001 | VMS Cu-Zn-Ag |
| Mountain Pass, CA | 12.6% | 17.4% (60/345) | 0.006 | 0.002 | Carbonatite REE |

¹ Zone coverage = strong anomaly zone area / TIR valid-pixel footprint area (not full site bbox).

Hit rates use **critical-mineral deposits only** (Non-Critical earth_mri_category excluded — see [Non-critical inflation](#non-critical-inflation) below). Deposits are clipped to the actual TIR valid-pixel footprint polygon before counting.

**Near-significant sites** (one test below 0.05, other borderline): McDermitt NV/OR (binom 0.052, perm 0.009), Thacker Pass NV (binom 0.075, perm 0.043), Yerington NV (binom 0.109, perm 0.0002), Hanover-Fierro NM (binom 0.117, perm 0.071). These sites have real signal but do not clear the dual-test threshold after footprint clipping tightened the null model.

**11 sites are anti-correlated** (p ≈ 1 — deposits actively avoid anomaly zones): Climax CO, Goldfield NV, Oatman AZ, Steamboat Springs NV, Marysvale UT, Bear Lodge WY, Pea Ridge MO, Bingham Canyon UT, Bagdad AZ, Elk Creek NE, Green River WY.

Anti-correlations are geologically coherent: Climax is a deep porphyry Mo with no surface alteration expression; Goldfield and Oatman are epithermal Au systems where the alteration cap is eroded or buried; Steamboat Springs is an active geothermal surface whose MRDS entries represent subsurface economics, not surface alteration. Green River has zero critical-mineral hits (all 53 MRDS records are non-critical); Elk Creek has one critical deposit and zero hits.

---

## Non-critical inflation

Non-critical deposits (stone, sand, gravel, aggregate) have a **higher** pooled hit rate (9.8%) than critical minerals (7.0%). ASTER TIR detects silica and carbonate enrichment — exactly the mineralogy that makes rock commercially useful as aggregate. This inflates headline all-deposits hit rates at sites with dense bulk-mineral MRDS coverage.

**Green River, WY** is the clearest example: 52 of 53 deposits are non-critical, producing a 15.1% all-deposits hit rate (p = 0.005, appeared significant). Critical-only: 1 deposit, 0 hits, p = 1.0. The significance was illusory.

**Thacker Pass, NV** is a genuine borderline case: 21.7% critical hit rate (5/23 deposits) but too few deposits for the dual-test threshold (critical-only binom p = 0.075, perm p = 0.043 vs. all-deposits binom p = 0.0002).

No sites gained significance when the Non-Critical filter was applied — confirming that aggregate deposits inflate rather than dilute.

**Aggregate quarry application:** The stone/sand/gravel correlation is a real economic signal for aggregate quarry siting. For that use case, the all-deposits and commodity-level rows in the CSVs (filtered to `earth_mri_category = 'Non-Critical'`) are the right lens. Green River, Silver Peak, and Marysvale are the strongest candidates.

---

## By Earth MRI category

| Category | Avg hit rate | TIR-detectable? | Notes |
|---|---|---|---|
| REE | ~27% | Yes | Carbonatite + vein-hosted REE co-spatial with silicic/carbonate alteration |
| Battery Metals – Co/Ni | ~13% | Yes | Mafic/ultramafic host rock shows mafic ratio anomalies |
| Base Metals | ~9% | Yes | Porphyry/skarn/VMS alteration halos |
| PGM | ~12% | Yes | Mafic ratio; small n |
| Gold/Silver | ~6% | Partial | Caldera/VMS-associated Au yes (Jerome, McDermitt); Carlin/placer no |
| Energy | ~4% | No | Roll-front U — no surface alteration expression |
| Industrial | ~5% | No | Structurally/sediment-hosted, no alteration |
| Non-Critical (aggregate) | ~10% | Yes — different application | See above |

No category reaches national significance individually; signal is site-specific rather than category-wide.

**Mineral system standouts:** Marine Chemocline and Magmatic REE systems have the highest pooled hit rates — both are associated with hydrothermally altered host rocks that ASTER TIR resolves well. Porphyry Cu-Mo-Au and skarn systems are moderate performers. Sediment-hosted (Carlin, MVT, roll-front) and placer systems score at or below chance.

---

## Structure proximity

The 6 significant sites span a wide range of mean deposit–fault distances (1–20 km on log scale), so structural proximity alone does not predict TIR detectability. Sites with the tightest deposit–fault clustering (Globe-Miami, Pea Ridge, Thacker Pass, mean ~1 km) are not the same set as those with the highest hit rates. The signal is more strongly driven by deposit type and host-rock alteration style than by mapped fault proximity.

---

## Figures

| Figure | Content |
|---|---|
| `figures/05_national_hit_rates.png` | Stacked bar by Earth MRI category, critical minerals only; sites sorted by critical hit rate |
| `figures/06_structure_hit_rate.png` | Log-scale scatter: mean fault distance vs critical-mineral hit rate; red = significant |
| `figures/07_hitrate_comparison.png` | Paired dot plot: all-deposits vs critical-only hit rate per site; exposes aggregate inflation |
| `figures/sites/{id}/00_composite_rgb.png` | False-color TIR composite |
| `figures/sites/{id}/01_tir_band_ratios.png` | Three band ratio maps with global colorbars |
| `figures/sites/{id}/02_classification.png` | Per-ratio classification + combined score |
| `figures/sites/{id}/03_deposit_overlay.png` | Anomaly zones, MRDS deposits, fault corridors, scale bar |
| `figures/sites/{id}/05_structure_proximity.png` | Strip chart: deposit distance to nearest fault by commodity group |

---

## Spatial coverage notes

ASTER TIR granules are ~60 × 60 km swaths at an oblique angle and do not always fill the full configured site bbox. MRDS deposits are clipped to the **valid-pixel footprint polygon** (the actual diagonal scene boundary, not the rectangular bbox) before counting, so deposit totals reflect only deposits within TIR coverage. The zone coverage fraction used in significance testing (zone area / footprint area) also uses this polygon. The TIR data boundary is shown as a dashed polygon outline on each figure 03.

**Sites with partial TIR coverage** (footprint polygon area / site bbox area):

| Site | TIR coverage | Notes |
|---|---|---|
| Green River, WY | ~50% of bbox | Granule covers roughly the SW half; 0 critical hits — confirmed false positive on all-deposits test |
| Bisbee, AZ | ~64% of bbox | Parallelogram ASTER swath; NE corner absent |
| Yerington, NV | ~55% of bbox | Diagonal upper-left cutoff |
| Gas Hills, WY | ~72% of bbox | Diagonal lower-right corner clipped |
| Pea Ridge, MO | ~86% of bbox | Diagonal upper-left cutoff; main deposit cluster within covered area |

**Bingham Canyon spatial mismatch:** The large TIR anomaly zone at this site (figure 03, upper-right cluster) corresponds to **Wasatch Front carbonate formations**, not the Bingham Canyon porphyry Cu-Mo-Au system. The actual mine (Oquirrh Mountains) is in the center-left of the bbox where MRDS deposit points cluster, but produces a weak TIR signal. The anti-correlation reflects geometric displacement between the dominant anomaly zone and the target deposit, not failure of alteration detection.

**Structural corridor density at Ozark sites:** Pea Ridge and Viburnum Trend show structural corridor buffers covering 40–50% of the map area, reflecting the density of mapped faults/lineaments in the Missouri Ozarks rather than real structural control on mineralization. The "% of deposits within 500 m of structure" metric is inflated at these sites.

---

## Interpretation limits

- **TIR-only:** SWIR clay/argillic mapping (B04–B09) is not available in LP DAAC v004 for these areas. The method maps silica, carbonate, and mafic mineralogy but misses argillic/phyllic alteration halos detectable with SWIR.
- **Scene-relative thresholds:** Percentile classification is per-scene; raw anomaly scores are not comparable across sites.
- **MRDS uncertainty:** Deposit locations are report-derived and may be offset from true outcrop or mineralized footprint.
- **Significance null model:** Both tests assume uniform deposit distribution within the bbox. Spatial clustering of real deposits makes p-values conservative — actual null distributions are not uniform.
- **Partial TIR coverage:** ASTER granule swaths do not always fill the full site bbox. At Green River (~50% covered), Bisbee (~64%), Yerington (~55%), Gas Hills (~72%), and Pea Ridge (~86%), deposit counts and significance tests reflect only the covered portion. Deposits and significance testing use the valid-pixel polygon (not a bounding rectangle); the TIR data boundary is shown on each figure 03.
- **Anti-correlations are informative:** p ≈ 1 means the method is physically incapable of detecting the dominant deposit type at that site, not that the zones are wrong.
