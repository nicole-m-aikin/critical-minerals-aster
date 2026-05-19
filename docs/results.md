# Results

ASTER TIR alteration mapping across 30 US critical mineral sites. Full figures: [`figures/index.html`](../figures/index.html).

---

## Site-level hit rates

**6 of 30 sites** show critical-mineral hit rates significantly above chance (binomial + spatial permutation, p < 0.05 on both tests):

| Site | Coverage | Critical hit rate | Binomial p | Permutation p | Deposit type |
|---|---|---|---|---|---|
| Bisbee, AZ | 7.3% | 26.1% | < 0.001 | < 0.001 | Skarn / VMS Cu-Zn |
| McDermitt, NV/OR | 8.3% | 21.6% | 0.010 | 0.010 | Caldera (Li, Hg, Au) |
| Mountain Pass, CA | 12.6% | 17.1% | 0.008 | 0.008 | Carbonatite REE |
| Hanover-Fierro, NM | 4.2% | 16.0% | 0.019 | 0.019 | Cu-Mo-Zn skarn |
| Yerington, NV | 5.3% | 12.8% | < 0.001 | < 0.001 | Porphyry Cu |
| Jerome, AZ | 5.7% | 11.6% | 0.003 | 0.003 | VMS Cu-Zn-Ag |

Hit rates use **critical-mineral deposits only** (Non-Critical earth_mri_category excluded — see [Non-critical inflation](#non-critical-inflation) below).

**10 sites are anti-correlated** (p ≈ 1 — deposits actively avoid anomaly zones): Climax CO, Goldfield NV, Oatman NV/AZ, Steamboat Springs NV, Marysvale UT, Bear Lodge WY, Pea Ridge MO, Bingham Canyon UT, Bagdad AZ, Silver Peak NV.

Anti-correlations are geologically coherent: Climax is a deep porphyry Mo with no surface alteration expression; Goldfield and Oatman are epithermal Au systems where the alteration cap is eroded or buried; Steamboat Springs is an active geothermal surface whose MRDS entries represent subsurface economics, not surface alteration.

---

## Non-critical inflation

Non-critical deposits (stone, sand, gravel, aggregate) have a **higher** pooled hit rate (9.8%) than critical minerals (7.0%). ASTER TIR detects silica and carbonate enrichment — exactly the mineralogy that makes rock commercially useful as aggregate. This inflates headline all-deposits hit rates at sites with dense bulk-mineral MRDS coverage.

**Green River, WY** is the clearest example: 52 of 53 deposits are non-critical, producing a 15.1% all-deposits hit rate (p = 0.005, appeared significant). Critical-only: 1 deposit, 0 hits, p = 1.0. The significance was illusory.

**Thacker Pass, NV** is a genuine borderline case: 19.2% critical hit rate but only 26 critical deposits, leaving it underpowered (critical-only p = 0.080 vs. all-deposits p = 0.0002).

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

## Interpretation limits

- **TIR-only:** SWIR clay/argillic mapping (B04–B09) is not available in LP DAAC v004 for these areas. The method maps silica, carbonate, and mafic mineralogy but misses argillic/phyllic alteration halos detectable with SWIR.
- **Scene-relative thresholds:** Percentile classification is per-scene; raw anomaly scores are not comparable across sites.
- **MRDS uncertainty:** Deposit locations are report-derived and may be offset from true outcrop or mineralized footprint.
- **Significance null model:** Both tests assume uniform deposit distribution within the bbox. Spatial clustering of real deposits makes p-values conservative — actual null distributions are not uniform.
- **Anti-correlations are informative:** p ≈ 1 means the method is physically incapable of detecting the dominant deposit type at that site, not that the zones are wrong.
