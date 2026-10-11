# Cloud masking and night-time scenes in ASTER TIR — methods review

Literature and data-access review compiled 2026-09-27. Scope: can published methods fix the two known mosaic problems in this pipeline (cloud-blended mosaics at 51/109 sites; 295/780 granules are night acquisitions, often mixed with day scenes)? Every reference below was checked against Crossref or the publisher page. Items I could not verify are marked **[unverified]**.

**Implemented (2026-10-09).** Recommendations adopted: AST_05 v004 day emissivity, QA_DataPlane
per-pixel mask with QA-aware reselection, Ninomiya QI/CI/MI, priority stacking, night AST_05 as a
separate replication layer, and AG100 as a multi-year benchmark. `MOD35_L2` was not needed. Results
are in `docs/data_verification.md` §2.5–2.6.

## Bottom line

- **The archive is closed and on-demand ordering is gone.** The ASTER TIR subsystem was permanently switched off on 2026-01-16 (per the LP DAAC CMR collection metadata). The LP DAAC on-demand system (S4PM) was decommissioned on 2025-12-15, and all V003 products are deprecated. The **V004 "final collections" are pre-computed for the whole mission** and delivered as per-band COGs on Earthdata Cloud, so **no ordering is needed**: `AST_05` v004 (emissivity), `AST_08` v004 (kinetic temperature), `AST_09T` v004 (surface radiance) and `AST_L1T` v004 are all directly searchable and downloadable with `earthaccess`.
- **Radiance ratios are temperature-dependent at a level comparable to the lithologic signal.** A purely blackbody (Planck) calculation at the ASTER band centres gives:
  - B13/B14: +0.084 %/K, so a 30 K day–night difference shifts it by about +2.5 %.
  - B13/B12 and B12/B13: ∓0.25 %/K, so 30 K shifts them by about 7–8 %.

  For comparison, Ninomiya's published Carbonate Index (CI) stretch spans only 1.005–1.055. Mixing day and night radiance in one percentile classification is therefore not defensible, and I found no published study that does it.
- **Emissivity (from temperature/emissivity separation, TES) is the published route to temperature-independent indices.** In principle it makes day and night comparable. In practice, night AST_05 retrievals are measurably less precise (Mushkin et al. 2020), so the safe choice is still daytime scenes, with night used only as a separate layer.
- **For a TIR-only, pixel-level cloud mask, the practical options are:**
  1. the QA_DataPlane shipped with V004 `AST_05`/`AST_08`;
  2. the coincident Terra MODIS cloud mask `MOD35_L2` v6.1, which is exactly what the ASTER project itself has used since 2009;
  3. a thermal "too cold" test against a reference temperature (a heuristic, not a published ASTER algorithm).

## 1. Pixel-level cloud detection for ASTER TIR

| Approach | Works at night / TIR-only? | Practical in Python for L1T? |
|---|---|---|
| Hulley & Hook (2008), two-pass Landsat-7-style algorithm with MODIS-derived cirrus and shadow tests | **No.** The pass-1 tests use VNIR, SWIR B4 and TIR, so it is daytime-only. It is degraded after April 2008 (SWIR loss); the ASTER GED team stopped adding post-2008 data "due to problems with cloud masking without the SWIR bands" (GED v3 User Guide). | Partially. An open-source pass-1 implementation exists (`Mining-for-the-Future/ASTER_preprocessing`, `aster_cloud_mask`), but it needs VNIR (and SWIR). |
| ASTER L1 ACCAA and scene `CloudCover` (the metadata the pipeline already uses) | Scene-level only. CC has been reassessed from MOD35 since June 2009. | Already in use. It is a scene percentage, not a mask. |
| **Coincident Terra MODIS `MOD35_L2` v6.1** (1 km; day and night). Tonooka & Tachikawa (2019) evaluate MOD35-derived ASTER cloud masks for 99.4 % of the archive. | **Yes.** It is simultaneous (same platform) and has night tests. Caveats: 12 % of night scenes have ≥ 30 % "uncertain" pixels, versus 2.2 % by day; about 2 % of masks show mislabelling or artefacts. C6.1 is recommended over C5/C6. | **Yes, medium effort.** Read with `satpy` (`modis_l2` reader) or `pyhdf`, geolocate with MOD03, resample to the L1T grid with `pyresample`, and buffer 1–2 km because 1 km cells are coarser than 90 m pixels. |
| **V004 `AST_05`/`AST_08` `QA_DataPlane` COGs.** The V4 User Guide says "each pixel will be marked as Good, Bad, Suspect and its cloudiness is noted." | Yes (it ships with night granules too). | **Yes, low effort.** It is one extra COG per granule. **[unverified: exact bit layout; see the ASTER Higher-Level Data Product QA Plan, JPL D-13841. Inspect one file before relying on it.]** |
| Ibaraki University "ASTER Cloud Mask Database" (per-scene MOD35 masks; http://tonolab.cis.ibaraki.ac.jp/ASTER/cloud/) | Yes. | Web form only; the site is live, but I found no API. |
| Thermal "cold-pixel" screen: flag B13 brightness temperature ≫ 10 K below a reference, such as the AG100 `Temperature Mean` or same-day MODIS LST | Yes. | Easy. **It is a heuristic, not a published ASTER method**, so thresholds would need tuning. Useful because Ninomiya & Fu (2016) note that **clouds and mist produce high B13/B14 (CI)**, which would inflate this pipeline's "silica" ratio. |
| Landsat Fmask (Zhu & Woodcock 2012) | No. It relies on VNIR/SWIR spectral tests. | Not transferable. **I found no published Fmask adaptation for ASTER.** |

## 2. ASTER Global Emissivity Dataset (GED): a cloud-free composite

**What it is.** Hulley et al. (2015), GRL. `AG100` v003 (100 m) and `AG1km` v003 (1 km) are the mean emissivity for **B10–B14** plus standard deviation, mean LST and its standard deviation, mean NDVI, a land/water mask, a resampled GDEM and an **observation count**. They come as 1°×1° HDF5 tiles, built by stacking **all clear-sky ASTER scenes from 2000 to 2008**. Processing used TES with Water Vapour Scaling (WVS) atmospheric correction and the Hulley & Hook (2008) cloud mask, and scenes with more than about 5 cm of precipitable water were discarded. `AG5KMMOH` v041 is the monthly 0.05° variant (too coarse for this work). Accuracy is about 1 % absolute band error at validated desert sites. Its predecessor, NAALSED (Hulley & Hook 2009), provided separate summer and winter means for North America.

**Access.**
- `earthaccess.search_data(short_name="AG100", version="003", bounding_box=...)` works; the tiles are cloud-hosted on LPCLOUD. For example, the McDermitt bbox returns 4 tiles.
- `AG1km` / `003` is available the same way.
- AppEEARS does **not** serve AG100/AST_05/AST_08. Its ASTER holdings are only the GDEM and water-body products.

**Geologic use.** The GED authors note that decorrelation stretches of GEDv3 "are easier to interpret for geologic studies." Ninomiya (2019, IGARSS) computed QI/CI/MI globally on GED and concluded the global map shows "the rough trend of lithology… however, the quality is not sufficient," recommending it be combined with scene-based regional maps. I found no peer-reviewed alteration-mapping study built solely on AG100. One Goldschmidt 2021 abstract (Ito et al.) maps silica content from GED.

**Pros for this pipeline:**
- Removes both problems at once: it is cloud-screened, and it is emissivity rather than radiance, so it is temperature-independent.
- Uniform processing across all 109 sites, which allows **fixed (absolute) index thresholds and cross-site comparison**, something the current scene-relative percentiles explicitly cannot offer.
- Tiny download and no mosaicking.
- The observation-count layer provides a per-pixel confidence weight.

**Cons:**
- 100 m grid, and averaging over 2000–2008 blurs features.
- Vegetation and soil-moisture mixing is baked into the mean.
- Low observation counts (< 3) produce artefacts (User Guide).
- Ninomiya found it inferior to curated single scenes.
- It is not the per-scene product the rest of the pipeline and its provenance assume.

**Assessment:** it is best adopted as a **parallel benchmark layer** (and a candidate replacement if enrichment results hold up), not as a silent swap.

## 3. Temperature/emissivity separation (TES) and why radiance ratios drift

**Why radiance ratios drift.** At-sensor radiance is roughly ε·B(λ,T) plus atmospheric terms, so

  ln(Lᵢ/Lⱼ) ≈ ln(εᵢ/εⱼ) + ln(B(λᵢ,T)/B(λⱼ,T)).

The second term depends on T approximately as c₂/T²·(1/λᵢ − 1/λⱼ), so the further apart the two bands, the larger the drift. That is why B12/B13-type ratios (9.1 vs 10.7 µm) drift about three times faster than B13/B14. On top of this, atmospheric water vapour affects B10–B12 and B14 differently, and daytime slope/aspect heating adds within-scene temperature variance.

**Published transforms, from cheapest to most rigorous:**
1. **Ninomiya brightness-temperature normalisation.** Ninomiya et al. (2005) and Ninomiya & Fu (2016) rescale each band's radiance so that the B13 brightness temperature equals **300 K**, applying the Planck function per pixel, and then compute the indices. They show QI and MI are fairly temperature-stable in good atmospheres, but **CI (B13/B14) "exhibits a strong dependency on the surface temperature,"** which this normalisation largely removes. It is cheap (pure numpy on L1T radiance, with no atmospheric correction), but it is only a "type of temperature separation" and not full TES.
2. **Emissivity normalisation (NEM)** and **alpha residuals.** Kealy & Hook (1993) and Hook et al. (1992) are simple TES variants: assume a fixed ε_max (for example 0.96), or remove the temperature term via Wien's approximation. Both need surface-leaving radiance, meaning atmospherically corrected input such as `AST_09T` v004.
3. **Full TES (Gillespie et al. 1998).** It combines NEM, band ratios and the maximum–minimum difference (MMD) empirical relation, with an iterative correction for reflected sky radiance. Expected accuracy is about ±1.5 K and ±0.015 emissivity. Residual errors are documented in Gillespie et al. (2011), and WVS correction (Tonooka 2005) reduces the atmospheric residuals. This is what `AST_05`/`AST_08` contain.

**Product status (checked in CMR, 2026-09-27):**

| Short name | Version | Status | Contents |
|---|---|---|---|
| `AST_05` | `004` | Live on LPCLOUD, 2000-03-04 → 2026-01-16 | `SRE_TIR_B10…B14.tif` plus `QA_DataPlane[2].tif` |
| `AST_08` | `004` | Live on LPCLOUD, same period | `SKT.tif` plus two QA planes |
| `AST_09T` | `004` | Live on LPCLOUD | Surface radiance TIR |
| `AST_L1T` | `004` | Live on LPCLOUD | Radiance COGs |
| `AST_05`, `AST_08` | `003` | **Deprecated 2025-12-15** | — |

- V004 Level-2 products use MERRA-2 atmospheric profiles (falling back to GDAS, then climatology) with MODTRAN; the TES algorithm itself is unchanged.
- V004 AST_05 includes **night granules**: the McDermitt bbox has 594 AST_05 granules in total, 334 of them flagged night.
- **[unverified]** Whether V004 applies WVS.

**Does emissivity make day and night comparable?** In principle yes, because emissivity is an intrinsic surface property. Empirically it is only partly true:
- Mushkin et al. (2020) found daytime AST_05 accurate to better than 0.01, but **night-time retrievals only to better than 0.03**. The differences concentrated in B10–B12, the quartz-diagnostic bands, and only the WVS correction brought night accuracy below 0.01.
- Dew adsorption raised apparent emissivity by up to 0.04 before dawn.
- Conversely, Kurata & Yamaguchi (2019) deliberately used **summer night** AST_05 at Cuprite and Goldfield because night data show less slope-orientation dependence and a better TES separation.

The conclusion is that night emissivity can be used, but should be kept as its own night-only layer, filtered to warm and dry nights, and not blended with day data.

## 4. Night-time ASTER TIR in geology

- **Thermal inertia and kinetic temperature.** This is the main published use. Examples:
  - Coolbaugh et al. (2007) used day and night AST_08 for geothermal anomalies at Brady's, Nevada.
  - Hewson et al. (2017) used night kinetic temperature and apparent thermal inertia (ATI) to separate outcrop from cover, noting that ASTER's 10:30/22:30 overpasses miss the diurnal extremes.
  - Nowicki et al. (2019) built Mojave LST mosaics for ATI.

  This is a *different* physical signal (density, porosity, induration), not mineralogy.
- **Night emissivity for lithology.** Kurata & Yamaguchi (2019) used summer night AST_05 (see §3).
- **Ninomiya indices.** The 2005 and 2016 papers select "atmospherically well-conditioned" scenes, check CI against multi-temporal scenes because water vapour and cloud mimic carbonate, and apply the 300 K normalisation. **[unverified]** Whether the 2005 paper explicitly compares day versus night index robustness; the full text is paywalled, and I read only the abstract, introduction and later summaries.
- **Rockwell & Hofstra (2008)** mapped quartz and carbonate across northern Nevada from 108 **daytime, May–September, "ostensibly cloud-free"** AST_05 scenes, chosen for maximum solar heating and no snow. Even so, they report index breaks at scene boundaries from uncorrected atmospheric and soil-moisture differences, and they had to mask cumulus manually from the carbonate map.
- **Mixing day and night in one mosaic:** I found no published precedent for lithologic or alteration indices.

## 5. Mosaicking best practice (ASTER TIR)

1. **Select scenes before mosaicking.** Match date, time, season and sun angle; choose low cloud; avoid snow and wet periods (Scheidt et al. 2008; Rockwell & Hofstra 2008 used May–September only).
2. **Compute indices per granule, then mosaic, rather than mosaicking raw radiance.** Use either histogram matching of index maps (Rockwell & Hofstra 2008) or linear gain/offset normalisation on pseudo-invariant features in the overlaps (Scheidt et al. 2008). Scheidt et al. found that normalising radiance before emissivity extraction worked, whereas normalising AST_05 emissivity afterwards "was not successful," and that B13/B14 normalisation failed because those bands have little variance.
3. **Use priority stacking (best scene on top), not blending** (Ninomiya & Fu 2016 rank each scene manually). Replace cloud gaps from another clear date rather than averaging over them; the GED User Guide warns that "filling in" cloud gaps produces artefacts.
4. **Prefer emissivity over radiance** when combining dates (Ninomiya & Fu 2016 show L1B-based and L2 emissivity-based index mosaics look "almost identical" when scenes are well chosen).
5. **Use orthorectified inputs.** L1T (L3A in Ninomiya's terms) avoids parallax misregistration.
6. **Expect residual seams.** Hewson et al. (2005) and Rockwell & Hofstra (2008) both document scene-boundary steps, so treat seam zones cautiously in the enrichment statistics.

## Recommended for this pipeline (ranked)

1. **Hard scene filter: daytime only, low cloud, same season** (already decided; hours of work). Use `AST_L1T` v004 with CMR `day_night_flag=day`, a granule `cloud_cover` cap (for example ≤ 10–20 %), and a May–September window. Never combine day and night granules in one site. Record the rejected granules and the reason in the provenance JSON.
2. **Ninomiya 300 K brightness-temperature normalisation of the L1T radiance before computing ratios** (hours). It is a published, cheap fix for the temperature term and matters most for B13/B14. It needs no new data.
3. **Per-pixel cloud mask** (0.5–1 day to wire in the QA plane; 2–4 days for MOD35). If the pipeline moves to item 4, use the AST_05/AST_08 V004 `QA_DataPlane` (after verifying the bit layout). If it stays on L1T radiance, use coincident `MOD35_L2` v6.1 (LAADS, cloud-hosted, same Terra overpass). Add the B13 cold-pixel screen as a cheap backstop.
4. **Switch the per-granule input from L1T radiance to `AST_05` v004 emissivity** (1–3 days). This is the Rockwell & Hofstra Nevada recipe: temperature-independent, with a QA mask included, pre-computed, and no ordering. Call `earthaccess.search_data(short_name="AST_05", version="004", bounding_box=..., temporal=...)` and filter to day. **[unverified: the `day_night_flag` / `cloud_cover` keyword passthrough in earthaccess was not tested live; the underlying CMR parameters were.]**
5. **AG100 v003 benchmark run** (1–2 days). Compute the same indices on GED mean emissivity for all 109 sites with `earthaccess` (`short_name="AG100"`, `version="003"`), weight by the observation-count layer, and compare deposit enrichment against the per-scene results. This is the only option that makes cross-site thresholds meaningful.
6. **Mosaic hygiene** (1 day). Mosaic per-granule indices with priority stacking and histogram or PIF matching in overlaps, and flag seam buffers.
7. **Optional, future: night data as a separate layer.** Use night `AST_08` v004 plus a day pair for apparent thermal inertia, or night-only `AST_05` for comparison. Do not feed night data into the alteration score.

**An issue this review surfaced in the current ratios (outside the brief, but it matters):** in the published literature, **B13/B14 is the Carbonate Index** (Ninomiya et al. 2005; Rockwell & Hofstra 2008; Kurata & Yamaguchi 2019). **B13/B12 is a quartz/felsic term** (it is part (b) of the Rockwell & Hofstra quartz index), and quartz proper is usually B11²/(B10·B12) (Ninomiya QI). The pipeline labels B13/B14 as "silica" and B13/B12 as "carbonate". In addition, `carbonate = B13/B12` and `mafic = B12/B13` in `spectral.alteration_ratios()` are **exact reciprocals**. Their per-scene percentile classes are therefore mutually exclusive, so together they can contribute at most 2 points to the 0–6 score, never 4. This deserves a separate review before any recalibration.

## Unverified or not found

- The bit layout of the V004 `QA_DataPlane` / `QA_DataPlane2` (JPL D-13841 QA Plan not retrieved).
- Whether V004 L2 uses WVS, and whether ASTER GED includes night scenes. Its cloud mask uses visible bands, which suggests a daytime bias, but I found no explicit statement.
- Whether Ninomiya et al. (2005) explicitly discusses night-time data (full text not accessed).
- No published Fmask adaptation for ASTER; no published lithologic mosaic that mixes day and night TIR.
- earthaccess `day_night_flag` / `cloud_cover` keyword passthrough (not executed; CMR parameters verified directly).

## References

All DOIs were checked via Crossref or the publisher page.

- Ackerman, S.A., Strabala, K.I., Menzel, W.P., Frey, R.A., Moeller, C.C., Gumley, L.E. (1998). Discriminating clear sky from clouds with MODIS. *J. Geophys. Res.* 103(D24), 32141–32157. https://doi.org/10.1029/1998JD200032
- Coolbaugh, M.F., Kratt, C., Fallacaro, A., Calvin, W.M., Taranik, J.V. (2007). Detection of geothermal anomalies using ASTER thermal infrared images at Bradys Hot Springs, Nevada, USA. *Remote Sens. Environ.* 106, 350–359. https://doi.org/10.1016/j.rse.2006.09.001
- Gillespie, A., Rokugawa, S., Matsunaga, T., Cothern, J.S., Hook, S., Kahle, A.B. (1998). A temperature and emissivity separation algorithm for ASTER images. *IEEE TGRS* 36(4), 1113–1126. https://doi.org/10.1109/36.700995
- Gillespie, A.R., Abbott, E.A., Gilson, L., Hulley, G., Jiménez-Muñoz, J.C., Sobrino, J.A. (2011). Residual errors in ASTER temperature and emissivity standard products AST08 and AST05. *Remote Sens. Environ.* 115, 3681–3694. https://doi.org/10.1016/j.rse.2011.09.007
- Hewson, R.D., Cudahy, T.J., Mizuhiko, S., Ueda, K., Mauger, A.J. (2005). Seamless geological map generation using ASTER in the Broken Hill–Curnamona province of Australia. *Remote Sens. Environ.* 99, 159–172. https://doi.org/10.1016/j.rse.2005.04.025
- Hewson, R.D., Robson, D., Carlton, A., Gilmore, P. (2017). Geological application of ASTER remote sensing within sparsely outcropping terrain, Central New South Wales, Australia. *Cogent Geoscience* 3, 1319259. https://doi.org/10.1080/23312041.2017.1319259
- Hook, S.J., Gabell, A.R., Green, A.A., Kealy, P.S. (1992). A comparison of techniques for extracting emissivity information from thermal infrared data for geologic studies. *Remote Sens. Environ.* 42, 123–135. https://doi.org/10.1016/0034-4257(92)90096-3
- Hulley, G.C., Hook, S.J. (2008). A new methodology for cloud detection and classification with ASTER data. *Geophys. Res. Lett.* 35, L16812. https://doi.org/10.1029/2008GL034644
- Hulley, G.C., Hook, S.J., Baldridge, A.M. (2008). ASTER land surface emissivity database of California and Nevada. *Geophys. Res. Lett.* 35, L13401. https://doi.org/10.1029/2008GL034507
- Hulley, G.C., Hook, S.J. (2009). The North American ASTER Land Surface Emissivity Database (NAALSED) Version 2.0. *Remote Sens. Environ.* 113, 1967–1975. https://doi.org/10.1016/j.rse.2009.05.005
- Hulley, G.C., Hook, S.J. (2011). Generating consistent land surface temperature and emissivity products between ASTER and MODIS data for Earth science research. *IEEE TGRS* 49(4), 1304–1315. https://doi.org/10.1109/TGRS.2010.2063034
- Hulley, G.C., Hook, S.J., Abbott, E., Malakar, N., Islam, T., Abrams, M. (2015). The ASTER Global Emissivity Dataset (ASTER GED): Mapping Earth's emissivity at 100 meter spatial scale. *Geophys. Res. Lett.* 42, 7966–7976. https://doi.org/10.1002/2015GL065564
- Kealy, P.S., Hook, S.J. (1993). Separating temperature and emissivity in thermal infrared multispectral scanner data: implications for recovering land surface temperatures. *IEEE TGRS* 31(6), 1155–1164. https://doi.org/10.1109/36.317447
- Kurata, K., Yamaguchi, Y. (2019). Integration and visualization of mineralogical and topographical information derived from ASTER and DEM data. *Remote Sens.* 11(2), 162. https://doi.org/10.3390/rs11020162
- Mushkin, A., Gillespie, A.R., Abbott, E.A., Batbaatar, J., Hulley, G., Tan, H., Tratt, D.M., Buckland, K.N. (2020). Validation of ASTER emissivity retrieval using the Mako airborne TIR imaging spectrometer at the Algodones Dune Field in Southern California, USA. *Remote Sens.* 12(5), 815. https://doi.org/10.3390/rs12050815
- Ninomiya, Y., Fu, B., Cudahy, T.J. (2005). Detecting lithology with ASTER multispectral thermal infrared "radiance-at-sensor" data. *Remote Sens. Environ.* 99, 127–139. https://doi.org/10.1016/j.rse.2005.06.009 (corrigendum: *RSE* 101(4), 567, 2006)
- Ninomiya, Y., Fu, B. (2016). Regional lithological mapping using ASTER-TIR data: case study for the Tibetan Plateau and the surrounding area. *Geosciences* 6(3), 39. https://doi.org/10.3390/geosciences6030039
- Ninomiya, Y. (2019). Global mapping of mineralo-lithological indices derived with ASTER multispectral thermal infrared data. *Proc. IGARSS 2019*, 5575–5578. https://doi.org/10.1109/IGARSS.2019.8900028
- Ninomiya, Y., Fu, B. (2019). Thermal infrared multispectral remote sensing of lithology and mineralogy based on spectral properties of materials. *Ore Geol. Rev.* 108, 54–72. https://doi.org/10.1016/j.oregeorev.2018.03.012
- Nowicki, S.A., Inman, D., Esque, T.C., Nussear, K.E., Edwards, C.S. (2019). Spatially consistent high-resolution land surface temperature mosaics for thermophysical mapping of the Mojave Desert. *Sensors* 19(12), 2669. https://doi.org/10.3390/s19122669
- Rockwell, B.W., Hofstra, A.H. (2008). Identification of quartz and carbonate minerals across northern Nevada using ASTER thermal infrared emissivity data — implications for geologic mapping and mineral resource investigations in well-studied and frontier areas. *Geosphere* 4(1), 218–246. https://doi.org/10.1130/GES00126.1
- Rowan, L.C., Mars, J.C. (2003). Lithologic mapping in the Mountain Pass, California area using ASTER data. *Remote Sens. Environ.* 84, 350–366. https://doi.org/10.1016/S0034-4257(02)00127-X
- Rowan, L.C., Mars, J.C., Simpson, C.J. (2005). Lithologic mapping of the Mordor, NT, Australia ultramafic complex by using ASTER. *Remote Sens. Environ.* 99, 105–126. https://doi.org/10.1016/j.rse.2004.11.021
- Scheidt, S., Ramsey, M., Lancaster, N. (2008). Radiometric normalization and image mosaic generation of ASTER thermal infrared data: an application to extensive sand sheets and dune fields. *Remote Sens. Environ.* 112, 920–933. https://doi.org/10.1016/j.rse.2007.06.020
- Tonooka, H. (2005). Accurate atmospheric correction of ASTER thermal infrared imagery using the WVS method. *IEEE TGRS* 43(12), 2778–2792. https://doi.org/10.1109/TGRS.2005.857886
- Tonooka, H., Tachikawa, T. (2019). ASTER cloud coverage assessment and mission operations analysis using Terra/MODIS cloud mask products. *Remote Sens.* 11(23), 2798. https://doi.org/10.3390/rs11232798
- Zhu, Z., Woodcock, C.E. (2012). Object-based cloud and cloud shadow detection in Landsat imagery. *Remote Sens. Environ.* 118, 83–94. https://doi.org/10.1016/j.rse.2011.10.028

**Data and product documentation:**
- LP DAAC (2025). *ASTER Version 4 Product User Guide* v1.1: https://lpdaac.usgs.gov/documents/2265/ASTER_User_Guide_V4_pcP80n5.pdf
- Hulley & Hook (2015). *ASTER GED User Guide* v3: https://lpdaac.usgs.gov/documents/120/ASTERGED_User_Guide_V3.pdf
- NASA Earthdata (2025). *ASTER On-Demand Processing To End Dec. 15, 2025*: https://www.earthdata.nasa.gov/data/alerts-outages/aster-demand-processing-end-dec.-15-2025
- Collection DOIs: `AST_05` v004 https://doi.org/10.5067/ASTER/AST_05.004 · `AST_L1T` v004 https://doi.org/10.5067/ASTER/AST_L1T.004
