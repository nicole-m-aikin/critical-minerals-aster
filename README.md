# critical-minerals-aster

[![CI](https://github.com/nicole-m-aikin/critical-minerals-aster/actions/workflows/ci.yml/badge.svg)](https://github.com/nicole-m-aikin/critical-minerals-aster/actions/workflows/ci.yml)
![Python 3.11](https://img.shields.io/badge/python-3.11-blue)
![109 sites](https://img.shields.io/badge/sites-109-green)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

Systematic TIR alteration screening across **109 US critical-mineral districts** from ASTER surface emissivity, validated against USGS MRDS deposit records under a corrected statistical framework: a **bedrock/cover-stratified site-specific null**, **BH-FDR** multiple-testing correction across both a binomial and a threshold-free continuous-score test, **spatial-clustering declustering** of MRDS records, and classification-threshold sensitivity analysis. 58 of the 109 sites were added in a pre-registered expansion — each carries a written detectability expectation (`results/deposit_system_categories.csv`) recorded before its imagery was processed.

**Zero-credential demo** (no EarthData account or rasters required):
```bash
git clone https://github.com/nicole-m-aikin/critical-minerals-aster.git
cd critical-minerals-aster && pip install -e .
python -m critical_minerals_aster demo
```

**Full results and figures:** [`docs/results.md`](docs/results.md) · [`figures/index.html`](figures/index.html) · GEE parity: [`docs/gee_parity.md`](docs/gee_parity.md)

---

## Scientific questions

1. Do ASTER-derived TIR alteration zones spatially correlate with known mineral occurrences in the USGS MRDS database?
2. Which deposit types are TIR-detectable, and which are systematically invisible?
3. Does the correlation vary by Earth MRI commodity category in geologically meaningful ways?
4. Are observed hit rates statistically significant above the null hypothesis of random deposit distribution?
5. Do non-critical bulk deposits (stone, sand, gravel) inflate headline hit rates — and is that signal geologically meaningful in its own right?
6. Does a cloud-native Google Earth Engine implementation of the same TIR recipe reproduce the local rasterio pipeline? (8-site parity module — [`docs/gee_parity.md`](docs/gee_parity.md); not a 109-site rewrite.)

---

## Key findings (109 sites, AST_05 emissivity, bedrock/cover null)

The input is clear-sky, daytime **AST_05 surface emissivity** with per-pixel QA and AG100
water/vegetation masks, scored with Ninomiya quartz/carbonate/mafic indices. The primary null is
**bedrock-vs-cover stratified**. Earlier radiance-ratio results turned out to depend on acquisition
temperature ([`docs/data_verification.md`](docs/data_verification.md) §2.4–2.6).

- **19 / 108** sites are FDR-significant on the primary binomial test (23 under the geometric
  null); **10** of them are also robust at every declustering radius.
- **Three sites replicate on every independent check** (declustering, night AST_05, and the AG100
  multi-year mean emissivity): **Courtland-Gleeson (AZ), Pioche (NV) and Providence Mountains
  (CA)**, all arid carbonate-replacement districts. Five more pass two of the three: Eureka,
  Lordsburg, Shafter, Sierrita and Tombstone.
- **Only Providence Mountains also survives a finer 8-class lithology null.** Courtland-Gleeson
  and Pioche collapse to ≤ 1× under it, so their zones are repeatable maps of the carbonate host
  package rather than evidence of alteration beyond it.
- **14 of the 28 skarn/carbonate-replacement sites are FDR-significant**, more than any other
  deposit system, but the category's median-enrichment CI now reaches 1.0 (1.44×, 1.00–1.63).
- **Pre-registered nulls held:** exposed-cap epithermal 0/4, Carlin-type 0/6, uranium 0/9. The
  radiance-era "unexplained" roll-front and VMS sites (crooks_gap, karnes_county, holden) are no
  longer significant.
- **Humid-terrain sites are now untestable rather than null:** vegetation masking removes 77–98%
  of each humid skarn control.
- **Single-scene results are acquisition-sensitive:** day–night site enrichment agrees only
  moderately (ρ 0.35), and 11 of the 20 radiance-era FDR sites did not reproduce.
- **GEMx / AVIRIS coverage (inventory, not a mineral result):** 71 / 109 sites, including 17 of
  the 19 primary sites, so SWIR can test alteration-specificity there. See
  [`docs/results.md`](docs/results.md#gemx--aviris-coverage).
- **Discovery-era split:** both cohorts are enriched by similar amounts (pre-1950 1.39×, post-1950
  1.32×). This is secondary, with ~14% dated coverage.

The method acts as a **deposit-type selector and host-package targeting tool, not a demonstrated
alteration detector** ([`docs/results.md`](docs/results.md#limitations)).

---

## Data sources

| Dataset | Source | Notes |
|---|---|---|
| ASTER AST_05 (v004) | NASA EarthData / LP DAAC | TES surface emissivity B10–B14 + per-pixel QA, 90 m (day primary; night replication) |
| ASTER GED AG100 (v003) | NASA EarthData / LP DAAC | Water / vegetation masks; multi-year mean-emissivity benchmark |
| USGS SGMC geology | USGS ScienceBase | Bedrock / cover strata for the null |
| ASTER L1T (v004) | NASA EarthData / LP DAAC | Legacy radiance input; GEE parity module |
| MRDS national deposit database | USGS mrdata.usgs.gov | deposit records across 109 bboxes |
| USGS Quaternary Faults | USGS QFAULTS REST API | Most sites |
| USGS SGMC fault data | USGS FeatureServer | All-age faults for select sites |
| GEMx / AVIRIS / MASTER | NASA CMR / ORNL DAAC (Earth MRI) | Footprints at 71 / 109 sites; one Goldfield AVIRIS-5 scene has band-depth indices. No USGS derived mineral maps ingested. |

**Note on SWIR availability:** ASTER SWIR bands (B04–B09), standard for clay/argillic mapping, are not available in LP DAAC v004 for these areas. TIR bands (B10–B14, 8–12 µm) are used instead. USGS Earth MRI's GEMx campaign (AVIRIS-C / AVIRIS-5, 2023–) supplies the SWIR that ASTER cannot: 71 sites have usable flight-line coverage, including 24 of the 28 epithermal and porphyry nulls. A one-scene Goldfield pilot now maps Al-OH plus kaolinite-like and alunite-like doublets; that is a band-depth index, not Tetracorder, and not a 24-site test.

---

## Methods

### Input and indices

AST_05 v004 TES surface emissivity from clear-sky (≤ 30% cloud), daytime, May–September granules, never mixed with night. Per pixel, the QA plane must report good quality, clear sky and no adjacent cloud. AG100 v003 water and vegetation (mean NDVI > 0.3) are removed before classification.

| Index (slot name) | Formula | Target mineralogy |
|---|---|---|
| Quartz QI (`silica`) | B11² / (B10·B12) | Quartz-rich rock, silicification |
| Carbonate CI (`carbonate`) | B13 / B14 | Calcite / dolomite |
| Mafic MI (`mafic`) | B12 / B13 | Mafic / ultramafic rock |

(Ninomiya 2003; Ninomiya et al. 2005. The legacy radiance ratios are kept in `spectral.alteration_ratios` for the GEE parity module.)

### Classification

Percentile-based thresholds (70th/90th) applied to the bbox index mosaic produce a 3-class map per index. An additive combined score (0–6) identifies pixels anomalous across multiple indicators. Strong anomaly zones (score ≥ 3) are vectorized to polygons. Thresholds are scene-relative; cross-site comparison of raw scores is not meaningful.

Multi-granule sites compute indices per granule and priority-stack them onto one 90 m grid. The reference granule wins, and later granules only fill gaps after mean/std matching on the overlap; nothing is blended.

### Statistical significance

**Bedrock/cover stratified null (primary):** the footprint is split into bedrock, cover (SGMC "Unconsolidated") and unmapped. Each deposit's chance of a hit is the zone coverage of its own stratum, with an exact Poisson-binomial test. The **geometric null** (zone area ÷ valid-pixel footprint area) is reported alongside. Two families of tests, FDR-corrected independently (Benjamini–Hochberg, α = 0.05): a one-sided exact **binomial** test and a threshold-free **Mann–Whitney U** test on the continuous 0–6 combined score. Robustness layers: footprint-aware Monte Carlo permutation (10,000 iters, seed 42), **DBSCAN spatial-clustering declustering** at 0/250/500/1000 m in each site's local metric CRS, and a classification-threshold sweep (65/85–75/95 percentiles; score cutoffs 2–4). MRDS deposits are clipped to the footprint polygon before counting. The full analysis chain is `scripts/regenerate_site_summaries.py` → `site_specific_null_significance.py` → `phase3…8` → `fig01–06`.

**Pre-registration:** the 58 expansion sites each carry a written NULL / POSITIVE / POSITIVE-ISH / NULL-TO-WEAK expectation in `results/deposit_system_categories.csv` (`[EXPANSION +58]` rows), recorded before the pipeline was run.

---

## Reproducing this analysis

### 1. Clone and set up

```bash
git clone git@github.com:nicole-m-aikin/critical-minerals-aster.git
cd critical-minerals-aster
conda env create -f environment.yml
conda activate aster-minerals
pip install -e .
```

### 2. EarthData credentials

Create a free account at [urs.earthdata.nasa.gov](https://urs.earthdata.nasa.gov). The pipeline uses `earthaccess.login(strategy="netrc")` on first run; credentials are cached in `~/.netrc`.

### 3. Run the pipeline

```bash
# Single site (data must already be in data/sites/{id}/aster/)
python -m critical_minerals_aster run --site mcdermitt

# Download from EarthData, build mosaic, then process
python -m critical_minerals_aster run --site mcdermitt --mosaic

# All sites (skip existing outputs)
python -m critical_minerals_aster run-batch --all-sites --skip-existing

# Regenerate national summary + synthesis figures
python -m critical_minerals_aster synthesize

# GEMx / AVIRIS coverage inventory (CMR metadata; no EarthData login)
python scripts/gemx_coverage_report.py
python scripts/fig_gemx_coverage.py

# Goldfield Al-OH / kaolinite / alunite band-depth pilot (~2 GB AVIRIS-5 L2A)
python scripts/gemx_goldfield_pilot.py

# Local rasterio vs Google Earth Engine TIR parity (8-site panel; no GEE login)
python -m critical_minerals_aster gee-parity
python scripts/fig_gee_parity.py
# Live Earth Engine (optional): pip install -e ".[gee]" && earthengine authenticate
python -m critical_minerals_aster gee-parity --gee --project YOUR_GCP_PROJECT
```

### 4. Query results

```bash
python -c "
import duckdb
con = duckdb.connect('results/results.duckdb')
print(con.execute(\"\"\"
    SELECT site_id, hit_rate_pct, n_deposits_bbox, n_deposits_in_zones
    FROM site_summaries WHERE row_type='site'
    ORDER BY hit_rate_pct DESC
\"\"\").fetchdf())
"
```

---

## Repo structure

```
critical-minerals-aster/
├── sites/
│   ├── index.yaml                   # list of 109 site IDs
│   └── {site_id}.yaml               # bbox, granule, classification params, structure layers
├── src/
│   └── critical_minerals_aster/
│       ├── config.py                # SiteConfig, ClassificationParams, StructureLayer
│       ├── paths.py                 # SitePaths — all file/dir paths
│       ├── spectral.py              # TIR I/O, granule selection, AST_05 QA, Ninomiya indices
│       ├── ag100.py                 # ASTER GED AG100 water / NDVI masks
│       ├── lithology.py             # SGMC lithology + bedrock/cover strata for the null
│       ├── classification.py        # percentile classification, vectorization
│       ├── metrics.py               # MRDS spatial join, per-site summary CSV
│       ├── mrds.py                  # MRDS CSV → GeoDataFrame, Earth MRI / mineral-system classifiers
│       ├── structure.py             # distance-to-fault annotation, buffer flags
│       ├── synthesis.py             # national summary CSV + figures (05, 06, 07)
│       ├── hyperspectral.py         # GEMx / AVIRIS CMR coverage lookup
│       ├── gee_parity.py            # local rasterio vs GEE TIR parity (8-site panel)
│       ├── terrain.py               # hillshade DEM overlay
│       └── pipeline.py             # run_site() / run_batch() orchestration, mosaic builder
├── docs/
│   └── results.md                   # scientific findings, figures, interpretation
├── scripts/                         # standalone analysis scripts (significance, structure fetch)
├── results/                         # generated per-site CSVs + results.duckdb (gitignored)
├── figures/                         # synthesis figures + per-site gallery
├── data/                            # not committed (ASTER rasters, MRDS CSV, structure GeoJSONs)
├── tests/
├── environment.yml
├── pyproject.toml
└── README.md
```

---

## Dependencies

- `rasterio` — raster I/O, feature extraction, mosaic blending
- `geopandas` / `shapely` — vector operations and spatial joins
- `scipy` — binomial significance tests
- `earthaccess` — NASA EarthData authentication and download
- `earthengine-api` — optional; live GEE half of the TIR parity benchmark (`pip install -e ".[gee]"`)
- `contextily` — Esri World Imagery basemap tiles for figure 03 satellite panel
- `duckdb` — SQL-queryable national results
- `numpy` / `pandas` — array and tabular operations
- `matplotlib` — visualization

---

## Author

**Nicole Aikin** — MS Earth & Space Sciences, University of Washington (2025)
Metamorphic petrology · geochronology · ML pipelines for geoscience
[github.com/nicole-m-aikin](https://github.com/nicole-m-aikin)
