# critical-minerals-aster

[![CI](https://github.com/nicole-m-aikin/critical-minerals-aster/actions/workflows/ci.yml/badge.svg)](https://github.com/nicole-m-aikin/critical-minerals-aster/actions/workflows/ci.yml)
![Python 3.11](https://img.shields.io/badge/python-3.11-blue)
![109 sites](https://img.shields.io/badge/sites-109-green)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

Systematic TIR alteration screening across **109 US critical-mineral districts**, validated against USGS MRDS deposit records under a corrected statistical framework: a **site-specific geometric null**, **BH-FDR** multiple-testing correction across both a binomial and a threshold-free continuous-score test, **spatial-clustering declustering** of MRDS records, and classification-threshold sensitivity analysis. 58 of the 109 sites were added in a pre-registered expansion — each carries a written detectability expectation (`results/deposit_system_categories.csv`) recorded before its imagery was processed.

**Zero-credential demo** (no EarthData account or rasters required):
```bash
git clone https://github.com/nicole-m-aikin/critical-minerals-aster.git
cd critical-minerals-aster && pip install -e .
python -m critical_minerals_aster demo
```

**Full results and figures:** [`docs/results.md`](docs/results.md) · [`figures/index.html`](figures/index.html)

---

## Scientific questions

1. Do ASTER-derived TIR alteration zones spatially correlate with known mineral occurrences in the USGS MRDS database?
2. Which deposit types are TIR-detectable, and which are systematically invisible?
3. Does the correlation vary by Earth MRI commodity category in geologically meaningful ways?
4. Are observed hit rates statistically significant above the null hypothesis of random deposit distribution?
5. Do non-critical bulk deposits (stone, sand, gravel) inflate headline hit rates — and is that signal geologically meaningful in its own right?

---

## Key findings (109 sites, corrected framework)

- **20 / 109** sites are FDR-significant on the primary binomial test; **31 / 109** on the continuous-score test (union 33).
- **13 sites** survive *every* check, including spatial-clustering declustering. **Nine are arid skarn / carbonate-replacement districts spanning AZ, NM, TX and CA** — Bisbee, Courtland-Gleeson, Tombstone, Christmas, Magdalena-Kelly, Organ Mountains, Lake Valley, Shafter, Eagle Mountain. The signal is a Southwest-wide regime, not one Arizona belt.
- **Skarn/Carbonate Replacement is the only deposit-system category whose median-enrichment 95% CI excludes 1.0** (1.64×, CI 1.02–2.25, n=28).
- A **systematic humid-terrain control** — 6 skarn districts of the same genesis in vegetated terrain (PA/NY/WA) — shows **0/6 FDR-significant**, supporting the vegetation/exposure hypothesis.
- **Exposed alteration caps do not rescue epithermal detectability**: 4 well-exposed-cap epithermal sites added, 0/4 significant; the category stays at median 0.44× (n=11).
- **Not yet demonstrated to be alteration-specific:** three sites — `crooks_gap`, `karnes_county` (roll-front uranium) and `holden` (VMS) — pass every filter despite deposit genesis that predicts no intrusion-related alteration halo (`holden` is a former mine with an exposed tailings field; its signal may be mine waste). These are 2/9 uranium and 1/9 VMS sites, not a category-wide effect. See [`docs/results.md`](docs/results.md#limitations).
- The **discovery-bias** stratification is unchanged: no era signal (post-1950 enrichment 0.92×, p=0.73); the original circularity claim is retracted.

The method acts as a **deposit-type selector, not a universal detector**. Two sites significant in an earlier 51-site version (Eureka, Lordsburg) are not in the robust set — both lose significance once spatially-clustered MRDS records are deduplicated (Eureka p=0.76, Lordsburg p=0.088 at the 1000 m radius).

---

## Data sources

| Dataset | Source | Notes |
|---|---|---|
| ASTER L1T (v004) | NASA EarthData / LP DAAC | TIR bands B10–B14, 90 m resolution |
| MRDS national deposit database | USGS mrdata.usgs.gov | deposit records across 109 bboxes |
| USGS Quaternary Faults | USGS QFAULTS REST API | Most sites |
| USGS SGMC fault data | USGS FeatureServer | All-age faults for select sites |

**Note on SWIR availability:** ASTER SWIR bands (B04–B09), standard for clay/argillic mapping, are not available in LP DAAC v004 for these areas. TIR bands (B10–B14, 8–12 µm) are used instead, suited for silica, carbonate, and mafic mineral mapping in arid terrain.

---

## Methods

### Band ratios

| Ratio | Formula | Target mineralogy |
|---|---|---|
| Silica/quartz | B13/B14 | Silicic alteration, rhyolite |
| Carbonate/dolomite | B13/B12 | Hydrothermal carbonate |
| Mafic | B12/B13 | Mafic volcanic rocks |

### Classification

Percentile-based thresholds (70th/90th) applied to the bbox-clipped ratio mosaic produce a 3-class map per ratio. An additive combined score (0–6) identifies pixels anomalous across multiple indicators. Strong anomaly zones (score ≥ 3) are vectorized to polygons. Thresholds are scene-relative; cross-site comparison of raw scores is not meaningful.

Multi-granule sites use feathered ratio mosaics with per-granule histogram normalization, classified once on the merged mosaic. This guarantees a 70/20/10 background/moderate/strong pixel distribution and avoids the per-granule max-merge inflation that biases results when N > 1 granules independently contribute classifications.

### Statistical significance

**Site-specific geometric null:** each site's chance hit rate is its own strong-anomaly-zone area ÷ valid-pixel TIR footprint area (not a pooled survey-wide rate). Two families of tests, FDR-corrected independently (Benjamini–Hochberg, α = 0.05): a one-sided exact **binomial** test and a threshold-free **Mann–Whitney U** test on the continuous 0–6 combined score. Robustness layers: footprint-aware Monte Carlo permutation (10,000 iters, seed 42), **DBSCAN spatial-clustering declustering** at 0/250/500/1000 m in each site's local metric CRS, and a classification-threshold sweep (65/85–75/95 percentiles; score cutoffs 2–4). MRDS deposits are clipped to the footprint polygon before counting. The full analysis chain is `scripts/regenerate_site_summaries.py` → `site_specific_null_significance.py` → `phase3…8` → `fig01–06`.

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
│       ├── spectral.py              # TIR I/O, granule selection, band ratios, mosaic tools
│       ├── classification.py        # percentile classification, vectorization
│       ├── metrics.py               # MRDS spatial join, per-site summary CSV
│       ├── mrds.py                  # MRDS CSV → GeoDataFrame, Earth MRI / mineral-system classifiers
│       ├── structure.py             # distance-to-fault annotation, buffer flags
│       ├── synthesis.py             # national summary CSV + figures (05, 06, 07)
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
- `contextily` — Esri World Imagery basemap tiles for figure 03 satellite panel
- `duckdb` — SQL-queryable national results
- `numpy` / `pandas` — array and tabular operations
- `matplotlib` — visualization

---

## Author

**Nicole Aikin** — MS Earth & Space Sciences, University of Washington (2025)
Metamorphic petrology · geochronology · ML pipelines for geoscience
[github.com/nicole-m-aikin](https://github.com/nicole-m-aikin)
