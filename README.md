# critical-minerals-aster

Spectral alteration mapping across 45 US critical mineral sites using ASTER thermal infrared (TIR) band ratio analysis. The pipeline identifies surface alteration zones, validates them against USGS MRDS mineral deposit data, and runs binomial significance tests to determine where — and for which deposit types — the method produces signal above chance.

**Results and figures:** [`docs/results.md`](docs/results.md) · [`figures/index.html`](figures/index.html)

---

## Scientific questions

1. Do ASTER-derived TIR alteration zones spatially correlate with known mineral occurrences in the USGS MRDS database?
2. Which deposit types are TIR-detectable, and which are systematically invisible?
3. Does the correlation vary by Earth MRI commodity category in geologically meaningful ways?
4. Are observed hit rates statistically significant above the null hypothesis of random deposit distribution?
5. Do non-critical bulk deposits (stone, sand, gravel) inflate headline hit rates — and is that signal geologically meaningful in its own right?

---

## Key findings (45 sites)

**12 sites** are significant on the critical-mineral-only binomial test (null = 8.5%); **3 additional sites** are significant on the all-deposit test only (non-critical or base-metals driven).

| Site | State | Deposit type | Critical hit rate | p |
|---|---|---|---|---|
| Bisbee | AZ | Skarn / carbonate-hosted Cu-Ag | 43.2% | < 0.001 |
| Eureka | NV | Pb-Zn-Au-Ag skarn / carbonate | 37.2% | < 0.001 |
| Magnet Cove | AR | Alkalic igneous complex REE/Ti | 23.1% | 0.019 |
| McDermitt Caldera | NV/OR | Li-Cs-REE sedimentary/caldera | 22.2% | 0.009 |
| Lordsburg | NM | Porphyry Cu-Mo | 21.8% | < 0.001 |
| Thacker Pass | NV | Lithium brine | 21.7% | 0.040 |
| Ely (Robinson) | NV | Porphyry Cu-Au | 20.8% | 0.048 |
| Sierrita | AZ | Porphyry Cu-Mo | 17.9% | 0.006 |
| Steamboat Springs | NV | Epithermal / geothermal Au-Ag | 14.5% | 0.003 |
| Yerington | NV | Porphyry Cu / skarn | 13.7% | 0.028 |
| Randsburg | CA | Epithermal Au-Ag-W | 11.8% | 0.004 |
| Mountain Pass | CA | Carbonatite REE | 11.5% | 0.038 |

The method acts as a **deposit-type selector, not a universal detector**: skarn/carbonate-hosted and arid porphyry systems produce the clearest signal; sediment-hosted (Carlin, MVT), placer, and vegetated-terrain deposits anti-correlate or score at chance. See [`docs/results.md`](docs/results.md) for full discussion including anti-correlations, non-critical inflation, and interpretation limits.

---

## Data sources

| Dataset | Source | Notes |
|---|---|---|
| ASTER L1T (v004) | NASA EarthData / LP DAAC | TIR bands B10–B14, 90 m resolution |
| MRDS national deposit database | USGS mrdata.usgs.gov | ~30,000 deposit records across 45 bboxes |
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

Binomial test (`scipy.stats.binomtest`, one-sided, greater) against a pooled null rate derived from the full 45-site survey. Two tests per site: critical-mineral-only (null = 8.5%) and all-deposit (null = 9.5%). MRDS deposits are clipped to the valid-pixel TIR footprint polygon before counting; zone coverage uses footprint area, not bbox area.

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

# All 45 sites (skip existing outputs)
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
│   ├── index.yaml                   # list of 45 site IDs
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
- `duckdb` — SQL-queryable national results
- `numpy` / `pandas` — array and tabular operations
- `matplotlib` — visualization

---

## Author

**Nicole Aikin** — MS Earth & Space Sciences, University of Washington (2025)
Metamorphic petrology · geochronology · ML pipelines for geoscience
[github.com/nicole-m-aikin](https://github.com/nicole-m-aikin)
