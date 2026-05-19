# critical-minerals-aster

Spectral alteration mapping across 30 critical mineral sites in the US using ASTER thermal infrared (TIR) band ratio analysis. The pipeline identifies surface alteration zones, validates them against USGS MRDS mineral deposit data, and runs spatial significance tests to determine where — and for which deposit types — the method produces signal above chance.

**Results and figures:** [`docs/results.md`](docs/results.md) · [`figures/index.html`](figures/index.html)

---

## Scientific questions

1. Do ASTER-derived TIR alteration zones spatially correlate with known mineral occurrences in the USGS MRDS database?
2. Which TIR band ratio combinations best distinguish silica, carbonate, and mafic alteration?
3. Does the correlation vary by commodity type — and does that pattern make geological sense?
4. Are observed hit rates statistically significant above the null hypothesis of random deposit distribution?
5. Which Earth MRI deposit categories are TIR-detectable, and which are systematically invisible to this method?
6. Do non-critical bulk deposits (stone, sand, gravel) inflate headline hit rates — and is that signal geologically meaningful in its own right?

---

## Study sites (30)

| Site | State | Primary Deposit Type | Critical-only significance |
|---|---|---|---|
| Bisbee | AZ | Skarn / VMS Cu-Zn | p < 0.001 ** |
| McDermitt Caldera | NV/OR | Caldera (Li, Hg, Au) | p = 0.010 * |
| Mountain Pass | CA | Carbonatite REE | p = 0.008 * |
| Hanover-Fierro | NM | Cu-Mo-Zn skarn | p = 0.019 * |
| Yerington | NV | Porphyry Cu | p < 0.001 ** |
| Jerome | AZ | VMS Cu-Zn-Ag | p = 0.003 * |
| Thacker Pass | NV | Li clay | borderline (p = 0.080, n = 26) |
| Jerritt Canyon | NV | Carlin-type Au | — |
| Viburnum Trend | MO | MVT Pb-Zn | — |
| Darwin | CA | Polymetallic skarn | — |
| Lemhi Pass | ID/MT | REE phosphate veins | — |
| Carlin Trend | NV | Carlin-type Au | — |
| Ajo | AZ | Porphyry Cu | — |
| Gas Hills | WY | Roll-front U | — |
| Stillwater Complex | MT | PGM layered intrusion | — |
| Globe-Miami | AZ | Porphyry Cu | — |
| Tonopah–Manhattan | NV | Epithermal Au-Ag | — |
| Ely (Robinson) | NV | Porphyry Cu-Mo | — |
| Elk Creek | NE | Carbonatite (Nb, REE) | — |
| Green River Basin | WY | Trona / industrial minerals | false positive (aggregate only) |
| Silver Peak | NV | Li brine / epithermal | anti-correlated |
| Marysvale | UT | Uranium / epithermal | anti-correlated |
| Bear Lodge | WY | Carbonatite REE | anti-correlated |
| Bagdad | AZ | Porphyry Cu / skarn | anti-correlated |
| Oatman | NV/AZ | Low-sulfidation Au | anti-correlated |
| Steamboat Springs | NV | Geothermal / Au-Ag | anti-correlated |
| Pea Ridge | MO | IOCG (Fe-Cu-Au) | anti-correlated |
| Goldfield–Cuprite | NV | Epithermal Au | anti-correlated |
| Bingham Canyon | UT | Porphyry Cu-Mo-Au | anti-correlated |
| Climax | CO | Porphyry Mo | anti-correlated |

Significance = critical-mineral deposits only (Non-Critical excluded). See [`docs/results.md`](docs/results.md) for full discussion.

---

## Data sources

| Dataset | Source | Notes |
|---|---|---|
| ASTER L1T (v004) | NASA EarthData / LP DAAC | TIR bands B10–B14, 90 m resolution |
| MRDS national deposit database | USGS mrdata.usgs.gov | ~8,500 deposits across 30 bboxes |
| USGS Quaternary Faults | USGS QFAULTS REST API | Most sites |
| USGS SGMC fault data | USGS FeatureServer | Bear Lodge, Jerome, and others (all-age faults) |

**Note on SWIR availability:** ASTER SWIR bands (B04–B09), standard for clay/argillic mapping, are not available in LP DAAC v004 for these areas. TIR bands (B10–B14, 8–12 µm) are used instead, suited for silica, carbonate, and mafic mineral mapping in arid volcanic terranes.

---

## Methods

### Band ratios

| Ratio | Formula | Target mineralogy |
|---|---|---|
| Silica/quartz | B13/B14 | Silicic alteration, rhyolite |
| Carbonate/dolomite | B13/B12 | Hydrothermal carbonate |
| Mafic | B12/B13 | Mafic volcanic rocks |

### Classification

Percentile-based thresholds (70th/90th) applied per scene produce a 3-class anomaly map per ratio. An additive combined score (0–6) identifies pixels anomalous across multiple indicators. Strong anomaly zones (score ≥ 3) are vectorized to polygons via `rasterio.features.shapes`. Thresholds are scene-relative; cross-site comparison of raw scores is not meaningful.

### Statistical significance

Two complementary tests per site, run on critical-mineral deposits only:

- **Binomial test** — exact one-sided `scipy.stats.binomtest`; H₀: each deposit has probability p = zone area / bbox area of falling in a zone.
- **Spatial permutation test** — Monte Carlo (10,000 iterations); anomaly zone rasterised onto a 1,000×1,000 grid, random deposit placements counted each iteration.

Both tests are also run on all deposits combined; results are compared in [`docs/results.md`](docs/results.md).

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

Create a free account at [urs.earthdata.nasa.gov](https://urs.earthdata.nasa.gov). The pipeline uses `earthaccess.login(strategy="interactive")` on first run; credentials are cached.

### 3. Run the pipeline

```bash
# Single site
python -m critical_minerals_aster run --site mcdermitt

# Download from EarthData then process
python -m critical_minerals_aster run --site mcdermitt --download

# All 30 sites
python -m critical_minerals_aster run-batch --all-sites

# Regenerate national summary + synthesis figures
python -m critical_minerals_aster synthesize
```

### 4. Significance tests

```bash
# Critical-mineral-only (primary result)
python scripts/significance_critical_only.py

# Whole-catalog baseline
python scripts/compute_significance.py

# TIR-detectable mineral systems only
python scripts/compute_significance_filtered.py

# Per-(site × Earth MRI category)
python scripts/compute_significance_by_category.py
```

### 5. Query results

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
│   ├── index.yaml                   # list of 30 site IDs
│   └── {site_id}.yaml               # bbox, granule, classification params, structure layers
├── src/
│   └── critical_minerals_aster/
│       ├── config.py                # SiteConfig, ClassificationParams, StructureLayer
│       ├── paths.py                 # SitePaths — all file/dir paths
│       ├── spectral.py              # TIR I/O, granule selection, band ratios
│       ├── classification.py        # percentile classification, vectorization
│       ├── metrics.py               # MRDS spatial join, per-site summary CSV
│       ├── mrds.py                  # MRDS CSV → GeoDataFrame, Earth MRI / mineral-system classifiers
│       ├── structure.py             # distance-to-fault annotation, buffer flags
│       ├── significance.py          # binomial + spatial permutation p-values
│       ├── synthesis.py             # national summary CSV + figures (05, 06, 07)
│       ├── terrain.py               # hillshade DEM overlay
│       └── pipeline.py             # run_site() / run_batch() orchestration
├── docs/
│   ├── results.md                   # scientific findings, figures, interpretation
│   ├── roadmap.md                   # live planning doc
│   ├── architecture.md              # design rationale
│   └── structure_layers.md          # fault overlay config reference
├── scripts/
│   ├── significance_critical_only.py
│   ├── compute_significance.py
│   ├── compute_significance_filtered.py
│   ├── compute_significance_by_category.py
│   ├── fetch_new_site_structures.py
│   ├── download_usgs_faults.py
│   └── download_sgmc_structures.py
├── results/                         # generated per-site CSVs + results.duckdb
├── figures/                         # synthesis figures + per-site gallery
├── data/                            # not committed (ASTER rasters, MRDS CSV, structure GeoJSONs)
├── tests/
├── environment.yml
├── pyproject.toml
└── README.md
```

---

## Dependencies

- `rasterio` — raster I/O, feature extraction, rasterization for permutation test
- `geopandas` / `shapely` — vector operations and spatial joins
- `scipy` — binomial significance tests
- `earthaccess` — NASA EarthData authentication and download
- `duckdb` — SQL-queryable national results
- `numpy` / `pandas` — array and tabular operations
- `matplotlib` / `contextily` — visualization and basemap tiles

---

## Author

**Nicole Aikin** — MS Earth & Space Sciences, University of Washington (2025)
Metamorphic petrology · geochronology · ML pipelines for geoscience
[github.com/nicole-m-aikin](https://github.com/nicole-m-aikin)
