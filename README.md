# critical-minerals-aster

Spectral alteration mapping across 15 critical mineral sites in the western US using ASTER thermal infrared (TIR) band ratio analysis. The pipeline identifies surface alteration zones, validates them against USGS MRDS mineral deposit data, and runs spatial significance tests to determine where — and for which deposit types — the method produces signal above chance.

---

## Scientific questions

1. Do ASTER-derived TIR alteration zones spatially correlate with known mineral occurrences in the USGS MRDS database?
2. Which TIR band ratio combinations best distinguish silica, carbonate, and mafic alteration?
3. Does the correlation vary by commodity type — and does that pattern make geological sense?
4. Are observed hit rates statistically significant above the null hypothesis of random deposit distribution?
5. Which Earth MRI deposit categories are TIR-detectable, and which are systematically invisible to this method?

---

## Study sites (15)

| Site | State | Primary Deposit Type | Significance |
|---|---|---|---|
| Yerington | NV | Porphyry Cu | p = 0.000048 ** |
| Mountain Pass | CA | Carbonatite REE | p = 0.001 ** |
| Jerome | AZ | VMS Cu-Zn-Ag | p = 0.010 * |
| McDermitt | NV/OR | Caldera (Li, Hg, Au) | p = 0.013 * |
| Silver Peak | NV | Li brine / epithermal | — |
| Jerritt Canyon | NV | Carlin-type Au | — |
| Tonopah | NV | Epithermal Au-Ag | — |
| Darwin | CA | Polymetallic skarn | — |
| Stillwater | MT | PGM layered intrusion | — |
| Marysvale | UT | Uranium / epithermal | — |
| Bear Lodge | WY | Carbonatite REE | — |
| Oatman | AZ | Low-sulfidation Au | anti-correlated |
| Goldfield | NV | Epithermal Au | anti-correlated |
| Climax | CO | Porphyry Mo | anti-correlated |
| Steamboat Springs | NV | Geothermal / Au-Ag | anti-correlated |

Significance = one-sided binomial test, H₀: hit rate ≤ zone coverage fraction. Anti-correlated sites have deposits actively avoiding anomaly zones (p ≈ 1).

---

## Data sources

| Dataset | Source | Notes |
|---|---|---|
| ASTER L1T (v004) | NASA EarthData / LP DAAC | TIR bands B10–B14, 90 m resolution |
| MRDS national deposit database | USGS mrdata.usgs.gov | ~5,000 deposits across 15 bboxes |
| USGS Quaternary Faults | USGS QFAULTS REST API | 13 sites |
| USGS SGMC fault data | USGS FeatureServer | Bear Lodge, Jerome (all-age faults) |

**Note on SWIR availability:** ASTER SWIR bands (B04–B09), standard for clay/argillic mapping, are not available in LP DAAC v004 for these areas. TIR bands (B10–B14, 8–12 µm) are used instead, which are well-suited for silica, carbonate, and mafic mineral mapping in arid volcanic terranes.

---

## Methods

### Band ratios

| Ratio | Formula | Target mineralogy |
|---|---|---|
| Silica/quartz | B13/B14 | Silicic alteration, rhyolite |
| Carbonate/dolomite | B13/B12 | Hydrothermal carbonate |
| Mafic | B12/B13 | Mafic volcanic rocks |

### Classification

Percentile-based thresholds (70th/90th) applied per scene produce a 3-class anomaly map per ratio. An additive combined score (0–6) identifies pixels anomalous across multiple indicators. Strong anomaly zones (score ≥ 3) are vectorized to polygons via `rasterio.features.shapes`.

Classification thresholds are scene-relative; cross-site comparison of raw scores is not meaningful — use hit rates and zone coverage fractions instead.

### Deposit validation

MRDS deposits within each scene bbox are spatially joined to strong anomaly zones. Hit rate = fraction of deposits falling inside a zone. Results are broken down by commodity group, Earth MRI category, and USGS mineral system.

### Statistical significance

Two complementary tests evaluate whether observed hit rates exceed chance:

**Binomial test** — exact one-sided test via `scipy.stats.binomtest`. Under H₀, each deposit has probability p = zone area / bbox area of falling in a zone. Tests whether observed hits significantly exceed that expectation.

**Spatial permutation test** — Monte Carlo (10,000 iterations). The anomaly zone union is rasterised onto a 1,000×1,000 grid; each iteration samples n_deposit random grid cells and counts zone hits. Returns P(random hits ≥ observed). Mathematically equivalent to placing random points uniformly in the bbox.

Both tests agree throughout (p-values never diverge by more than 0.01), confirming consistency.

### Structure proximity

MRDS deposits are annotated with distance to the nearest mapped fault. Deposits within 500 m of a fault are flagged as structurally controlled. Per-site structure GeoJSONs are fetched automatically from USGS QFAULTS (and SGMC as fallback) on first run.

---

## Key results

### Site-level

4 of 15 sites show hit rates significantly above chance. 5 sites are anti-correlated — deposits actively avoid the anomaly zones. The anti-correlations are geologically coherent: Climax (deep porphyry Mo, no surface expression), Goldfield/Oatman (epithermal Au with eroded/covered alteration), Steamboat Springs (active geothermal surface ≠ MRDS deposit locations).

### By Earth MRI category (national pooled)

| Category | Deposits | Hits | Expected | Ratio | TIR-detectable? |
|---|---|---|---|---|---|
| Base Metals | 951 | 92 | 85 | 1.08 | Yes — porphyry/skarn/VMS halos |
| Battery Metals – Co/Ni | 202 | 20 | 18 | 1.12 | Maybe |
| PGM | 25 | 3 | 2.4 | 1.26 | Yes — mafic ratio |
| REE | 28 | 3 | 3.5 | 0.87 | Yes (carbonatite), but small n |
| Gold/Silver | 2145 | 147 | 212 | 0.69 | No — placer/Carlin/epithermal |
| Energy | 289 | 11 | 30 | 0.37 | No — roll-front U |
| Battery Metals – Li/Brine | 9 | 0 | 0.9 | 0.00 | No — brine-hosted |
| Industrial | 214 | 11 | 24 | 0.47 | No — sand/gravel/stone |

No category reaches national significance individually — the signal is site-specific rather than category-wide. The strongest per-site results are Base Metals at Mountain Pass (ratio 2.53, p < 0.001) and Yerington (ratio 2.50, p = 0.001).

### Dilution finding

The whole MRDS catalog dilutes the signal, but the main diluter is **within-category heterogeneity** rather than cross-category mixing. "Gold/Silver" spans placer/Carlin (TIR-invisible) and caldera-hosted/VMS-associated Au (TIR-visible): Jerome and McDermitt both show significant Gold/Silver hits because those deposits are co-spatial with the alteration system. Filtering by mineral system helps but reduces n too much for most individual sites.

---

## Figures

Each site generates five figures:

| Figure | Content |
|---|---|
| `00_composite_rgb.png` | False-color TIR composite |
| `01_tir_band_ratios.png` | Three band ratio maps with global colorbars |
| `02_classification.png` | Per-ratio classification + combined score |
| `03_deposit_overlay.png` | Anomaly zones, MRDS deposits, fault corridors, scale bar |
| `05_structure_proximity.png` | Strip chart: deposit distance to nearest fault by commodity group |

National synthesis figures in `figures/`:

| Figure | Content |
|---|---|
| `05_national_hit_rates.png` | Stacked bar by Earth MRI category across all sites |
| `06_structure_hit_rate.png` | Log-scale scatter: mean fault distance vs hit rate |
| `index.html` | Sortable site gallery (no external deps) |

---

## Repo structure

```
critical-minerals-aster/
├── sites/
│   ├── index.yaml                   # list of 15 site IDs
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
│       ├── synthesis.py             # national summary CSV + figures
│       ├── terrain.py               # hillshade DEM overlay
│       └── pipeline.py             # run_site() / run_batch() orchestration
├── docs/
│   ├── roadmap.md                   # live planning doc — current status, backlog
│   ├── architecture.md              # design rationale
│   └── structure_layers.md          # fault overlay config reference
├── notebooks/
│   ├── 00_verify_setup.ipynb
│   ├── 01_data_download.ipynb
│   ├── 02_band_ratios.ipynb
│   ├── 03_classification.ipynb
│   ├── 04_deposit_overlay.ipynb
│   └── 05_national_synthesis.ipynb
├── scripts/
│   ├── synthesize_national.py
│   ├── compute_significance.py            # whole-catalog binomial + permutation, all sites
│   ├── compute_significance_filtered.py   # TIR-detectable systems only
│   ├── compute_significance_by_category.py # per-(site × Earth MRI category)
│   ├── download_usgs_faults.py
│   └── download_sgmc_structures.py
├── results/                         # generated per-site CSVs + results.duckdb
├── tests/
├── data/                            # not committed (ASTER rasters, MRDS CSV, structure GeoJSONs)
├── figures/
├── environment.yml
├── pyproject.toml
└── README.md
```

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
# Single site (uses cached ASTER rasters)
python -m critical_minerals_aster run --site mcdermitt

# Download from EarthData then process
python -m critical_minerals_aster run --site mcdermitt --download

# All 15 sites in parallel (4 workers, ~90 s)
python -m critical_minerals_aster run-batch --all --workers 4

# Skip already-processed sites
python -m critical_minerals_aster run-batch --all --workers 4 --skip-existing

# Regenerate national summary + synthesis figures
python scripts/synthesize_national.py
```

### 4. Significance tests

```bash
# Whole-catalog: binomial + permutation for all 15 sites
python scripts/compute_significance.py

# TIR-detectable mineral systems only (filters MRDS before testing)
python scripts/compute_significance_filtered.py

# Per-(site × Earth MRI category) binomial test + national pooled test
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

## Interpretation limits

- **TIR-only:** SWIR clay/argillic mapping (B04–B09) is not available in LP DAAC v004 for these areas.
- **Scene-relative thresholds:** Percentile classification is per-scene; don't compare raw scores across sites.
- **MRDS uncertainty:** Deposit locations are report-derived and may be offset from true geology.
- **Significance null model:** Both tests assume deposits are uniformly distributed within the bbox. Spatial clustering of real deposits means p-values are conservative — the true null distribution is not uniform.
- **Anti-correlations are informative:** p ≈ 1 at a site means the method is physically incapable of detecting the dominant deposit type there, not that the zones are wrong.

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
