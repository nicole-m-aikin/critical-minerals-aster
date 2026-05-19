# critical-minerals-aster

Spectral alteration mapping across 30 critical mineral sites in the US using ASTER thermal infrared (TIR) band ratio analysis. The pipeline identifies surface alteration zones, validates them against USGS MRDS mineral deposit data, and runs spatial significance tests to determine where — and for which deposit types — the method produces signal above chance.

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

Significance = one-sided binomial + spatial permutation tests on **critical-mineral deposits only** (Non-Critical earth_mri_category excluded). Anti-correlated sites have p ≈ 1 — deposits actively avoid anomaly zones. Green River was significant in the all-deposits test (p = 0.005) but has only 1 critical deposit (0 hits); the headline rate was driven entirely by aggregate/industrial deposits.

---

## Data sources

| Dataset | Source | Notes |
|---|---|---|
| ASTER L1T (v004) | NASA EarthData / LP DAAC | TIR bands B10–B14, 90 m resolution |
| MRDS national deposit database | USGS mrdata.usgs.gov | ~8,500 deposits across 30 bboxes |
| USGS Quaternary Faults | USGS QFAULTS REST API | Most sites |
| USGS SGMC fault data | USGS FeatureServer | Bear Lodge, Jerome, and others (all-age faults) |

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

MRDS deposits within each scene bbox are spatially joined to strong anomaly zones. Hit rate = fraction of deposits falling inside a zone. Results are broken down by commodity group, Earth MRI category, and USGS mineral system. The primary significance metric uses **critical-mineral deposits only** (Non-Critical excluded).

### Statistical significance

Two complementary tests evaluate whether observed hit rates exceed chance:

**Binomial test** — exact one-sided test via `scipy.stats.binomtest`. Under H₀, each deposit has probability p = zone area / bbox area of falling in a zone. Tests whether observed hits significantly exceed that expectation.

**Spatial permutation test** — Monte Carlo (10,000 iterations). The anomaly zone union is rasterised onto a 1,000×1,000 grid; each iteration samples n_deposit random grid cells and counts zone hits. Returns P(random hits ≥ observed). Mathematically equivalent to placing random points uniformly in the bbox.

Both tests agree throughout (p-values never diverge by more than 0.02), confirming consistency. Tests are run twice: once on all MRDS deposits and once on critical-mineral deposits only (Non-Critical excluded).

### Structure proximity

MRDS deposits are annotated with distance to the nearest mapped fault. Deposits within 500 m of a fault are flagged as structurally controlled. Per-site structure GeoJSONs are fetched automatically from USGS QFAULTS (and SGMC as fallback) on first run.

---

## Key results

### Site-level

**6 of 30 sites** show critical-mineral hit rates significantly above chance (binomial + permutation, p < 0.05). 10 sites are anti-correlated — deposits actively avoid the anomaly zones. The anti-correlations are geologically coherent: Climax (deep porphyry Mo, no surface expression), Goldfield/Oatman (epithermal Au with eroded/covered alteration), Steamboat Springs (active geothermal surface ≠ MRDS deposit locations).

### Non-critical inflation

Non-critical deposits (stone, sand, gravel) have a **higher** pooled hit rate (9.8%) than critical minerals (7.0%). ASTER TIR detects silica and carbonate enrichment — exactly the mineralogy that makes rock commercially useful as aggregate. This inflates headline hit rates at sites with dense bulk-mineral MRDS coverage.

**Green River** is the clearest example: 52 of 53 deposits are non-critical, producing a 15.1% all-deposits hit rate (p = 0.005). Critical-only: 1 deposit, 0 hits, p = 1.0. The apparent significance was illusory. **Thacker Pass** is a genuine borderline case: 19.2% critical hit rate but only 26 critical deposits, leaving it underpowered (p = 0.080).

The stone/sand/gravel correlation is also a **real economic signal** for aggregate quarry siting — ASTER TIR can screen for silicified host rock. For that application, use the all-deposits and commodity-level rows in the CSVs rather than the critical-only results.

### By Earth MRI category (national pooled, 30 sites)

| Category | Deposits | Hits | Hit rate | TIR-detectable? |
|---|---|---|---|---|
| REE | ~120 | ~30 | ~27% avg | Yes — carbonatite + veins |
| Battery Metals – Co/Ni | ~200 | ~25 | ~13% avg | Yes — mafic/ultramafic halos |
| Base Metals | ~1,200 | ~100 | ~9% avg | Yes — porphyry/skarn/VMS halos |
| PGM | ~30 | ~4 | ~12% avg | Yes — mafic ratio |
| Gold/Silver | ~2,500 | ~160 | ~6% avg | Partial — caldera/VMS yes, Carlin/placer no |
| Energy | ~400 | ~15 | ~4% avg | No — roll-front U, no surface expression |
| Industrial | ~400 | ~20 | ~5% avg | No |
| Non-Critical (aggregate) | ~1,000 | ~100 | ~10% avg | Yes — but different use case |

No category reaches national significance individually; the signal is site-specific.

### Structure proximity

The 6 significant sites span a wide range of mean deposit–fault distances (1–20 km), so structural proximity alone does not predict TIR detectability. The signal is more strongly driven by deposit type and host-rock alteration style.

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
| `05_national_hit_rates.png` | Stacked bar by Earth MRI category, critical minerals only, sites ordered by critical hit rate |
| `06_structure_hit_rate.png` | Log-scale scatter: mean fault distance vs critical-mineral hit rate; red = significant sites |
| `07_hitrate_comparison.png` | Paired dot plot: all-deposits vs critical-only hit rate per site; exposes aggregate inflation |
| `index.html` | Sortable site gallery (no external deps) |

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
│   ├── compute_significance.py                 # whole-catalog binomial + permutation, all sites
│   ├── compute_significance_filtered.py        # TIR-detectable systems only
│   ├── compute_significance_by_category.py     # per-(site × Earth MRI category)
│   ├── significance_critical_only.py           # critical-mineral-only retest, excludes Non-Critical
│   ├── fetch_new_site_structures.py            # batch USGS fault fetch for new sites
│   ├── download_usgs_faults.py
│   └── download_sgmc_structures.py
├── results/                         # generated per-site CSVs + results.duckdb
│   ├── {site_id}_summary.csv
│   ├── {site_id}_provenance.json
│   ├── significance_critical_only.csv          # critical-only p-values, all 30 sites
│   ├── national_summary.csv
│   └── results.duckdb
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

# All 30 sites in batch
python -m critical_minerals_aster run-batch --all-sites

# Regenerate national summary + synthesis figures (05, 06, 07)
python -m critical_minerals_aster synthesize
```

### 4. Significance tests

```bash
# Whole-catalog: binomial + permutation for all 30 sites
python scripts/compute_significance.py

# Critical-mineral-only retest (excludes stone/sand/gravel from denominator)
python scripts/significance_critical_only.py

# TIR-detectable mineral systems only
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
- **Non-critical inflation:** Sites with dense aggregate/stone MRDS coverage (Green River, Silver Peak, Marysvale) show inflated all-deposits hit rates. Use `significance_critical_only.csv` for critical-mineral conclusions.
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
