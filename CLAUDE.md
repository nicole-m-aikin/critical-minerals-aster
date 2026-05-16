# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project summary

Spectral alteration mapping of the McDermitt Caldera (NV/OR) using ASTER thermal infrared (TIR, B10–B14) band ratio analysis. Produces strong anomaly zone polygons and validates them against USGS MRDS mineral deposit data. Multi-site capable via per-site YAML configuration.

## Environment setup

```bash
conda env create -f environment.yml
conda activate aster-minerals
pip install -e .
```

The editable install (`pip install -e .`) is required for notebooks and tests to import `critical_minerals_aster` from `src/`.

## Commands

**Run pipeline (single site):**
```bash
python -m critical_minerals_aster run --site mcdermitt
python -m critical_minerals_aster run --site mcdermitt --download   # fetch from EarthData first
python -m critical_minerals_aster run --site mcdermitt --skip-figures
```

**Batch and synthesis:**
```bash
python -m critical_minerals_aster run-batch --all-sites
python -m critical_minerals_aster synthesize
```

**Tests:**
```bash
pytest tests/
pytest tests/test_classification.py   # single test file
pytest tests/test_classification.py::test_band_ratio_divide_by_zero  # single test
```

## Architecture

### Data flow

```
sites/{id}.yaml
    → SiteConfig (config.py)
    → SitePaths  (paths.py)       # resolves data/aster/, vectors/, figures/, results/
    → pipeline.run_site()         # orchestrates all steps
        ├── spectral.load_tir_bands_10_14()   → B10–B14 numpy arrays
        ├── spectral.alteration_ratios()       → silica, carbonate, mafic ratios
        ├── classification.classify_percentiles() → 3-class maps per ratio
        ├── classification.combined_score()    → additive score 0–6
        ├── classification.vectorize_strong_zones() → GeoDataFrame (GeoJSON)
        ├── metrics.compute_site_summary()    → MRDS spatial join + CSV
        └── pipeline.write_provenance()       → results/{id}_provenance.json
```

### Key modules (`src/critical_minerals_aster/`)

| Module | Purpose |
|---|---|
| `config.py` | `SiteConfig` dataclass + YAML loader; `ClassificationParams`, `StructureLayer` |
| `paths.py` | `SitePaths` — all file/dir paths, layout-aware (`flat` vs `nested`) |
| `spectral.py` | ASTER TIR I/O, granule selection scoring, `alteration_ratios()` |
| `classification.py` | Percentile classification, `combined_score()`, polygon vectorization |
| `metrics.py` | MRDS spatial join, per-site summary CSV |
| `mrds.py` | MRDS CSV → GeoDataFrame with CRS reprojection |
| `structure.py` | Distance-to-fault and buffer annotation for deposits |
| `synthesis.py` | Aggregate `results/*_summary.csv` → national summary |
| `pipeline.py` | `run_site()` / `run_batch()` orchestration + figure generation |
| `__main__.py` | `argparse` CLI entry point |

### Site configuration (`sites/`)

Each site has a YAML with: `id`, `name`, `bbox_wgs84`, `granule_id` (or `null` for auto-select), `layout` (`flat` for current McDermitt, `nested` for multi-site), `classification` params, `temporal` range, and optional `structure_layers`.

`sites/index.yaml` lists site ids for batch runs.

**Layout difference:** `flat` writes to `data/aster/` and `data/vectors/`; `nested` writes to `data/sites/{id}/aster/` and `data/sites/{id}/vectors/`. All path logic lives in `SitePaths`.

### Band ratios

- **Silica/quartz:** B13/B14
- **Carbonate/dolomite:** B13/B12
- **Mafic:** B12/B13

Classification uses per-scene percentile thresholds (default 70th/90th). Combined score ≥ 3 defines "strong anomaly" zones. Thresholds are scene-relative; cross-site comparison of raw scores is not meaningful.

### Outputs

- `data/vectors/strong_anomaly_zones.geojson` — vectorized strong anomaly polygons
- `figures/0*.png` — band ratio maps, classification maps, deposit overlay
- `results/{site_id}_summary.csv` — MRDS hit rates by commodity, Earth MRI category, and mineral system
- `results/{site_id}_provenance.json` — granule id, git commit, package versions

### Summary CSV row types

Each `results/*_summary.csv` contains four `row_type` values:

| `row_type` | Key column | Description |
|---|---|---|
| `site` | — | Site-level totals |
| `commodity` | `commodity_group` | Simplified commodity buckets (Gold/Silver, Uranium, etc.) |
| `earth_mri` | `earth_mri_category` | 8 Earth MRI priority groups (Energy, REE, Battery Metals…) |
| `mineral_system` | `mineral_system` | One of 24 OFR 2020-1042 mineral systems (Placer, Porphyry Cu-Mo-Au…) |

### Earth MRI / mineral-system classification (`mrds.py`)

Two independent classifiers operate on each MRDS row:

**`reclassify_mrds_earth_mri(df)`** — adds `earth_mri_category` from `commod1`.
Priority: Energy > REE > Battery Metals > PGM > Base Metals > Specialty/High-Tech > Gold/Silver > Industrial > Non-Critical.
Reflects the 2022 Final List of 50 Critical Minerals (includes scandium, cesium, rubidium, arsenic vs. the 2018 list).

**`reclassify_mrds_mineral_system(df)`** — adds `mineral_system` from a three-field cascade:
  1. `model` field (USGS deposit-model code string, e.g. `"53: Porphyry Cu"`) — most specific
  2. `dep_type` field (e.g. `"Placer"`, `"Stratabound"`)
  3. `commod1` field — broadest fallback

Coverage: ~84% of critical-mineral MRDS deposits are assigned a system; bulk non-critical (sand/gravel/stone) remain `"Unknown System"`.

**`is_critical_mineral(earth_mri_category)`** — returns `True` for all categories except `"Non-Critical"`.

### EarthData / ASTER data

Raw ASTER TIR rasters live under `data/aster/` (not committed). Files are named `{granule_id}_TIR_B{10-14}.tif`. The `--download` flag triggers `earthaccess.login(strategy="interactive")` on first run; subsequent runs use cached credentials. SWIR bands (B04–B09) are not available in LP DAAC v004 for this area.
