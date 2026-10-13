# Early-exploration prospectivity ML

## Claim (what this is / is not)

**Is:** A triage ranker for early exploration — “worth a closer look” — combining
thermal host-rock cues (ASTER AST_05 QI/CI/MI) with hyperspectral clay/alteration
cues (AVIRIS now; EMIT/EnMAP/PRISMA next), validated so the model cannot merely
memorize known mine districts or pre-1950 surface finds.

**Is not:** Drill targets, reserve estimates, a fused “ore detector,” or
Tetracorder-class mineral identification.

## Sensor hierarchy

| Role | Interface | Primary now | Next / global |
|------|-----------|-------------|----------------|
| Host filter | `ThermalScene` | ASTER AST_05 emissivity QI/CI/MI | AG100 mean; other TIR emissivity with band mapping |
| Alteration check | `HyperspectralScene` | AVIRIS-5 L2A + AVIRIS-Classic L2 (US GEMx census) | **EMIT** (~60 m, open LP DAAC), then EnMAP / PRISMA |

Every scene records `sensor_id`, GSD, wavelength sampling, QA mask, and coverage
fraction. Features are **normalized per scene** (percentile / equal-area ranks)
so ASTER 90 m and EMIT 60 m share a model after resampling to the study grid.
Sensor one-hots are diagnostics only — not headline predictors (they leak geography).

Airborne SWIR stays highest spectral fidelity but is campaign footprints, not
wall-to-wall global. MASTER remains partial-SWIR metadata, not a clay feature source.

## Feature contract (per 90 m cell)

- **TIR:** QI/CI/MI, combined score, within-scene ranks, distance to strong zone
- **HSI:** Al-OH / kaolinite / alunite depths + ranks, equal-area clay mask, `hsi_missing`
- **Context:** SGMC stratum code (bedrock/cover); structure distance when configured
- **Label:** binary near critical MRDS (default 180 m) + stratum-matched hard negatives

## Blocked validation (mandatory reporting rules)

1. **Primary — site-blocked CV:** `GroupKFold(groups=site_id)`. All pixels from a
   site stay together. Headline metrics are **site-means** (hit rate, enrichment,
   median rank gap, PR-AUC) — never pool pixels across sites as the claim.
2. **Secondary — era-blocked stress test:** train on `disc_yr < 1950` positives
   (+ backgrounds); evaluate on `disc_yr ≥ 1950`. Undated held out of both.
   Collapse on post-1950 ⇒ ascertainment leakage. Report **n dated** explicitly
   (~14% of critical MRDS have usable dates — underpowered).
3. **Do not average** site-CV and era-block into one number. Publish site-blocked
   as the primary claim; era-blocked as a robustness panel.
4. **Ablations:** TIR-only vs HSI-only vs TIR+HSI; sites with `hsi_missing=1`
   evaluate as TIR-only domain.

## Outputs

| Path | Role |
|------|------|
| `results/gemx_census_feature_ready.csv` | Frozen feature-ready status (TIR required; HSI optional) |
| `data/prospectivity/features/{site}.parquet` | Per-site feature/label store |
| `results/prospectivity_site_cv.csv` | Per-test-site metrics (primary) |
| `results/prospectivity_site_cv_summary.json` | Site-mean headlines + ablations |
| `results/prospectivity_era_block.json` | Era-block robustness panel |
| `results/emit_coverage.csv` | EMIT (+ EnMAP/PRISMA placeholder) inventory |
| `data/emit/{site}/` | EMIT adapter pilot Al-OH GeoTIFFs |
| `docs/prospectivity_ml.md` | This claim document |

## Census note (Classic size gate)

AVIRIS-Classic full cubes can exceed 100 GB per flight line. Downloads skip
cubes larger than 45 GB (`AVIRIS_CLASSIC_MAX_CUBE_GB`, default 45) and continue
with remaining greedy-selected lines. Seed selection keeps **at most one**
granule per must-contain point so town-overlapping strips are not all forced.

## Commands

```bash
# Freeze census → feature-ready table (safe to re-run as Classic finishes)
conda run -n aster-minerals python scripts/gemx_census_freeze.py

# Feature store (uses feature_ready=1 sites)
conda run -n aster-minerals python scripts/prospectivity_feature_store.py

# Train MVP (site-CV + era-block)
conda run -n aster-minerals python scripts/prospectivity_train.py --write-maps

# Global HSI inventory (EMIT CMR; EnMAP/PRISMA placeholders)
conda run -n aster-minerals python scripts/emit_coverage_inventory.py

# EMIT adapter contract (synthetic by default; --download for EarthData)
conda run -n aster-minerals python scripts/emit_adapter_pilot.py --synthetic
```

## Global gap (explicit)

Tasked satellites (EnMAP, PRISMA) and ISS EMIT’s arid acquisition mask are **not**
wall-to-wall global coverage. Expanding beyond the US GEMx census needs new
credentials/APIs and cloud screening; the feature contract stays fixed so new
adapters plug in without rewriting the model.
