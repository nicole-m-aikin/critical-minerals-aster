# ASTER TIR: local rasterio vs Google Earth Engine parity

A reproducibility / comparison module, not a rewrite of the 109-site
pipeline.  The question is implementation-level:

> Does a cloud-native Earth Engine recipe reproduce the local ASTER TIR
> workflow, and which preprocessing / statistical choices move the
> classified map?

That is the same class of question this project has already answered
once: a diagnostic script that skipped the bbox clip computed 70th/90th
percentiles over a several-fold-too-large raster and reversed three
"continuous-only" sites.  A GEE port has more places to do the same
thing silently.

## What is compared

For the 8-site panel in `DEFAULT_PARITY_SITES`
(`src/critical_minerals_aster/gee_parity.py`):

| Site | Role in the paper |
|---|---|
| `goldfield` | HS epithermal; pinned granule; SWIR pilot |
| `bisbee` | Robust arid skarn |
| `christmas` | Robust skarn replicate |
| `tombstone` | Robust CRD analog |
| `mcdermitt` | Original / first site |
| `eureka` | FDR then clustering failure |
| `mountain_pass` | Reversed by the site-specific null |
| `marysvale` | HS + alunite; SWIR complementary |

Both sides perform the same six steps:

1. Footprint / scene selection (coverage score, 30% gate, optional pin)
2. TIR B10–B14 retrieval
3. Silica (B13/B14), carbonate (B13/B12), mafic (B12/B13)
4. Clip to `site.bbox_wgs84` **before** percentiles
5. 70th / 90th percentile thresholds
6. Combined score 0–6; strong = score ≥ 3

Agreement is pixel-wise on the overlap: overall accuracy, Cohen's κ on
the 0–6 score, Jaccard of the strong-zone mask, MAE of the score, and
Δ of the 70th/90th silica thresholds.

## The GEE port

The Earth Engine recipe is generated as Code Editor JavaScript
(`results/gee_parity_recipes/{site}.js`) and as a Python `ee` graph
(`build_gee_image`).  Collection: `ASTER/AST_L1T_003`.  Local granules
are LP DAAC **v004** via earthaccess.  Same instrument, different
archive cut — scene ids will not always match.

Live execution (`--gee`) needs `earthengine-api` and a Cloud project:

```bash
pip install 'critical-minerals-aster[gee]'
earthengine authenticate
python -m critical_minerals_aster gee-parity --gee --project YOUR_GCP_PROJECT
```

CI and the default benchmark never import `ee`.

## Offline discrepancy probes

A live GEE pull is not required to learn *why* the two stacks disagree.
Each probe changes one implementation choice and re-classifies the same
local DNs:

| Probe | Local contract | Cloud-native analogue |
|---|---|---|
| `clip_order` | `clip_bands_to_bbox` then `nanpercentile` | `reduceRegion` on the unclipped scene footprint |
| `histogram_percentile` | `np.nanpercentile` (exact) | `ee.Reducer.percentile` (histogram bins) |
| `radiance_ucc` | raw DN ratios (0 = nodata) | radiance = (DN − 1) × UCC; UCC differs by TIR band |
| `reproject_4326` | native UTM ~90 m grid | GEE `reduceRegion` default CRS is EPSG:4326 |

`radiance_ucc` only runs when raw B12–B14 GeoTIFFs exist (not on
ratio-only mosaics).

Offline panel result (8 sites, local rasters):

| Probe | n | Mean strong-zone Jaccard | Mean \|Δ strong fraction\| |
|---|---|---|---|
| clip order | 8 | 0.676 | 0.024 |
| EPSG:4326 resample | 8 | 0.693 | 0.003 |
| histogram percentile | 8 | 1.000 | 0.000 |
| radiance UCC | 5 | 0.998 | 0.000 |

Two different kinds of disagreement. **Clip order** moves both which
pixels are strong *and* the area fraction that becomes `p_null`.
**EPSG:4326** (GEE's default `reduceRegion` CRS) remaps the zones
(Jaccard ~0.69) but the strong fraction stays ~10% because percentiles
are recomputed on the new grid. Histogram binning and the TIR UCC
offset are not material at these DN values: a per-band constant scale
preserves rank, so percentile classes do not move.

Goldfield is the loudest clip-order site (Jaccard 0.372, Δstrong −0.073).
That is the same bug class this project already found in its own
diagnostic scripts.

## Commands

```bash
# Offline panel (local rasters required)
python -m critical_minerals_aster gee-parity
python scripts/fig_gee_parity.py

# Subset
python -m critical_minerals_aster gee-parity --sites goldfield bisbee

# Live Earth Engine
python -m critical_minerals_aster gee-parity --gee --project YOUR_GCP_PROJECT
```

Outputs: `results/gee_parity_agreement.csv`,
`results/gee_parity_agreement.json`,
`results/gee_parity_recipes/*.js`,
`figures/gee_parity_agreement.png`.

## What this is not

- Not a 109-site GEE rewrite.
- Not a claim that GEE and LP DAAC v004 are the same product.
- Not mineral mapping.  Same TIR ratios as the rest of the repo.
- Live GEE numbers are absent until `--gee` is run.  The figure's
  default story is the probes: which implementation choice moves the map.
