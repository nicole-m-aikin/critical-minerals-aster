# Structural geology layers (optional)

Add line or polygon structure data to a site YAML under `structure_layers`:

```yaml
structure_layers:
  - path: data/structures/mcdermitt_faults.geojson
    type: faults
    label: "Quaternary Faults (QFaults)"   # optional; defaults to type.title()
    buffer_m: 500
```

Supported `type` values: `faults`, `contacts`, `folds`.

When layers are configured, the CLI records in `results/{site_id}_provenance.json`:

- `n_deposits_on_structure` — MRDS points within the buffer of any structure
- `mean_nearest_structure_m` — mean distance from deposits to structure linework

Use `critical_minerals_aster.structure.annotate_deposits_with_structure()` in notebooks for per-deposit `nearest_structure_m` and `on_structure` columns.

**Data sources in use:**
- USGS Quaternary Faults REST API — 13 of 15 sites
- USGS SGMC FeatureServer (all geological ages) — Bear Lodge (WY) and Jerome (AZ), where Quaternary faults are absent

Download scripts: `scripts/download_usgs_faults.py` (Quaternary), `scripts/download_sgmc_structures.py` (SGMC + WSGS fallback).
