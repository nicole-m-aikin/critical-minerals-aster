# Structural geology layers (optional)

Add line or polygon structure data to a site YAML under `structure_layers`:

```yaml
structure_layers:
  - path: data/structures/mcdermitt_faults.geojson
    type: faults
    buffer_m: 500
```

Supported `type` values: `faults`, `contacts`, `folds`.

When layers are configured, the CLI records in `results/{site_id}_provenance.json`:

- `n_deposits_on_structure` — MRDS points within the buffer of any structure
- `mean_nearest_structure_m` — mean distance from deposits to structure linework

Use `critical_minerals_aster.structure.annotate_deposits_with_structure()` in notebooks for per-deposit `nearest_structure_m` and `on_structure` columns.

**Data sources (future):** USGS Quaternary faults, state geologic map WFS, SGMC contacts.
