#!/usr/bin/env python
"""
Consistency check of the hand-curated deposit-system categories.

results/deposit_system_categories.csv is typed by hand (one analyst, not blind
to outcome). CIM (2018) §2.6.6 / §2.7.4 asks for a data-verification step on
manually entered data. A true double-entry check needs a second person; this
script does the automatable half and produces the review list for that person:

1. Structural validation — every site in sites/index.yaml has exactly one row,
   the category is from the allowed set, confidence is high/medium/low, and the
   basis text is non-empty.
2. Cross-check against the independent MRDS signal — mineral systems of the
   validated MRDS records inside each site's valid TIR footprint, assigned
   from the ``model`` and ``dep_type`` fields ONLY. The commodity fallback in
   ``reclassify_mrds_mineral_system`` is deliberately not used here: it labels
   any gold record "Orogenic" and any Pb/Zn record "Basin Brine Path", which
   would make almost every district disagree. A site is flagged ``review``
   when none of its top-2 geology-derived systems is compatible with the hand
   category, and ``no_check`` when fewer than MIN_TYPED records carry model or
   dep_type information. A flag is a prompt for a second geologist, not an
   error (docs/learning_log.md 1.5).

Outputs
-------
    results/deposit_category_review.csv

Usage
-----
    conda run -n aster-minerals python scripts/check_deposit_categories.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from critical_minerals_aster.config import load_site_by_id
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.mrds import _classify_mineral_system, filter_mrds_bbox, mrds_to_points_gdf
from critical_minerals_aster.paths import site_paths_for

COMPATIBLE: dict[str, set[str]] = {
    "Skarn/Carbonate Replacement": {"Skarn / Carbonate Replacement", "Coeur d'Alene-type"},
    "Porphyry Cu-Mo-Au": {"Porphyry Cu-Mo-Au", "Climax-type", "Skarn / Carbonate Replacement"},
    "Sediment-hosted": {"Carlin-type", "Basin Brine Path", "Marine Chemocline"},
    # Epithermal-vein model strings are captured by the Porphyry rule in mrds.py.
    "Epithermal": {"Meteoric Convection", "Porphyry Cu-Mo-Au"},
    "Alkaline/Carbonatite": {"Magmatic REE", "Hybrid Magmatic REE"},
    "Volcanogenic/VMS": {"Volcanogenic Seafloor"},
    "Uranium/Energy": {"Meteoric Recharge", "Petroleum"},
    "Mafic/Ultramafic": {"Mafic Magmatic"},
    "Mixed/Other": set(),
}
CONFIDENCE = {"high", "medium", "low"}
MIN_TYPED = 3


def geology_typed_systems(site_id: str) -> pd.Series:
    """Mineral system per validated MRDS record in the site footprint, from
    model/dep_type only (records with neither field come back Unknown)."""
    site = load_site_by_id(site_id, REPO_ROOT / "sites")
    paths = site_paths_for(site, REPO_ROOT)
    prov_path = REPO_ROOT / "results" / f"{site_id}_provenance.json"
    prov = json.loads(prov_path.read_text()) if prov_path.exists() else {}
    bbox = tuple(prov.get("raster_bbox_wgs84") or site.bbox_wgs84)
    local = filter_mrds_bbox(read_mrds_national(paths), bbox)
    if prov.get("tir_footprint_wgs84") is not None and len(local):
        fp = gpd.GeoDataFrame.from_features(prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326")
        local = gpd.clip(mrds_to_points_gdf(local, "EPSG:4326"), fp)
    return pd.Series(
        [_classify_mineral_system(m, d, "") for m, d in zip(local.get("model", []), local.get("dep_type", []))],
        dtype=object,
    )


def main() -> int:
    site_ids = yaml.safe_load((REPO_ROOT / "sites" / "index.yaml").read_text())["sites"]
    cats = pd.read_csv(REPO_ROOT / "results" / "deposit_system_categories.csv")

    problems: list[str] = []
    missing = sorted(set(site_ids) - set(cats["site_id"]))
    extra = sorted(set(cats["site_id"]) - set(site_ids))
    dups = sorted(cats.loc[cats["site_id"].duplicated(), "site_id"])
    if missing:
        problems.append(f"sites without a category: {missing}")
    if extra:
        problems.append(f"categories for unknown sites: {extra}")
    if dups:
        problems.append(f"duplicate category rows: {dups}")
    bad_cat = cats.loc[~cats["category"].isin(COMPATIBLE), "site_id"].tolist()
    bad_conf = cats.loc[~cats["confidence"].str.lower().isin(CONFIDENCE), "site_id"].tolist()
    no_basis = cats.loc[cats["basis"].fillna("").str.strip() == "", "site_id"].tolist()
    if bad_cat:
        problems.append(f"category not in allowed set: {bad_cat}")
    if bad_conf:
        problems.append(f"confidence not high/medium/low: {bad_conf}")
    if no_basis:
        problems.append(f"empty basis: {no_basis}")

    rows = []
    for _, r in cats.iterrows():
        systems = geology_typed_systems(r.site_id)
        counts = systems[systems != "Unknown System"].value_counts()
        top = list(zip(counts.index, counts.astype(int)))[:2]
        n_known = int(counts.sum())
        compatible = COMPATIBLE.get(r.category, set())
        if r.category == "Mixed/Other" or n_known < MIN_TYPED:
            status = "no_check"
        elif any(sys_name in compatible for sys_name, _ in top):
            status = "consistent"
        else:
            status = "review"
        rows.append({
            "site_id": r.site_id,
            "category": r.category,
            "confidence": r.confidence,
            "mrds_top1": top[0][0] if top else "",
            "mrds_top1_n": top[0][1] if top else 0,
            "mrds_top2": top[1][0] if len(top) > 1 else "",
            "mrds_top2_n": top[1][1] if len(top) > 1 else 0,
            "mrds_n_geology_typed": n_known,
            "status": status,
        })

    out = pd.DataFrame(rows).sort_values(["status", "site_id"])
    out_path = REPO_ROOT / "results" / "deposit_category_review.csv"
    out.to_csv(out_path, index=False)

    print("Structural validation:", "OK" if not problems else "")
    for p in problems:
        print("  PROBLEM:", p)
    print(f"\nCross-check vs MRDS mineral_system -> {out_path.relative_to(REPO_ROOT)}")
    print(out["status"].value_counts().to_string())
    review = out[out["status"] == "review"]
    if not review.empty:
        print("\nFor second-geologist review:")
        print(review[["site_id", "category", "confidence", "mrds_top1", "mrds_top2"]].to_string(index=False))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
