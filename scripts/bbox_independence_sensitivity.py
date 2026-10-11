#!/usr/bin/env python
"""
Box-independence sensitivity for the FDR-significant sites (publication audit:
"MRDS ascertainment and target leakage").

Each site's classification thresholds (70th/90th percentile of each band
ratio) and its deposit sample are both defined by the analyst-drawn
``bbox_wgs84``. If boxes were drawn or nudged toward MRDS clusters, the
validation domain is outcome-informed. This script:

1. Documents how each box was derived (results/bbox_provenance.csv), from the
   site YAML comments, scripts/expansion_add_sites.py, the expansion docs and
   `git log -p -- sites/{id}.yaml` (number of bbox revisions).

2. Re-runs the exact classification recipe of ``pipeline.run_classification``
   (ratio-mosaic branch) and the critical-mineral hit test of
   ``scripts/site_specific_null_significance.py`` on extents the analyst did
   not draw. The cached ``{id}_mosaic_ratio_*.tif`` files are written at the
   full granule-union extent (only the ``_cls_*`` files are bbox-clipped), so
   the whole-scene test runs offline. Variants:

     analyst_box              reproduction of the published box (sanity check)
     whole_scene              percentiles + deposits over the full valid footprint
     whole_scene_outside_box  whole-scene zones, scored only outside the analyst
                              raster window (territory the analyst never framed)
     scene_centred_box        0.46 deg template box centred on the scene-footprint
                              centroid (depends on granule geometry only)
     shift_*/shrink/expand    8 perturbations of the analyst box (+/-15%),
                              fully reclassified

   If a site's cached ratio mosaic does not extend past its box (whole-scene
   valid pixels < 1.1x box valid pixels), the whole-scene test is impossible
   offline; the 8 perturbations are then scored on the cached zones and
   footprint with NO reclassification (method column = "cached_no_reclass").

Outputs
-------
    results/bbox_provenance.csv
    results/bbox_independence_sensitivity.csv
    figures/verification_bbox_independence.png
    data/sites/{id}/rasters/bbox_sensitivity/whole_scene_zones.gpkg  (gitignored)

Usage:
    conda run -n aster-minerals python scripts/bbox_independence_sensitivity.py
    conda run -n aster-minerals python scripts/bbox_independence_sensitivity.py --sites bisbee tombstone
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from shapely.geometry import box

from critical_minerals_aster.bbox_sensitivity import (
    TEMPLATE_BOX_DEG,
    bbox_overlap_fraction,
    centred_box,
    classify_strong_zones,
    critical_hit_test,
    footprint_centroid_wgs84,
    perturbation_variants,
    raster_extent_polygon,
    restrict_to_region,
)
from critical_minerals_aster.config import ClassificationParams, load_site_by_id
from critical_minerals_aster.metrics import read_mrds_national
from critical_minerals_aster.paths import site_paths_for
from critical_minerals_aster.pipeline import resolve_granule_id
from critical_minerals_aster.significance import _bbox_in_crs
from critical_minerals_aster.spectral import (
    alteration_ratios,
    clip_bands_to_bbox,
    load_ratio_mosaic,
    load_tir_bands_10_14,
    raster_bbox_wgs84,
)

REPO_ROOT = Path(__file__).parent.parent
RESULTS = REPO_ROOT / "results"
FIG_PATH = REPO_ROOT / "figures" / "verification_bbox_independence.png"

FDR_SIG_20 = [
    "bisbee", "christmas", "contact_district", "courtland_gleeson", "crooks_gap",
    "eagle_mountain", "eureka", "holden", "karnes_county", "lake_valley",
    "lisbon_valley", "lordsburg", "magdalena_kelly", "mascot_jefferson_city",
    "organ_mountains", "providence_mountains", "shafter", "tombstone", "tri_state",
    "uravan",
]
ROBUST_13 = {
    "bisbee", "christmas", "courtland_gleeson", "crooks_gap", "eagle_mountain",
    "holden", "karnes_county", "lake_valley", "magdalena_kelly", "organ_mountains",
    "providence_mountains", "shafter", "tombstone",
}
PERTURB_FRAC = 0.15
FULL_SCENE_MIN_RATIO = 1.1

# From docs/expansion_candidates.md: "bbox?" centres flagged for an MRDS nudge,
# and the separate MRDS spot-check list. The build log there records that the
# flagged boxes were "spot-checked against USGS/mindat -- all within tolerance,
# kept as-is".
MRDS_NUDGE_FLAGGED = {"alligator_ridge", "cherry_creek", "contact_district", "lemitar", "arizona_strip"}
MRDS_SPOTCHECK_LISTED = {"magdalena_kelly", "lake_valley", "victorio", "providence_mountains",
                         "eagle_mountain", "josephine_ophiolite"}

# Hand-curated evidence for boxes that pre-date the scripted +58 expansion.
PRE_EXPANSION_EVIDENCE = {
    "bisbee": (
        "Added in b8f5e3a (30-site expansion). YAML comment gives geology only; no derivation "
        "recorded in docs, commit message or learning log. 0.50x0.45 deg box, never revised in git."
    ),
    "eureka": (
        "Added in ad4f470 (45-site batch) whose commit message says deposit TYPES were "
        "'prioritized based on highest observed hit rates' (outcome-informed at the site-selection "
        "level). No box derivation recorded. 0.70x0.55 deg box, never revised in git."
    ),
    "lordsburg": (
        "Added in ad4f470 (same outcome-prioritized batch as eureka). YAML comment locates the "
        "district '~30 km south of Lordsburg'; no box derivation recorded. 0.70x0.80 deg box, "
        "never revised in git."
    ),
    "christmas": (
        "Not in expansion_add_sites.py; first committed in 1b708a0 although added earlier (51-site "
        "stage), so any pre-commit edits are invisible to git. YAML comment: geology only. "
        "0.50x0.44 deg box."
    ),
    "courtland_gleeson": (
        "Not in expansion_add_sites.py; first committed in 1b708a0 (added at the 51-site stage; "
        "pre-commit edits invisible to git). YAML comment: 'adjacent district to Bisbee/Tombstone'. "
        "0.50x0.44 deg box overlapping the tombstone box."
    ),
    "tombstone": (
        "Not in expansion_add_sites.py; first committed in 1b708a0 (added at the 51-site stage, "
        "after Bisbee's result was known, explicitly 'to test whether Bisbee's fully robust signal "
        "replicates'). Pre-commit edits invisible to git. 0.50x0.44 deg box overlapping the bisbee "
        "and courtland_gleeson boxes."
    ),
}


# ---------------------------------------------------------------------------
# Task 1 -- box provenance
# ---------------------------------------------------------------------------


def _load_expansion_sites() -> dict[str, list[float]]:
    spec = importlib.util.spec_from_file_location(
        "expansion_add_sites", REPO_ROOT / "scripts" / "expansion_add_sites.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # module-level code only defines constants; main() is guarded
    return {row[0]: list(row[2]) for row in mod.SITES}


def _git_bbox_history(site_id: str) -> tuple[str, int]:
    """(first commit touching the YAML, number of commits that changed the bbox line)."""
    try:
        log = subprocess.check_output(
            ["git", "log", "--format=@@%h", "-p", "--", f"sites/{site_id}.yaml"],
            cwd=REPO_ROOT, stderr=subprocess.DEVNULL, text=True,
        )
    except Exception:
        return "", -1
    commits, n_changes, current, changed = [], 0, None, False
    for line in log.splitlines():
        if line.startswith("@@") and not line.startswith("@@ "):
            if current is not None and changed:
                n_changes += 1
            current, changed = line[2:], False
            commits.append(current)
        elif line.startswith("+bbox_wgs84"):
            changed = True
    if current is not None and changed:
        n_changes += 1
    # The creating commit adds the bbox line; only later changes are revisions.
    return (commits[-1] if commits else ""), max(n_changes - 1, 0)


def build_provenance(site_ids: list[str]) -> pd.DataFrame:
    expansion = _load_expansion_sites()
    boxes = {s: tuple(load_site_by_id(s, REPO_ROOT / "sites").bbox_wgs84) for s in site_ids}
    rows = []
    for s in site_ids:
        bbox = boxes[s]
        first_commit, n_rev = _git_bbox_history(s)
        in_exp = s in expansion
        matches = in_exp and np.allclose(expansion[s], bbox)
        overlaps = [
            f"{o}:{bbox_overlap_fraction(bbox, boxes[o]):.0%}"
            for o in site_ids if o != s and bbox_overlap_fraction(bbox, boxes[o]) > 0
        ]
        if in_exp:
            mrds_note = (
                "flagged 'bbox?' for MRDS nudge; build log: kept as-is"
                if s in MRDS_NUDGE_FLAGGED else
                "listed for MRDS spot-check; no change recorded"
                if s in MRDS_SPOTCHECK_LISTED else ""
            )
            derivation = "literature_or_fixed"
            evidence = (
                "Hardcoded in scripts/expansion_add_sites.py SITES (research-subagent estimate of the "
                "district centre, fixed ~0.46 deg template box per docs/expansion_plan_100_sites.md); "
                "YAML written with a pre-registered prediction before any pipeline run; "
                f"YAML bbox {'identical to' if matches else 'DIFFERS from'} script; "
                f"{n_rev} bbox revision(s) in git."
            )
            if mrds_note:
                evidence += f" MRDS check: {mrds_note}."
        else:
            mrds_note = ""
            derivation = "unknown"
            evidence = PRE_EXPANSION_EVIDENCE.get(s, "No derivation recorded.")
            evidence += f" {n_rev} bbox revision(s) in git."
        rows.append({
            "site_id": s,
            "derivation_class": derivation,
            "bbox_wgs84": json.dumps([round(v, 4) for v in bbox]),
            "width_deg": round(bbox[2] - bbox[0], 3),
            "height_deg": round(bbox[3] - bbox[1], 3),
            "first_commit": first_commit,
            "n_bbox_revisions_git": n_rev,
            "in_expansion_script": in_exp,
            "mrds_check_note": mrds_note,
            "overlapping_boxes": "; ".join(overlaps),
            "evidence": evidence,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Task 2/3 -- box-independent sensitivity
# ---------------------------------------------------------------------------


def load_scene(site, paths):
    """Full-extent silica/carbonate/mafic ratios + B10, as run_classification loads them."""
    granule_id = resolve_granule_id(site, paths)
    if (paths.aster_dir / f"{granule_id}_ratio_silica.tif").exists():
        silica, carbonate, mafic, b10, transform, crs = load_ratio_mosaic(paths.aster_dir, granule_id)
    else:
        b10, _, b12, b13, b14, _, transform, crs = load_tir_bands_10_14(paths.aster_dir, granule_id)
        silica, carbonate, mafic = alteration_ratios(b12, b13, b14)
    return granule_id, silica, carbonate, mafic, b10, transform, crs


def run_variant(bands, transform, crs, cp, mrds, bbox=None):
    """Clip (optional) -> classify -> critical hit test. Returns (result, zones, footprint, transform, shape)."""
    if bbox is not None:
        bands, transform = clip_bands_to_bbox(list(bands), transform, crs, bbox)
    silica, carbonate, mafic, b10 = bands
    zones, footprint, raster_bbox = classify_strong_zones(
        silica, carbonate, mafic, b10, transform, crs, cp
    )
    res = critical_hit_test(mrds, zones, footprint, raster_bbox)
    res["valid_px"] = int(np.isfinite(b10).sum())
    return res, zones, footprint, transform, silica.shape


def _output_dir(site, paths) -> Path:
    # Both locations are covered by existing .gitignore data patterns.
    if site.layout == "nested":
        return paths.aster_dir.parent / "rasters" / "bbox_sensitivity"
    return paths.aster_dir / "bbox_sensitivity"


def cached_perturbations(site, paths, mrds) -> list[dict]:
    """Fallback: perturb the box over cached zones/footprint without reclassifying."""
    zones = gpd.read_file(paths.strong_zones_geojson)
    prov = json.loads(paths.site_provenance_json.read_text())
    footprint = gpd.GeoDataFrame.from_features(
        prov["tir_footprint_wgs84"]["features"], crs="EPSG:4326"
    ).to_crs(zones.crs)
    rows = []
    for name, vb in perturbation_variants(tuple(site.bbox_wgs84), PERTURB_FRAC).items():
        region = box(*_bbox_in_crs(vb, zones.crs))
        z, fp = restrict_to_region(zones, footprint, region)
        res = critical_hit_test(mrds, z, fp, vb)
        rows.append({"variant": name, "method": "cached_no_reclass", "bbox": vb, **res})
    return rows


def analyse_site(site_id: str, mrds: pd.DataFrame) -> list[dict]:
    site = load_site_by_id(site_id, REPO_ROOT / "sites")
    paths = site_paths_for(site, REPO_ROOT)
    cp = site.classification or ClassificationParams()
    granule_id, silica, carbonate, mafic, b10, transform, crs = load_scene(site, paths)
    bands = (silica, carbonate, mafic, b10)
    analyst_bbox = tuple(site.bbox_wgs84)
    rows: list[dict] = []

    def add(variant, method, bbox, res):
        rows.append({"site_id": site_id, "variant": variant, "method": method,
                     "bbox": None if bbox is None else [round(v, 4) for v in bbox], **res})

    res_a, _, _, a_transform, a_shape = run_variant(bands, transform, crs, cp, mrds, analyst_bbox)
    add("analyst_box", "reclassified", analyst_bbox, res_a)

    scene_valid = int(np.isfinite(b10).sum())
    full_scene = scene_valid >= FULL_SCENE_MIN_RATIO * res_a["valid_px"]

    if full_scene:
        res_w, zones_w, fp_w, _, _ = run_variant(bands, transform, crs, cp, mrds, None)
        add("whole_scene", "reclassified", None, res_w)

        out_dir = _output_dir(site, paths)
        out_dir.mkdir(parents=True, exist_ok=True)
        zones_w.to_file(out_dir / "whole_scene_zones.gpkg", driver="GPKG")

        region = fp_w.geometry.union_all().difference(raster_extent_polygon(a_transform, a_shape))
        z_out, fp_out = restrict_to_region(zones_w, fp_w, region)
        res_o = critical_hit_test(mrds, z_out, fp_out, raster_bbox_wgs84(transform, silica.shape, crs))
        add("whole_scene_outside_box", "reclassified_whole_scene", None, res_o)

        lon, lat = footprint_centroid_wgs84(fp_w)
        cbox = centred_box(lon, lat, TEMPLATE_BOX_DEG)
        res_c, *_ = run_variant(bands, transform, crs, cp, mrds, cbox)
        add("scene_centred_box", "reclassified", cbox, res_c)

        for name, vb in perturbation_variants(analyst_bbox, PERTURB_FRAC).items():
            res_p, *_ = run_variant(bands, transform, crs, cp, mrds, vb)
            add(name, "reclassified", vb, res_p)
    else:
        for r in cached_perturbations(site, paths, mrds):
            add(r.pop("variant"), r.pop("method"), r.pop("bbox"), r)

    for r in rows:
        r["granule_id"] = granule_id
        r["scene_to_box_valid_px_ratio"] = round(scene_valid / res_a["valid_px"], 2)
        r["whole_scene_available"] = full_scene
    return rows


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


def make_figure(df: pd.DataFrame, prov: pd.DataFrame) -> None:
    prov_tag = {"literature_or_fixed": "L", "mrds_informed": "M", "unknown": "U"}
    tags = dict(zip(prov["site_id"], prov["derivation_class"].map(prov_tag)))
    wide_e = df.pivot(index="site_id", columns="variant", values="enrichment")
    wide_p = df.pivot(index="site_id", columns="variant", values="p_binomial")
    perturb_cols = [c for c in wide_e.columns if c.startswith(("shift_", "shrink_", "expand_"))]

    order = wide_e["analyst_box"].sort_values().index.tolist()
    y = np.arange(len(order))
    p_cap = 1e-30

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(14, 11), sharey=True, gridspec_kw={"width_ratios": [1.25, 1]}
    )
    series = [
        ("analyst_box", "analyst box (reproduces published)", "o", "black"),
        ("whole_scene", "whole scene", "D", "#c0392b"),
        ("whole_scene_outside_box", "whole scene, outside analyst box", "^", "#e67e22"),
        ("scene_centred_box", f"{TEMPLATE_BOX_DEG}\u00b0 box at scene centre", "s", "#2471a3"),
    ]
    offsets = {"whole_scene": 0.27, "analyst_box": 0.09,
               "scene_centred_box": -0.09, "whole_scene_outside_box": -0.27}

    for ax, kind in ((ax1, "e"), (ax2, "p")):
        for i, s in enumerate(order):
            if perturb_cols:
                vals = wide_e.loc[s, perturb_cols] if kind == "e" else \
                    -np.log10(wide_p.loc[s, perturb_cols].clip(lower=p_cap))
                ax.fill_betweenx([i - 0.36, i + 0.36], vals.min(), vals.max(),
                                 color="0.85", lw=0, zorder=1)
        for key, label, marker, color in series:
            if key not in wide_e.columns:
                continue
            e = wide_e.loc[order, key].to_numpy(float)
            p = wide_p.loc[order, key].to_numpy(float)
            xv = e if kind == "e" else -np.log10(np.clip(p, p_cap, 1.0))
            sig = p < 0.05
            yy = y + offsets[key]
            ax.scatter(xv[sig], yy[sig], marker=marker, s=38, color=color,
                       edgecolor=color, zorder=3)
            ax.scatter(xv[~sig], yy[~sig], marker=marker, s=38, facecolor="white",
                       edgecolor=color, linewidth=1.4, zorder=3)
        for i in range(1, len(order)):
            ax.axhline(i - 0.5, color="0.93", lw=0.8, zorder=0)

    ax1.axvline(1.0, color="0.3", ls="--", lw=1)
    ax1.set_xlabel("Critical-mineral enrichment (hit rate / site-specific null p0)")
    ax1.set_xlim(left=0)
    ax2.axvline(-np.log10(0.05), color="0.3", ls="--", lw=1)
    ax2.set_xlabel("\u2212log10 one-sided binomial p  (capped at 30; dashed = p 0.05)")
    ax2.set_xlim(left=0)

    labels = [f"{s} [{tags.get(s, '?')}]" for s in order]
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels)
    for tick, s in zip(ax1.get_yticklabels(), order):
        if s in ROBUST_13:
            tick.set_fontweight("bold")
    ax1.set_ylim(-0.7, len(order) - 0.3)
    for ax in (ax1, ax2):
        ax.grid(axis="x", color="0.9")
        ax.set_axisbelow(True)

    handles = [Line2D([], [], marker=m, color=c, ls="", markersize=8, label=l)
               for _, l, m, c in series]
    handles.append(Line2D([], [], color="0.85", lw=8,
                          label=f"range over 8 box perturbations (\u00b1{int(PERTURB_FRAC*100)}%, reclassified)"))
    handles.append(Line2D([], [], marker="o", color="0.3", markerfacecolor="white", ls="",
                          markersize=8, label="hollow marker = p \u2265 0.05"))
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.0), fontsize=9.5)
    fig.suptitle(
        "Does enrichment depend on the analyst-drawn study box?  20 FDR-significant sites\n"
        "bold = robust-13;  box provenance: [L] literature/fixed template, [U] undocumented",
        fontsize=12, y=0.985,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.955))
    FIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_PATH, dpi=150)
    plt.close(fig)
    print(f"Saved {FIG_PATH}")


# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--sites", nargs="*", default=FDR_SIG_20)
    ap.add_argument("--figure-only", action="store_true",
                    help="re-render the figure from the existing result CSVs")
    args = ap.parse_args()
    site_ids = args.sites

    if args.figure_only:
        make_figure(pd.read_csv(RESULTS / "bbox_independence_sensitivity.csv"),
                    pd.read_csv(RESULTS / "bbox_provenance.csv"))
        return

    prov = build_provenance(site_ids)
    prov.to_csv(RESULTS / "bbox_provenance.csv", index=False)
    print(f"Saved {RESULTS / 'bbox_provenance.csv'}")
    print(prov[["site_id", "derivation_class", "n_bbox_revisions_git", "mrds_check_note",
                "overlapping_boxes"]].to_string(index=False))

    first = load_site_by_id(site_ids[0], REPO_ROOT / "sites")
    mrds = read_mrds_national(site_paths_for(first, REPO_ROOT))

    rows: list[dict] = []
    for s in site_ids:
        t0 = time.time()
        site_rows = analyse_site(s, mrds)
        rows.extend(site_rows)
        by = {r["variant"]: r for r in site_rows}
        a, w = by["analyst_box"], by.get("whole_scene")
        msg = f"{s:24s} box n={a['n_crit']:4d} E={a['enrichment']:.2f} p={a['p_binomial']:.1e}"
        if w:
            msg += f" | scene n={w['n_crit']:4d} E={w['enrichment']:.2f} p={w['p_binomial']:.1e}"
        print(f"{msg}  ({time.time() - t0:.0f}s)", flush=True)

    df = pd.DataFrame(rows)
    df["bbox"] = df["bbox"].apply(lambda b: "" if b is None else json.dumps(b))
    df["sig_05"] = df["p_binomial"] < 0.05
    df["robust_13"] = df["site_id"].isin(ROBUST_13)

    published = pd.read_csv(RESULTS / "site_specific_null_significance.csv")
    pub = published.set_index("site_id")
    a = df[df["variant"] == "analyst_box"].set_index("site_id")
    common = a.index.intersection(pub.index)
    mism = common[
        (a.loc[common, "n_crit"] != pub.loc[common, "n_crit"])
        | (a.loc[common, "hits_crit"] != pub.loc[common, "hits_crit"])
        | ((a.loc[common, "p0_null"] - pub.loc[common, "p0_null"]).abs() > 0.002)
    ]
    df["analyst_box_matches_published"] = df["site_id"].map(lambda s: s in common and s not in mism)
    if len(mism):
        print(f"\nWARNING analyst_box reproduction differs from published for: {list(mism)}")
    else:
        print(f"\nanalyst_box reproduction matches published n/hits/p0 for all {len(common)} sites")

    cols = ["site_id", "variant", "method", "n_crit", "hits_crit", "hit_rate_crit_pct", "p0_null",
            "expected_hits_null", "enrichment", "p_binomial", "sig_05", "footprint_km2", "zone_km2",
            "n_zones", "valid_px", "bbox", "robust_13", "whole_scene_available",
            "scene_to_box_valid_px_ratio", "analyst_box_matches_published", "granule_id"]
    df = df[cols]
    df.to_csv(RESULTS / "bbox_independence_sensitivity.csv", index=False)
    print(f"Saved {RESULTS / 'bbox_independence_sensitivity.csv'}\n")

    wide = df.pivot(index="site_id", columns="variant", values=["enrichment", "p_binomial"])
    show = [v for v in ("analyst_box", "whole_scene", "whole_scene_outside_box", "scene_centred_box")
            if v in wide["enrichment"].columns]
    print(wide.loc[:, (slice(None), show)].round(4).to_string())

    if "whole_scene" in df["variant"].values:
        ws = df[df["variant"] == "whole_scene"].set_index("site_id")
        r13 = ws[ws.index.isin(ROBUST_13)]
        print(f"\nRobust-13 significant (p<0.05) on whole scene: "
              f"{int(r13['sig_05'].sum())}/{len(r13)}  -> {sorted(r13.index[r13['sig_05']])}")
        print(f"All sites significant on whole scene: {int(ws['sig_05'].sum())}/{len(ws)}")

    make_figure(df, prov)


if __name__ == "__main__":
    main()
