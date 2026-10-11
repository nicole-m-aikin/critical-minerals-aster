"""Retroactive cloud-cover audit of the ASTER L1T granules behind every site.

Granule selection (``spectral.score_granule`` / ``select_granule`` and the
mosaic candidate filter in ``pipeline.download_aster_mosaic``) ranked only by
bbox coverage and TIR band count, with no cloud screening.  This script
recovers the granule(s) each site was actually built from and pulls their
cloud metadata from NASA CMR, so cloud-contaminated inputs can be flagged.

Granule recovery, per site
--------------------------
1. ``results/{id}_provenance.json`` ``granule_id`` is a single AST_L1T id
   → that granule (``provenance``).
2. ``{id}_mosaic`` with per-band ``*_TIR_B10.tif`` files still in the site's
   aster dir → those granules, since the mosaic builder prefers local files
   (``local_files``).
3. ``{id}_mosaic`` without local granules (downloaded to a temp dir, then
   discarded) → replay the mosaic candidate search against CMR with the same
   parameters (AST_L1T, search bbox, temporal range, first 20 hits, coverage
   > 0.30, ≥3 TIR bands, bundle < max_bundle_mb) (``reconstructed``).

For 2 and 3 the union footprint of the recovered granules is compared with
the extent of the surviving ``{id}_mosaic_TIR_B10.tif`` (``mosaic_extent_check``)
to confirm the recovered set is the one the mosaic was built from.

Cloud metadata
--------------
* AST_L1T v004 UMM-G: scene ``CloudCover``, ``DayNightFlag``, acquisition time.
* AST_L1A v004 UMM-G (same acquisition): per-quadrant cloud coverage, scene
  orientation angle and rotated scene corners, from which an area-weighted
  cloud estimate over the site bbox is derived (``bbox_cloud_pct_est``).

Mosaics feather-blend every constituent granule into each overlapping pixel,
so a site-level ``site_blend_cloud_pct`` is also reported: the bbox-coverage-
weighted mean of per-granule cloud (bbox estimate where available, else
scene), approximating the cloud share of the blended signal.

Flags: ``ok`` < 10 %, ``review`` 10–30 %, ``high`` > 30 %, ``unknown``.

Outputs
-------
results/granule_cloud_audit.csv   — one row per (site, granule)

Usage
-----
    conda run -n aster-minerals python scripts/audit_granule_cloud_cover.py

Requires no EarthData credentials: CMR metadata search is public.  Raw CMR
responses are cached under data/cmr_cache/granule_cloud/ so re-runs are offline.
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from critical_minerals_aster.config import load_site_by_id, search_bbox  # noqa: E402
from critical_minerals_aster.paths import site_paths_for  # noqa: E402
from critical_minerals_aster.spectral import score_granule  # noqa: E402

SITES_DIR = ROOT / "sites"
RESULTS_DIR = ROOT / "results"
CACHE_DIR = ROOT / "data" / "cmr_cache" / "granule_cloud"
OUT_CSV = RESULTS_DIR / "granule_cloud_audit.csv"
CMR = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"

_GRANULE_ID_RE = re.compile(r"(AST_L1T_\d+_\d+)")

# Mirrors pipeline.download_aster_mosaic candidate filter.
MOSAIC_SEARCH_COUNT = 20
MOSAIC_MIN_COVERAGE = 0.30
MOSAIC_MIN_BANDS = 3
# Degrees; mosaic grid is snapped to the first granule's UTM pixel grid and
# L1T footprints are lon/lat rectangles, so exact equality is not expected.
EXTENT_TOL_DEG = 0.05

FDR_SIGNIFICANT = {
    "bisbee", "christmas", "contact_district", "courtland_gleeson", "crooks_gap",
    "eagle_mountain", "eureka", "holden", "karnes_county", "lake_valley",
    "lisbon_valley", "lordsburg", "magdalena_kelly", "mascot_jefferson_city",
    "organ_mountains", "providence_mountains", "shafter", "tombstone",
    "tri_state", "uravan",
}

# Visual inspection of figures/sites/{id}/01_tir_band_ratios.png and
# 02_classification.png for FDR-significant sites whose inputs were flagged.
VISUAL_CHECKS: dict[str, str] = {
    "tri_state": (
        "LIKELY AFFECTED: 16/20 blended granules >30% cloud. Ratio ranges very compressed "
        "(silica 0.927-0.934), pervasive scan striping, sharp-edged dark wedge in NW of "
        "coverage and broad non-geological bright arc mid-image that maps to the main "
        "combined-score≥3 cluster."
    ),
    "holden": (
        "POSSIBLE: ratios track topography (bright valley floors/Lake Chelan arm, dark "
        "ridges) rather than discrete cloud blobs; Mar-Jun alpine scenes, so high metadata "
        "cloud % may partly be snow. Strong zone follows the lake/valley, not cloud shapes."
    ),
    "mascot_jefferson_city": (
        "LIKELY AFFECTED: compact round sharp-edged patch at W edge (low silica, high "
        "carbonate, pale mafic) that classifies as carbonate class 2; other pale sharp-edged "
        "patches in mafic panel. Rest follows Valley-and-Ridge strike."
    ),
    "lisbon_valley": (
        "POSSIBLE: no discrete cloud blobs, but strong horizontal striping (all 12 night "
        "scenes) and a broad NW-SE regional gradient in all three ratios."
    ),
    "magdalena_kelly": (
        "LIKELY AFFECTED: lobate, crisp-edged cumulus-like blobs at S edge (dark silica, "
        "dark carbonate, white mafic) that classify as carbonate class 2 / elevated combined "
        "score; heavy striping throughout."
    ),
    "lake_valley": (
        "POSSIBLE: blotchy low-frequency patches (dark NE corner in silica+carbonate) and "
        "heavy striping; no crisp cloud outlines."
    ),
    "uravan": (
        "POSSIBLE: dark NE corner in silica and broad gradient plus striping; canyon "
        "network still dominates carbonate/mafic, no crisp cloud outlines."
    ),
    "eureka": (
        "LIKELY AFFECTED: several very dark, round, sharp-edged blobs (W-central and N-central) "
        "in silica with matching carbonate/mafic anomalies - classic cumulus shapes; W-central "
        "blob is carbonate+mafic class 2. Visible granule seam step in combined score."
    ),
    "tombstone": (
        "POSSIBLE: sharp oval (S-central) bright in silica/dark in mafic scores 4 in combined "
        "map; small dark round blobs NW in silica; rectangular tile edge in carbonate top. "
        "Could be cloud or real - check per-granule inputs."
    ),
    "karnes_county": (
        "UNLIKELY/INCONCLUSIVE: dense striping and speckled carbonate patches consistent with "
        "agricultural field mosaic; no crisp cloud outlines; dark SW corner."
    ),
    "courtland_gleeson": (
        "POSSIBLE: dark sharp-edged blob at top-centre in silica; rectangular field grids "
        "(real agriculture) in carbonate/mafic; large dark carbonate core looks geological."
    ),
    "crooks_gap": (
        "POSSIBLE: diffuse blotchy silica patches, vertical seam line and heavy striping "
        "(all night scenes); no crisp cloud outlines."
    ),
    "contact_district": (
        "UNLIKELY: ratios follow terrain/structure; one dark rounded patch in mafic (centre) "
        "but boundaries follow topography."
    ),
    "christmas": (
        "UNLIKELY: terrain-consistent; diagonal bright banding in silica looks like "
        "swath/feather seams rather than cloud."
    ),
    "shafter": (
        "UNLIKELY/INCONCLUSIVE: lobate bright-silica/pale-carbonate patch in SW quadrant with "
        "diffuse edges; otherwise terrain-consistent."
    ),
    "bisbee": (
        "UNLIKELY: single night granule (10% scene cloud); striping and scattered bright "
        "silica specks, no cloud-shaped features; one straight vertical edge NE in mafic."
    ),
}


def cloud_flag(pct: float | None) -> str:
    if pct is None or pd.isna(pct) or pct < 0:
        return "unknown"
    if pct < 10:
        return "ok"
    if pct <= 30:
        return "review"
    return "high"


def _cmr_get(params: list[tuple[str, str]], cache_path: Path) -> dict[str, Any]:
    if cache_path.is_file():
        return json.loads(cache_path.read_text())
    url = f"{CMR}?{urllib.parse.urlencode(params)}"
    last_exc: Exception | None = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            break
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            last_exc = exc
            time.sleep(2 ** attempt)
    else:
        raise RuntimeError(f"CMR request failed: {url}") from last_exc
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(payload))
    return payload


def fetch_l1t(granule_id: str) -> dict[str, Any] | None:
    payload = _cmr_get(
        [("short_name", "AST_L1T"), ("version", "004"),
         ("readable_granule_name", granule_id)],
        CACHE_DIR / f"{granule_id}__L1T.json",
    )
    items = payload.get("items") or []
    return items[0] if items else None


def fetch_l1a(granule_id: str) -> dict[str, Any] | None:
    """L1A record for the same acquisition (``AST_L1T_<17 digits>`` → ``AST_L1A_<17 digits>*``)."""
    acq = granule_id.split("_")[2]
    payload = _cmr_get(
        [("short_name", "AST_L1A"), ("version", "004"),
         ("readable_granule_name", f"AST_L1A_{acq}*"),
         ("options[readable_granule_name][pattern]", "true")],
        CACHE_DIR / f"{granule_id}__L1A.json",
    )
    items = payload.get("items") or []
    return items[0] if items else None


def search_site_candidates(site) -> list[dict[str, Any]]:
    """Replay the mosaic builder's earthaccess search as a raw CMR query."""
    lon0, lat0, lon1, lat1 = search_bbox(site)
    payload = _cmr_get(
        [("short_name", "AST_L1T"),
         ("bounding_box", f"{lon0},{lat0},{lon1},{lat1}"),
         ("temporal", f"{site.temporal_start}T00:00:00Z,{site.temporal_end}T23:59:59Z"),
         ("page_size", str(MOSAIC_SEARCH_COUNT))],
        CACHE_DIR / f"search__{site.id}.json",
    )
    return payload.get("items") or []


class UmmGranule:
    """Duck-types the earthaccess DataGranule surface used by spectral.score_granule."""

    def __init__(self, item: dict[str, Any]):
        self._item = item
        self._umm = item["umm"]

    def __getitem__(self, key: str) -> Any:
        if key == "umm":
            return self._umm
        return self._item[key]

    def data_links(self) -> list[str]:
        return [u["URL"] for u in self._umm.get("RelatedUrls", []) if u.get("Type") == "GET DATA"]

    def size(self) -> float:
        info = self._umm.get("DataGranule", {}).get("ArchiveAndDistributionInformation", [])
        return sum(float(s.get("SizeInBytes", 0)) for s in info) / (1024 * 1024)

    @property
    def granule_id(self) -> str:
        return self._umm.get("GranuleUR", "")


def _footprint(umm: dict[str, Any]) -> Polygon | None:
    try:
        pts = umm["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"]["GPolygons"][0][
            "Boundary"]["Points"]
    except (KeyError, IndexError, TypeError):
        return None
    poly = Polygon([(float(p["Longitude"]), float(p["Latitude"])) for p in pts])
    return poly if poly.is_valid else poly.buffer(0)


def _attr(umm: dict[str, Any], name: str) -> float | None:
    for a in umm.get("AdditionalAttributes", []) or []:
        if a.get("Name") == name:
            try:
                return float(a["Values"][0])
            except (KeyError, IndexError, TypeError, ValueError):
                return None
    return None


def _scene_quadrants(scene: Polygon, orientation_deg: float) -> dict[str, Polygon] | None:
    """Split an ASTER L1A scene footprint into UL/UR/LL/LR image quadrants.

    L1A footprints are the rotated scene corners.  Descending (day) passes have
    ASTERSceneOrientationAngle ≈ +10° (image-up ≈ north, image-left ≈ west);
    ascending (night) passes ≈ 170°, i.e. the image is rotated 180° and its
    upper-left corner is the south-east corner on the ground.
    """
    corners = list(dict.fromkeys(scene.exterior.coords))
    if len(corners) != 4:
        return None
    by_lat = sorted(corners, key=lambda c: c[1], reverse=True)
    nw, ne = sorted(by_lat[:2], key=lambda c: c[0])
    sw, se = sorted(by_lat[2:], key=lambda c: c[0])
    if abs(orientation_deg) > 90:
        ul, ur, ll, lr = se, sw, ne, nw
    else:
        ul, ur, ll, lr = nw, ne, sw, se

    def mid(a, b):
        return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)

    top, bottom, left, right = mid(ul, ur), mid(ll, lr), mid(ul, ll), mid(ur, lr)
    c = scene.centroid.coords[0]
    return {
        "ul": Polygon([ul, top, c, left]),
        "ur": Polygon([top, ur, right, c]),
        "lr": Polygon([c, right, lr, bottom]),
        "ll": Polygon([left, c, bottom, ll]),
    }


def granule_cloud_record(
    granule_id: str,
    bbox_wgs84: tuple[float, float, float, float],
    l1t_item: dict[str, Any] | None = None,
) -> dict[str, Any]:
    l1t = l1t_item or fetch_l1t(granule_id)
    l1a = fetch_l1a(granule_id)
    rec: dict[str, Any] = {
        "granule_id": granule_id,
        "acquisition_datetime": None,
        "day_night": None,
        "scene_cloud_pct": None,
        "quad_ul_pct": None,
        "quad_ur_pct": None,
        "quad_ll_pct": None,
        "quad_lr_pct": None,
        "bbox_cloud_pct_est": None,
        "bbox_coverage": None,
        "cmr_found": l1t is not None,
    }
    study = box(*bbox_wgs84)
    if l1t is not None:
        umm = l1t["umm"]
        temporal = umm.get("TemporalExtent", {})
        rec["acquisition_datetime"] = temporal.get("SingleDateTime") or temporal.get(
            "RangeDateTime", {}).get("BeginningDateTime")
        rec["day_night"] = umm.get("DataGranule", {}).get("DayNightFlag")
        cc = umm.get("CloudCover")
        rec["scene_cloud_pct"] = float(cc) if cc is not None else None
        fp = _footprint(umm)
        if fp is not None and study.area > 0:
            rec["bbox_coverage"] = round(study.intersection(fp).area / study.area, 3)

    if l1a is not None:
        umm_a = l1a["umm"]
        quads = {
            "ul": _attr(umm_a, "UpperLeftQuadCloudCoverage"),
            "ur": _attr(umm_a, "UpperRightQuadCloudCoverage"),
            "ll": _attr(umm_a, "LowerLeftQuadCloudCoverage"),
            "lr": _attr(umm_a, "LowerRightQuadCloudCoverage"),
        }
        for k, v in quads.items():
            rec[f"quad_{k}_pct"] = v if v is not None and v >= 0 else None
        if rec["day_night"] is None:
            rec["day_night"] = umm_a.get("DataGranule", {}).get("DayNightFlag")
        scene = _footprint(umm_a)
        orientation = _attr(umm_a, "ASTERSceneOrientationAngle")
        geoms = (
            _scene_quadrants(scene, orientation)
            if scene is not None and orientation is not None
            else None
        )
        if geoms and all(rec[f"quad_{k}_pct"] is not None for k in geoms):
            weighted, area = 0.0, 0.0
            for k, g in geoms.items():
                a = study.intersection(g).area
                weighted += a * rec[f"quad_{k}_pct"]
                area += a
            if area > 0:
                rec["bbox_cloud_pct_est"] = round(weighted / area, 1)
    return rec


def _mosaic_extent_check(aster_dir: Path, site_id: str, footprints: list[Polygon]) -> str:
    mosaic = aster_dir / f"{site_id}_mosaic_TIR_B10.tif"
    if not mosaic.is_file() or not footprints:
        return "no_mosaic_raster" if footprints else "no_footprints"
    import rasterio
    from rasterio.warp import transform_bounds

    with rasterio.open(mosaic) as ds:
        mb = transform_bounds(ds.crs, "EPSG:4326", *ds.bounds)
    ub = unary_union(footprints).bounds
    diff = max(abs(a - b) for a, b in zip(mb, ub))
    status = "match" if diff <= EXTENT_TOL_DEG else "mismatch"
    return f"{status} (max_diff_deg={diff:.3f})"


def recover_site_granules(site_id: str) -> tuple[list[tuple[str, dict | None]], str, str]:
    """Return ([(granule_id, l1t_item_or_None)], granule_source, mosaic_extent_check)."""
    site = load_site_by_id(site_id, SITES_DIR)
    paths = site_paths_for(site, ROOT)
    prov = json.loads(paths.site_provenance_json.read_text())
    gid = prov.get("granule_id") or ""

    if _GRANULE_ID_RE.fullmatch(gid):
        return [(gid, None)], "provenance", "not_mosaic"

    aster_dir = paths.aster_dir
    local = sorted({
        m.group(1)
        for f in (aster_dir.glob("*_TIR_B10.tif") if aster_dir.is_dir() else [])
        if "mosaic" not in f.name and (m := _GRANULE_ID_RE.search(f.name))
    })
    if local:
        items = [(g, fetch_l1t(g)) for g in local]
        fps = [fp for _, it in items if it and (fp := _footprint(it["umm"])) is not None]
        return items, "local_files", _mosaic_extent_check(aster_dir, site_id, fps)

    covering: list[tuple[str, dict]] = []
    for item in search_site_candidates(site):
        g = UmmGranule(item)
        try:
            coverage, _, band_count = score_granule(g, site.bbox_wgs84)
        except (KeyError, IndexError, TypeError, ValueError):
            continue
        if (
            coverage > MOSAIC_MIN_COVERAGE
            and band_count >= MOSAIC_MIN_BANDS
            and g.size() < site.max_bundle_mb
        ):
            covering.append((g.granule_id, item))
    fps = [fp for _, it in covering if (fp := _footprint(it["umm"])) is not None]
    return covering, "reconstructed", _mosaic_extent_check(aster_dir, site_id, fps)


def audit_site(site_id: str) -> list[dict[str, Any]]:
    site = load_site_by_id(site_id, SITES_DIR)
    granules, source, extent_check = recover_site_granules(site_id)
    if not granules:
        return [{
            "site_id": site_id, "granule_id": None, "granule_source": source,
            "mosaic_extent_check": extent_check, "flag": "unknown",
        }]
    rows = []
    for gid, item in granules:
        rec = granule_cloud_record(gid, site.bbox_wgs84, item)
        rows.append({
            "site_id": site_id,
            "fdr_significant": site_id in FDR_SIGNIFICANT,
            "granule_source": source,
            "n_granules": len(granules),
            "mosaic_extent_check": extent_check,
            **rec,
        })
    return rows


def main() -> None:
    site_ids = yaml.safe_load((SITES_DIR / "index.yaml").read_text())["sites"]
    with ThreadPoolExecutor(max_workers=8) as pool:
        per_site = list(pool.map(audit_site, site_ids))
    df = pd.DataFrame([r for rows in per_site for r in rows])
    df["fdr_significant"] = df["site_id"].isin(FDR_SIGNIFICANT)
    df["flag"] = df["scene_cloud_pct"].map(cloud_flag)
    df["bbox_flag"] = df["bbox_cloud_pct_est"].map(cloud_flag)

    best_cloud = df["bbox_cloud_pct_est"].fillna(df["scene_cloud_pct"])
    weight = df["bbox_coverage"].where(best_cloud.notna())
    blend = (
        (best_cloud * weight).groupby(df["site_id"]).sum(min_count=1)
        / weight.groupby(df["site_id"]).sum(min_count=1)
    ).round(1)
    df["site_blend_cloud_pct"] = df["site_id"].map(blend)
    df["site_flag"] = df["site_blend_cloud_pct"].map(cloud_flag)
    df["visual_check"] = df["site_id"].map(VISUAL_CHECKS).fillna("")

    cols = [
        "site_id", "fdr_significant", "granule_id", "granule_source", "n_granules",
        "mosaic_extent_check", "acquisition_datetime", "day_night", "scene_cloud_pct",
        "quad_ul_pct", "quad_ur_pct", "quad_ll_pct", "quad_lr_pct",
        "bbox_cloud_pct_est", "bbox_coverage", "flag", "bbox_flag",
        "site_blend_cloud_pct", "site_flag", "cmr_found", "visual_check",
    ]
    df = df.reindex(columns=cols)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)

    print(f"Wrote {OUT_CSV.relative_to(ROOT)}  ({len(df)} granule rows, "
          f"{df['site_id'].nunique()} sites)")
    print("\nGranule source:", df.drop_duplicates("site_id")["granule_source"]
          .value_counts().to_dict())
    ext = df.drop_duplicates("site_id")["mosaic_extent_check"].str.split(" ").str[0]
    print("Mosaic extent check:", ext.value_counts().to_dict())
    print("\nScene-level flag counts (granule rows):", df["flag"].value_counts().to_dict())
    print("Bbox-estimate flag counts (granule rows):", df["bbox_flag"].value_counts().to_dict())

    sites = (
        df.assign(
            _high=df["bbox_flag"].where(df["bbox_flag"] != "unknown", df["flag"]) == "high",
            _night=df["day_night"] == "Night",
        )
        .groupby("site_id")
        .agg(
            fdr=("fdr_significant", "first"),
            n=("n_granules", "first"),
            n_high=("_high", "sum"),
            n_night=("_night", "sum"),
            max_scene=("scene_cloud_pct", "max"),
            blend=("site_blend_cloud_pct", "first"),
            site_flag=("site_flag", "first"),
        )
        .sort_values(["fdr", "blend"], ascending=False)
    )
    print("\nSite-level flag (blended bbox cloud):", sites["site_flag"].value_counts().to_dict())

    def _print_sites(frame: pd.DataFrame) -> None:
        for sid, r in frame.iterrows():
            star = "*" if r["fdr"] else " "
            print(f" {star} {sid:<24} {r['site_flag']:<7} blend={r['blend']:>5}%  "
                  f"granules={r['n']:>2} (high-cloud {r['n_high']}, night {r['n_night']})  "
                  f"max_scene={r['max_scene']}%")

    flagged = sites[sites["site_flag"].isin(["review", "high"])]
    print(f"\nFlagged sites ({len(flagged)}; * = FDR-significant):")
    _print_sites(flagged)

    print("\nAll FDR-significant sites:")
    _print_sites(sites[sites["fdr"]])


if __name__ == "__main__":
    main()
