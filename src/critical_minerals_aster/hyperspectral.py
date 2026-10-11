"""GEMx / AVIRIS hyperspectral coverage lookup for study sites.

Motivation
----------
``docs/results.md`` limitation 2 states that ASTER SWIR bands are unavailable
in LP DAAC v004 for these areas, so argillic / phyllic alteration is invisible
to the TIR method "regardless of any statistical correction", and that the
below-chance epithermal and sediment-hosted medians may partly reflect that
detection gap.  That is a stated but *untested* hypothesis.

USGS Earth MRI's GEMx campaign (Geological Earth Mapping Experiment, PIs
Kokaly and Hoefen, 2023-2027) flies imaging spectrometers over the arid
Southwest specifically to map critical-mineral alteration.  AVIRIS resolves
the 2.2 um Al-OH absorption features that identify kaolinite, illite,
muscovite and alunite -- exactly the minerals ASTER TIR cannot see.

Before any of that data can be used, one question has to be answered: which
of the 109 sites actually sit under a flight line?  Coverage is commonly
described state-by-state, but GEMx flies discrete numbered boxes (AZ01, NV36,
NM17 ...), so "site is in Nevada" says very little.  This module answers the
question quantitatively.

Approach
--------
Granule footprints are queried from NASA's Common Metadata Repository (CMR).
CMR metadata search needs no credentials, so the coverage report stays
runnable by anyone who clones the repo -- matching the zero-credential
``demo`` subcommand.  ``earthaccess`` is used elsewhere in this package for
*downloads*, which do require an EarthData login.

Two caveats the returned table makes explicit:

1. The AVIRIS collections are facility-wide archives spanning every AVIRIS
   campaign since 2006, not GEMx-only.  A bbox hit may be an old vegetation
   or greenhouse-gas flight.  Counts are therefore reported both unfiltered
   and restricted to the GEMx era (``GEMX_EPOCH`` onward).
2. Granule count is a poor proxy for usable coverage, because a flight line
   can clip a bbox corner.  ``coverage_frac`` is the union of GEMx-era
   footprints intersected with the site bbox, computed in an equal-area
   projection.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, box
from shapely.ops import unary_union

CMR_GRANULE_URL = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"

#: First GEMx deployment (NASA Airborne Science Program mission list).
GEMX_EPOCH = "2023-04-25"

#: Equal-area projection for CONUS, used for coverage fractions.
EQUAL_AREA_CRS = "EPSG:5070"

#: Whether a sensor resolves the ~2.2 um Al-OH features needed to identify
#: argillic / phyllic minerals by absorption-feature fitting.
SwirCapability = Literal["full", "partial", "none"]


@dataclass(frozen=True)
class Collection:
    """A CMR collection carrying imagery (or flight metadata) of interest."""

    key: str
    concept_id: str
    sensor: str
    #: True when the collection is exclusively GEMx flights, so a bbox hit
    #: needs no date filtering to be attributable to the campaign.
    gemx_only: bool
    swir: SwirCapability
    n_bands: int | None
    wavelength_nm: tuple[int, int] | None
    gsd_m: str
    note: str = ""


#: MASTER GEMx seasons are campaign-specific collections; the AVIRIS entries
#: are facility-wide archives that must be date-filtered to the GEMx era.
COLLECTIONS: tuple[Collection, ...] = (
    Collection(
        "master_gemx_2023", "C2905518711-ORNL_CLOUD", "MASTER",
        gemx_only=True, swir="partial", n_bands=50, wavelength_nm=(460, 12879),
        gsd_m="~50", note="GEMx CA-AZ, Summer 2023; L1B radiance",
    ),
    Collection(
        "master_gemx_2024", "C3255928163-ORNL_CLOUD", "MASTER",
        gemx_only=True, swir="partial", n_bands=50, wavelength_nm=(460, 12879),
        gsd_m="~50", note="GEMx Spring 2024; L1B radiance",
    ),
    Collection(
        "master_gemx_2025", "C3896679508-ORNL_CLOUD", "MASTER",
        gemx_only=True, swir="partial", n_bands=50, wavelength_nm=(460, 12879),
        gsd_m="~50", note="GEMx Summer 2025; L1B + L2 emissivity/LST",
    ),
    Collection(
        "master_gemx_2026", "C4271136560-ORNL_CLOUD", "MASTER",
        gemx_only=True, swir="partial", n_bands=50, wavelength_nm=(460, 12879),
        gsd_m="~50", note="GEMx Spring 2026; L1B radiance",
    ),
    Collection(
        "aviris5_l2a_rfl", "C4184090548-ORNL_CLOUD", "AVIRIS-5",
        gemx_only=False, swir="full", n_bands=425, wavelength_nm=(380, 2500),
        gsd_m="12-20", note="L2A orthocorrected surface reflectance, 5 nm sampling",
    ),
    # AVIRIS-3 (C3369603199-ORNL_CLOUD) is omitted on purpose. A full CMR
    # campaign census of its 511 granules is only GHG + 2025 LA Fires; it
    # flies a B-200 King Air, not the GEMx ER-2, and matched zero of the
    # 109 site bboxes in the GEMx era.
    Collection(
        "avirisc_l2_rfl", "C2711871294-ORNL_CLOUD", "AVIRIS-Classic",
        gemx_only=False, swir="full", n_bands=224, wavelength_nm=(400, 2500),
        gsd_m="~20", note="L2 calibrated reflectance (ER-2 high altitude)",
    ),
)

#: Collections whose SWIR sampling supports absorption-feature mineral ID.
SWIR_COLLECTION_KEYS = tuple(c.key for c in COLLECTIONS if c.swir == "full")


# ---------------------------------------------------------------------------
# CMR query
# ---------------------------------------------------------------------------


def _get_json(url: str, retries: int = 3, backoff: float = 1.5) -> dict:
    """GET *url* and parse JSON, retrying transient CMR failures."""
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                return json.loads(resp.read().decode())
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            last = exc
            if attempt < retries - 1:
                time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"CMR request failed after {retries} attempts: {url}") from last


def _granule_polygon(umm: dict) -> Polygon | None:
    """Extract a footprint polygon from a granule's UMM spatial extent.

    CMR granules express footprints either as ``GPolygons`` (the flight-line
    quadrilateral, used by MASTER and AVIRIS) or as ``BoundingRectangles``.
    Only the exterior ring is used; interior holes are not meaningful for a
    coverage fraction.
    """
    geom = (
        umm.get("SpatialExtent", {})
        .get("HorizontalSpatialDomain", {})
        .get("Geometry", {})
    )

    polys: list[Polygon] = []
    for gp in geom.get("GPolygons", []) or []:
        pts = gp.get("Boundary", {}).get("Points", []) or []
        if len(pts) >= 4:
            ring = [(p["Longitude"], p["Latitude"]) for p in pts]
            poly = Polygon(ring)
            if poly.is_valid and poly.area > 0:
                polys.append(poly)

    for br in geom.get("BoundingRectangles", []) or []:
        try:
            polys.append(
                box(br["WestBoundingCoordinate"], br["SouthBoundingCoordinate"],
                    br["EastBoundingCoordinate"], br["NorthBoundingCoordinate"])
            )
        except KeyError:
            continue

    if not polys:
        return None
    merged = unary_union(polys)
    return merged if isinstance(merged, Polygon) else merged.convex_hull


def _granule_start(umm: dict) -> str | None:
    t = umm.get("TemporalExtent", {})
    rng = t.get("RangeDateTime")
    if rng:
        return rng.get("BeginningDateTime")
    single = t.get("SingleDateTime")
    return single


def query_granules(
    concept_id: str,
    bbox: tuple[float, float, float, float],
    page_size: int = 2000,
    max_pages: int = 10,
) -> list[dict]:
    """Return granules from *concept_id* whose footprint intersects *bbox*.

    Each record carries ``granule_ur``, ``start`` (ISO timestamp or None) and
    ``geometry`` (shapely polygon or None).  Spatial filtering is done by CMR
    via its ``bounding_box`` parameter, so only matching granules transfer.
    """
    out: list[dict] = []
    for page in range(1, max_pages + 1):
        query = urllib.parse.urlencode({
            "collection_concept_id": concept_id,
            "bounding_box": ",".join(str(v) for v in bbox),
            "page_size": page_size,
            "page_num": page,
        })
        payload = _get_json(f"{CMR_GRANULE_URL}?{query}")
        items = payload.get("items", [])
        for item in items:
            umm = item.get("umm", {})
            out.append({
                "granule_ur": umm.get("GranuleUR", ""),
                "start": _granule_start(umm),
                "geometry": _granule_polygon(umm),
            })
        if len(items) < page_size:
            break
    return out


# ---------------------------------------------------------------------------
# Coverage geometry
# ---------------------------------------------------------------------------


def coverage_fraction(
    geometries: Iterable[Polygon | None],
    bbox: tuple[float, float, float, float],
) -> float:
    """Fraction of *bbox* covered by the union of *geometries* (0-1).

    Computed in ``EQUAL_AREA_CRS`` rather than degrees, because a degree of
    longitude spans meaningfully different ground distance across the 31-49 N
    latitude range of the survey.
    """
    polys = [g for g in geometries if g is not None and not g.is_empty]
    if not polys:
        return 0.0

    site = gpd.GeoSeries([box(*bbox)], crs="EPSG:4326").to_crs(EQUAL_AREA_CRS)
    site_area = float(site.iloc[0].area)
    if site_area <= 0:
        return 0.0

    union = gpd.GeoSeries([unary_union(polys)], crs="EPSG:4326").to_crs(EQUAL_AREA_CRS)
    overlap = float(union.iloc[0].intersection(site.iloc[0]).area)
    return min(overlap / site_area, 1.0)


def _in_gemx_era(start: str | None, epoch: str = GEMX_EPOCH) -> bool:
    """True when a granule's acquisition date falls on or after *epoch*."""
    if not start:
        return False
    return start[:10] >= epoch


# ---------------------------------------------------------------------------
# Per-site coverage
# ---------------------------------------------------------------------------


def site_coverage(
    site_id: str,
    bbox: tuple[float, float, float, float],
    collections: tuple[Collection, ...] = COLLECTIONS,
    cache_dir: Path | None = None,
) -> list[dict]:
    """One coverage record per collection for a single site.

    When *cache_dir* is given, raw CMR responses are cached per
    (site, collection) so repeat runs need no network access.
    """
    records: list[dict] = []
    for coll in collections:
        granules = _cached_query(coll, site_id, bbox, cache_dir)

        gemx_era = [g for g in granules if _in_gemx_era(g["start"])]
        # GEMx-only collections need no date filter; keep every granule so a
        # missing timestamp cannot silently zero out real coverage.
        usable = granules if coll.gemx_only else gemx_era
        starts = sorted(g["start"][:10] for g in usable if g["start"])

        records.append({
            "site_id": site_id,
            "collection": coll.key,
            "sensor": coll.sensor,
            "swir": coll.swir,
            "gemx_only": coll.gemx_only,
            "n_bands": coll.n_bands,
            "gsd_m": coll.gsd_m,
            "n_granules_all": len(granules),
            "n_granules_gemx_era": len(usable),
            "coverage_frac": round(
                coverage_fraction((g["geometry"] for g in usable), bbox), 4
            ),
            "first_date": starts[0] if starts else None,
            "last_date": starts[-1] if starts else None,
        })
    return records


def _cached_query(
    coll: Collection,
    site_id: str,
    bbox: tuple[float, float, float, float],
    cache_dir: Path | None,
) -> list[dict]:
    """``query_granules`` with an on-disk JSON cache keyed by site/collection."""
    if cache_dir is None:
        return query_granules(coll.concept_id, bbox)

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{site_id}__{coll.key}.json"

    if path.is_file():
        raw = json.loads(path.read_text())
        for rec in raw:
            wkt = rec.pop("geometry_wkt", None)
            if wkt:
                from shapely import wkt as shapely_wkt
                rec["geometry"] = shapely_wkt.loads(wkt)
            else:
                rec["geometry"] = None
        return raw

    granules = query_granules(coll.concept_id, bbox)
    path.write_text(json.dumps([
        {
            "granule_ur": g["granule_ur"],
            "start": g["start"],
            "geometry_wkt": g["geometry"].wkt if g["geometry"] is not None else None,
        }
        for g in granules
    ], indent=1))
    return granules


def build_coverage_table(
    sites: dict[str, tuple[float, float, float, float]],
    collections: tuple[Collection, ...] = COLLECTIONS,
    cache_dir: Path | None = None,
    progress: bool = True,
) -> pd.DataFrame:
    """Long-form coverage table: one row per (site, collection)."""
    rows: list[dict] = []
    total = len(sites)
    for i, (site_id, bbox) in enumerate(sites.items(), start=1):
        if progress:
            print(f"  [{i:>3}/{total}] {site_id}", flush=True)
        rows.extend(site_coverage(site_id, bbox, collections, cache_dir))
    return pd.DataFrame(rows)


def summarize_by_site(coverage: pd.DataFrame) -> pd.DataFrame:
    """Collapse the long coverage table to one decision row per site.

    ``swir_coverage_frac`` is the headline column: the best fractional bbox
    coverage from any sensor that can resolve the 2.2 um Al-OH features, and
    therefore the column that decides whether a site can be used to test the
    SWIR detection gap.
    """
    swir = coverage[coverage["collection"].isin(SWIR_COLLECTION_KEYS)]
    master = coverage[coverage["sensor"] == "MASTER"]

    def _agg(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
        grouped = frame.groupby("site_id")
        out = pd.DataFrame({
            f"{prefix}_coverage_frac": grouped["coverage_frac"].max(),
            f"{prefix}_granules": grouped["n_granules_gemx_era"].sum(),
        })
        best = (
            frame.sort_values("coverage_frac", ascending=False)
            .groupby("site_id")["sensor"]
            .first()
        )
        out[f"{prefix}_best_sensor"] = best
        return out

    site_ids = pd.Index(coverage["site_id"].unique(), name="site_id")
    out = pd.DataFrame(index=site_ids)
    out = out.join(_agg(swir, "swir")).join(_agg(master, "master"))
    out = out.fillna({
        "swir_coverage_frac": 0.0, "swir_granules": 0,
        "master_coverage_frac": 0.0, "master_granules": 0,
    })
    out["has_swir"] = out["swir_coverage_frac"] > 0
    out["swir_tier"] = pd.cut(
        out["swir_coverage_frac"],
        bins=[-0.001, 0.0, 0.25, 0.75, 1.0],
        labels=["none", "marginal", "partial", "good"],
    )
    return out.reset_index()


# ---------------------------------------------------------------------------
# Al-OH band depth — the SWIR feature ASTER TIR cannot see
# ---------------------------------------------------------------------------

# Clark-style single-feature continuum.  The 2.20 um Al-OH absorption is the
# diagnostic that identifies kaolinite, illite, muscovite and alunite; ASTER
# TIR ratios have no access to it.  Shoulders avoid the 1.9 um water and
# 2.3 um Mg-OH regions.
ALOH_FEATURE_NM = 2200.0
ALOH_CONTINUUM_NM = (2120.0, 2270.0)


def nearest_wavelength_index(wavelengths, target_nm: float) -> int:
    """Index of the band whose centre is closest to *target_nm*."""
    import numpy as np

    w = np.asarray(wavelengths, dtype=float)
    if w.size == 0:
        raise ValueError("wavelength array is empty")
    return int(np.argmin(np.abs(w - target_nm)))


def continuum_removed_band_depth(
    r_lo,
    r_feat,
    r_hi,
    w_lo: float,
    w_feat: float,
    w_hi: float,
):
    """Clark band depth: ``1 - R_feature / R_continuum``.

    The continuum is a straight line between the short- and long-wavelength
    shoulders.  Returns NaN where either reflectance is non-positive or the
    interpolated continuum is non-positive (nodata, water, shade).
    """
    import numpy as np

    if not (w_lo < w_feat < w_hi):
        raise ValueError(
            f"feature {w_feat} nm must lie between shoulders {w_lo} and {w_hi}"
        )
    t = (w_feat - w_lo) / (w_hi - w_lo)
    r_cont = r_lo + t * (r_hi - r_lo)
    with np.errstate(divide="ignore", invalid="ignore"):
        depth = 1.0 - np.asarray(r_feat, dtype=float) / r_cont
    valid = (
        (np.asarray(r_lo) > 0)
        & (np.asarray(r_feat) > 0)
        & (np.asarray(r_hi) > 0)
        & (r_cont > 0)
        & np.isfinite(depth)
    )
    return np.where(valid, depth, np.nan)


def aloh_band_depth(reflectance, wavelengths, nodata: float = -9999.0):
    """Al-OH absorption depth at ~2200 nm from a 3-D (band, y, x) cube.

    *reflectance* may also be (y, x, band); the wavelength axis is the one
    whose length matches ``len(wavelengths)``.
    """
    import numpy as np

    w = np.asarray(wavelengths, dtype=float)
    cube = np.asarray(reflectance, dtype=float)
    band_axis = int(np.argmin(np.abs(np.array(cube.shape) - w.size)))
    if cube.shape[band_axis] != w.size:
        raise ValueError(
            f"no axis of shape {cube.shape} matches {w.size} wavelengths"
        )
    cube = np.moveaxis(cube, band_axis, 0)
    cube = np.where(cube == nodata, np.nan, cube)

    i_lo = nearest_wavelength_index(w, ALOH_CONTINUUM_NM[0])
    i_feat = nearest_wavelength_index(w, ALOH_FEATURE_NM)
    i_hi = nearest_wavelength_index(w, ALOH_CONTINUUM_NM[1])
    return continuum_removed_band_depth(
        cube[i_lo], cube[i_feat], cube[i_hi],
        float(w[i_lo]), float(w[i_feat]), float(w[i_hi]),
    )


# Diagnostic feature positions (nm) after Clark / Swayze / Rowan & Mars.
# Kaolinite has a *pair* of Al-OH dips; alunite's pair is shifted a few nm
# longer.  A single 2200 nm depth cannot tell them apart; requiring *both*
# members of the pair is what makes the index mineral-specific-ish.
KAOLINITE_FEATURES_NM = (2165.0, 2206.0)
ALUNITE_FEATURES_NM = (2175.0, 2210.0)


def doublet_score(depth_a, depth_b):
    """Score that is high only when *both* absorption features are present.

    Uses the element-wise minimum: a pixel that has a deep 2.16 um dip but
    no 2.21 um dip (or the reverse) scores near zero.  That is the cheap
    version of "this looks like a kaolinite/alunite doublet, not just clay."
    """
    import numpy as np

    a = np.asarray(depth_a, dtype=float)
    b = np.asarray(depth_b, dtype=float)
    score = np.minimum(a, b)
    return np.where(np.isfinite(score) & (score > 0), score, np.nan)
