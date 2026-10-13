"""
Finer-scale geologic-unit loaders beside USGS SGMC.

Provides ``load_fine_units`` for site-level lithology-conditioned nulls when
1:500k–1:1M SGMC is too coarse (e.g. skarn/CRD districts that collapse under
the 8-class SGMC null — see docs/data_verification.md §3.5).

Sources
-------
sgmc
    Existing ``load_sgmc_units`` path (FeatureServer layer 3).
nbmg_100k
    Nevada Bureau of Mines and Geology geology layers (MapServer). Tries
    denser layers first (24k → 100k → 62k → 48k) because coverage is patchy:
    Pioche has 24k features but zero 100k features in this service. Lithology
    from the ``lithology`` / ``description`` fields. **Caveat:** even 24k may
    cover only a fraction of a site bbox (Pioche ~19% of footprint); deposits
    outside mapped polygons fall in ``unmapped`` and the stratified null
    reverts toward geometric for those points.
tx_250k
    Texas Geologic Database 1:250k ``RockUnitPoly250K`` (MapServer/9), with
    attribute join to ``tblRockUnit250K`` (MapServer/11) for
    ``description`` / ``long_nm``. Keyword-mapped to SGMC GENERALIZED_LITH.
azgs_local
    Local GeoPackage or GeoJSON (e.g. AZGS OFR-20-02 Courtland–Gleeson).
    Requires ``geology_path=``. See module docstring section below.
auto
    Prefer a state-appropriate fine source when available; fall back to SGMC
    if the fine fetch returns no polygons.

Cache layout
------------
``data/lithology/{site_id}_{source}_units.geojson`` (same folder as SGMC).

Attribute → GENERALIZED_LITH mapping
------------------------------------
Downstream ``generalize_lithology`` / ``add_lith_class`` / ``build_strata``
expect an SGMC-style ``GENERALIZED_LITH`` string. Fine sources rarely carry
that field, so ``map_row_to_generalized_lith`` / ``normalize_units_lithology``
fill it:

1. If ``GENERALIZED_LITH`` is already non-empty, keep it.
2. Else scan text fields in order: lithology, LITH, LITHNAME, LITHNAME1,
   description, long_nm, UNIT_NAME, name, rockunit_cd (+ short_nm).
3. Keyword rules (case-insensitive substring) map common lithology words onto
   the same 33-value SGMC vocabulary used by ``LITH_CLASS_MAP`` in
   lithology.py. First match wins (more specific tokens listed first).

NBMG field note
---------------
Layer 17 exposes OneGeology strings in ``lithology`` (and
``geologicUnitType``). Values are free text (e.g. "limestone", "alluvium");
they are not the SGMC controlled vocabulary. Mapping is therefore
keyword-based and may leave uncommon units as ``Unknown`` → lith class
``other``.

Texas field note
----------------
Polygon layer 9 has only ``rockunit_cd`` / ``rockunit_id`` (no lithology).
Unit names and descriptions live in table layer 11; the loader joins on
``rockunit_id`` when the table query succeeds. If the join fails, only
``rockunit_cd`` keywords (e.g. Qal → alluvium heuristics) are available —
many Texas units then map to ``Unknown``.

Courtland–Gleeson (AZGS OFR-20-02)
----------------------------------
Download the GeoPackage from AZGS Open-File Report 20-02:

  https://library.azgs.arizona.edu/item/AOFR-1603324919493-818

Direct zip (``CourtlandGleeson_geopkg.zip``) is linked from that item page
(also mirrored at https://hdl.handle.net/10150/647725). Unzip into
``data/lithology/`` so the ``.gpkg`` is findable, then run::

    python scripts/lithology_matched_null.py --sites courtland_gleeson \\
        --geology-source azgs_local \\
        --geology-path data/lithology/CourtlandGleeson.gpkg

(Exact ``.gpkg`` filename inside the zip may vary; pass the real path.)
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd
from shapely.geometry import box
from shapely.validation import make_valid

from critical_minerals_aster.lithology import (
    SGMC_FIELDS,
    add_lith_class,
    load_sgmc_units,
)

NBMG_MAPSERVER = (
    "https://gisweb.unr.edu/nbmg/rest/services/Geology/NBMG_Geology/MapServer"
)
# Prefer denser layers: Pioche has 24k coverage but not 100k in this service.
NBMG_LAYER_CANDIDATES: tuple[tuple[str, int], ...] = (
    ("24k", 5),
    ("100k", 17),
    ("62k", 14),
    ("48k", 11),
)
NBMG_100K_URL = f"{NBMG_MAPSERVER}/17/query"
TX_250K_POLY_URL = (
    "https://feature.geographic.texas.gov/arcgis/rest/services"
    "/Geologic_Database/GeologicDatabaseofTexas_250k/MapServer/9/query"
)
TX_250K_TABLE_URL = (
    "https://feature.geographic.texas.gov/arcgis/rest/services"
    "/Geologic_Database/GeologicDatabaseofTexas_250k/MapServer/11/query"
)

GEOLOGY_SOURCES = frozenset({"sgmc", "nbmg_100k", "tx_250k", "azgs_local", "auto"})

SOURCE_LABELS = {
    "sgmc": (
        "USGS SGMC (Horton et al. 2017, doi:10.5066/F7WH2N65), "
        "FeatureServer layer 3 SGMC_Geology"
    ),
    "nbmg_100k": (
        "NBMG 1:100k geology (OneGeology), "
        "https://gisweb.unr.edu/nbmg/rest/services/Geology/NBMG_Geology/MapServer/17"
    ),
    "tx_250k": (
        "Texas Geologic Database 1:250k RockUnitPoly250K + tblRockUnit250K, "
        "https://feature.geographic.texas.gov/arcgis/rest/services/"
        "Geologic_Database/GeologicDatabaseofTexas_250k/MapServer"
    ),
    "azgs_local": "Local geology file (e.g. AZGS OFR-20-02 Courtland–Gleeson GeoPackage)",
}

# Text fields consulted (in order) when GENERALIZED_LITH is absent.
_LITH_TEXT_FIELDS = (
    "lithology",
    "LITH",
    "LITHNAME",
    "LITHNAME1",
    "description",
    "long_nm",
    "UNIT_NAME",
    "unit_name",
    "name",
    "short_nm",
    "rockunit_cd",
    "geologicUnitType",
)

# (substring, SGMC GENERALIZED_LITH) — first match wins; keep specific tokens first.
_KEYWORD_TO_SGMC: tuple[tuple[str, str], ...] = (
    ("serpentinite", "Metamorphic, serpentinite"),
    ("amphibolite", "Metamorphic, amphibolite"),
    ("granulite", "Metamorphic, granulite"),
    ("gneiss", "Metamorphic, gneiss"),
    ("schist", "Metamorphic, schist"),
    ("marble", "Metamorphic, carbonate"),
    ("quartzite", "Metamorphic, sedimentary clastic"),
    ("meta", "Metamorphic, undifferentiated"),
    ("granite", "Igneous, intrusive"),
    ("granodiorite", "Igneous, intrusive"),
    ("diorite", "Igneous, intrusive"),
    ("gabbro", "Igneous, intrusive"),
    ("monzonite", "Igneous, intrusive"),
    ("syenite", "Igneous, intrusive"),
    ("intrusive", "Igneous, intrusive"),
    ("pluton", "Igneous, intrusive"),
    ("rhyolite", "Igneous, volcanic"),
    ("dacite", "Igneous, volcanic"),
    ("andesite", "Igneous, volcanic"),
    ("basalt", "Igneous, volcanic"),
    ("tuff", "Igneous, volcanic"),
    ("ignimbrite", "Igneous, volcanic"),
    ("volcanic", "Igneous, volcanic"),
    ("lava", "Igneous, volcanic"),
    ("ash-flow", "Igneous, volcanic"),
    ("dolomite", "Sedimentary, carbonate"),
    ("limestone", "Sedimentary, carbonate"),
    ("carbonate", "Sedimentary, carbonate"),
    ("evaporite", "Sedimentary, evaporite"),
    ("gypsum", "Sedimentary, evaporite"),
    ("halite", "Sedimentary, evaporite"),
    ("chert", "Sedimentary, chemical"),
    ("sandstone", "Sedimentary, clastic"),
    ("siltstone", "Sedimentary, clastic"),
    ("mudstone", "Sedimentary, clastic"),
    ("shale", "Sedimentary, clastic"),
    ("conglomerate", "Sedimentary, clastic"),
    ("clastic", "Sedimentary, clastic"),
    ("alluvium", "Unconsolidated, undifferentiated"),
    ("alluvial", "Unconsolidated, undifferentiated"),
    ("colluvium", "Unconsolidated, undifferentiated"),
    ("playa", "Unconsolidated, undifferentiated"),
    ("dune", "Unconsolidated, undifferentiated"),
    ("loess", "Unconsolidated, undifferentiated"),
    ("unconsolidated", "Unconsolidated, undifferentiated"),
    ("gravel", "Unconsolidated, undifferentiated"),
    # Texas Quaternary codes often lack a prose description after a failed join.
    ("qal", "Unconsolidated, undifferentiated"),
    ("water", "Water"),
)


def fine_units_cache_path(site_id: str, source: str, repo_root: Path) -> Path:
    """``data/lithology/{site_id}_{source}_units.geojson`` under ``repo_root``."""
    return Path(repo_root) / "data" / "lithology" / f"{site_id}_{source}_units.geojson"


def preferred_fine_source(site_id: str, site_name: str | None = None) -> str:
    """Best non-SGMC source for a site from id/name heuristics; else ``sgmc``.

    NV (and NV/OR) → nbmg_100k; TX → tx_250k; everything else → sgmc.
    ``azgs_local`` is never auto-selected (requires an explicit local path).
    """
    name = (site_name or "").lower()
    sid = site_id.lower()
    blob = f"{sid} {name}"
    if (
        ", nv" in name
        or " nv/" in name
        or name.endswith(" nv")
        or "/nv" in name
        or "nevada" in blob
        or sid.endswith("_nv")
        or "_nv_" in sid
    ):
        return "nbmg_100k"
    if ", tx" in name or "texas" in blob or sid.endswith("_tx") or "_tx_" in sid:
        return "tx_250k"
    return "sgmc"


def map_row_to_generalized_lith(props: dict[str, Any] | pd.Series) -> str:
    """Map a feature's attributes to an SGMC ``GENERALIZED_LITH`` string."""
    get = props.get if isinstance(props, dict) else lambda k, default=None: (
        props[k] if k in props.index else default
    )
    existing = get("GENERALIZED_LITH")
    if existing is not None and str(existing).strip() and str(existing).strip().lower() not in (
        "nan",
        "none",
        "null",
    ):
        return str(existing).strip()

    parts: list[str] = []
    for field in _LITH_TEXT_FIELDS:
        val = get(field)
        if val is None or (isinstance(val, float) and pd.isna(val)):
            continue
        text = str(val).strip()
        if text and text.lower() not in ("nan", "none", "null"):
            parts.append(text)
    blob = " ".join(parts).lower()
    if not blob:
        return "Unknown"
    for needle, sgmc in _KEYWORD_TO_SGMC:
        if needle in blob:
            return sgmc
    return "Unknown"


def normalize_units_lithology(units: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Ensure ``GENERALIZED_LITH`` and ``lith_class`` columns for ``build_strata``."""
    out = units.copy()
    if len(out) == 0:
        if "GENERALIZED_LITH" not in out.columns:
            out["GENERALIZED_LITH"] = pd.Series(dtype=object)
        return add_lith_class(out) if "lith_class" not in out.columns else out
    out["GENERALIZED_LITH"] = [
        map_row_to_generalized_lith(row) for _, row in out.iterrows()
    ]
    return add_lith_class(out)


def _http_get_json(url: str, timeout: int = 60, max_retries: int = 4) -> dict:
    last_exc: Exception | None = None
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "critical-minerals-aster/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            if "error" in data:
                raise RuntimeError(f"ArcGIS error: {data['error']}")
            return data
        except Exception as exc:
            last_exc = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Geology query failed after {max_retries} attempts: {last_exc!r}")


def _fetch_arcgis_geojson(
    query_url: str,
    bbox_wgs84: tuple[float, float, float, float],
    *,
    out_fields: str = "*",
    page_size: int = 2000,
    pause_s: float = 0.25,
    order_by: str | None = "OBJECTID",
) -> dict:
    """Paged bbox GeoJSON query against an ArcGIS MapServer/FeatureServer layer."""
    west, south, east, north = bbox_wgs84
    features: list[dict] = []
    offset = 0
    while True:
        params: dict[str, Any] = {
            "where": "1=1",
            "geometry": f"{west},{south},{east},{north}",
            "geometryType": "esriGeometryEnvelope",
            "inSR": "4326",
            "outSR": "4326",
            "spatialRel": "esriSpatialRelIntersects",
            "outFields": out_fields,
            "resultOffset": offset,
            "resultRecordCount": page_size,
            "f": "geojson",
        }
        if order_by:
            params["orderByFields"] = order_by
        data = _http_get_json(query_url + "?" + urllib.parse.urlencode(params))
        page = data.get("features", [])
        features.extend(page)
        exceeded = data.get("exceededTransferLimit") or (data.get("properties") or {}).get(
            "exceededTransferLimit"
        )
        if not page or not exceeded:
            break
        offset += len(page)
        time.sleep(pause_s)
    return {"type": "FeatureCollection", "features": features}


def fetch_nbmg_units(
    bbox_wgs84: tuple[float, float, float, float],
    *,
    layers: tuple[tuple[str, int], ...] = NBMG_LAYER_CANDIDATES,
) -> gpd.GeoDataFrame:
    """Fetch NBMG geology polygons; try denser layers first until one has features.

    Coverage is patchy by scale: e.g. Pioche returns features on layer 5 (24k)
    but zero on layer 17 (100k). The returned frame carries
    ``attrs['nbmg_layer']`` / ``attrs['nbmg_scale']`` for provenance.
    """
    last_empty = gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    for scale, layer_id in layers:
        url = f"{NBMG_MAPSERVER}/{layer_id}/query"
        try:
            fc = _fetch_arcgis_geojson(url, bbox_wgs84, order_by="OBJECTID")
        except RuntimeError:
            try:
                fc = _fetch_arcgis_geojson(url, bbox_wgs84, order_by=None)
            except RuntimeError:
                continue
        gdf = _features_to_gdf(fc, bbox_wgs84)
        if len(gdf) == 0:
            last_empty = gdf
            continue
        gdf.attrs["nbmg_layer"] = layer_id
        gdf.attrs["nbmg_scale"] = scale
        return gdf
    return last_empty


def fetch_nbmg_100k_units(
    bbox_wgs84: tuple[float, float, float, float],
) -> gpd.GeoDataFrame:
    """Back-compat alias: densest available NBMG layer for the bbox."""
    return fetch_nbmg_units(bbox_wgs84)


def _fetch_tx_unit_table() -> pd.DataFrame:
    """Download tblRockUnit250K attribute rows (no geometry)."""
    rows: list[dict] = []
    offset = 0
    page_size = 2000
    while True:
        params = {
            "where": "1=1",
            "outFields": "rockunit_id,rockunit_cd,short_nm,long_nm,description",
            "resultOffset": offset,
            "resultRecordCount": page_size,
            "f": "json",
        }
        data = _http_get_json(TX_250K_TABLE_URL + "?" + urllib.parse.urlencode(params))
        page = data.get("features", [])
        for feat in page:
            rows.append(feat.get("attributes") or {})
        if not page or len(page) < page_size:
            break
        offset += len(page)
        time.sleep(0.25)
    if not rows:
        return pd.DataFrame(columns=["rockunit_id", "rockunit_cd", "short_nm", "long_nm", "description"])
    return pd.DataFrame(rows)


def fetch_tx_250k_units(
    bbox_wgs84: tuple[float, float, float, float],
) -> gpd.GeoDataFrame:
    """Fetch Texas 1:250k rock-unit polygons and join description table when possible."""
    try:
        fc = _fetch_arcgis_geojson(TX_250K_POLY_URL, bbox_wgs84, order_by="objectid")
    except RuntimeError:
        fc = _fetch_arcgis_geojson(TX_250K_POLY_URL, bbox_wgs84, order_by=None)
    gdf = _features_to_gdf(fc, bbox_wgs84)
    if len(gdf) == 0:
        return gdf
    # Normalize join key case.
    id_col = next((c for c in gdf.columns if c.lower() == "rockunit_id"), None)
    if id_col is None:
        return gdf
    try:
        table = _fetch_tx_unit_table()
    except Exception:
        return gdf
    if table.empty or "rockunit_id" not in table.columns:
        return gdf
    table = table.drop_duplicates(subset=["rockunit_id"])
    merged = gdf.merge(
        table,
        left_on=id_col,
        right_on="rockunit_id",
        how="left",
        suffixes=("", "_tbl"),
    )
    return gpd.GeoDataFrame(merged, geometry="geometry", crs=gdf.crs)


def load_azgs_local_units(
    geology_path: Path,
    bbox_wgs84: tuple[float, float, float, float],
) -> gpd.GeoDataFrame:
    """Read a local GeoPackage/GeoJSON and clip to ``bbox_wgs84``."""
    path = Path(geology_path)
    if not path.is_file():
        raise FileNotFoundError(f"azgs_local geology file not found: {path}")
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    else:
        gdf = gdf.to_crs("EPSG:4326")
    gdf["geometry"] = gdf.geometry.apply(lambda g: make_valid(g) if g is not None else g)
    gdf = gpd.clip(gdf, box(*bbox_wgs84))
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()]
    gdf = gdf[gdf.geom_type.isin(["Polygon", "MultiPolygon", "GeometryCollection"])]
    return gdf


def _features_to_gdf(
    fc: dict,
    bbox_wgs84: tuple[float, float, float, float],
) -> gpd.GeoDataFrame:
    if not fc.get("features"):
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    gdf = gpd.GeoDataFrame.from_features(fc["features"], crs="EPSG:4326")
    gdf["geometry"] = gdf.geometry.apply(lambda g: make_valid(g) if g is not None else g)
    gdf = gpd.clip(gdf, box(*bbox_wgs84))
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()]
    gdf = gdf[gdf.geom_type.isin(["Polygon", "MultiPolygon", "GeometryCollection"])]
    return gdf


def _write_cache(gdf: gpd.GeoDataFrame, cache: Path) -> None:
    cache.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(cache, driver="GeoJSON")


def _load_or_fetch(
    site_id: str,
    source: str,
    bbox_wgs84: tuple[float, float, float, float],
    repo_root: Path,
    *,
    refresh: bool,
    fetch_fn,
) -> gpd.GeoDataFrame:
    cache = fine_units_cache_path(site_id, source, repo_root)
    if cache.is_file() and not refresh:
        gdf = gpd.read_file(cache)
        return gdf.set_crs("EPSG:4326") if gdf.crs is None else gdf
    gdf = fetch_fn()
    if len(gdf):
        gdf = normalize_units_lithology(gdf)
        _write_cache(gdf, cache)
        return gdf
    return gpd.GeoDataFrame(
        columns=[*SGMC_FIELDS, "GENERALIZED_LITH", "lith_class", "geometry"],
        geometry="geometry",
        crs="EPSG:4326",
    )


def load_fine_units(
    site_id: str,
    bbox_wgs84: tuple[float, float, float, float],
    source: str,
    repo_root: Path,
    *,
    refresh: bool = False,
    geology_path: Path | str | None = None,
    site_name: str | None = None,
) -> gpd.GeoDataFrame:
    """Cache-first finer geology loader (EPSG:4326), clipped to ``bbox_wgs84``.

    Returns a GeoDataFrame with ``GENERALIZED_LITH`` and ``lith_class`` suitable
    for ``build_strata``. For ``source='auto'``, tries the preferred fine source
    for the site and falls back to SGMC when that source yields no polygons
    (empty fine results are not written to the fine-source cache).
    """
    source = (source or "sgmc").strip().lower()
    if source not in GEOLOGY_SOURCES:
        raise ValueError(f"Unknown geology source {source!r}; choose from {sorted(GEOLOGY_SOURCES)}")

    root = Path(repo_root)

    def _tag(gdf: gpd.GeoDataFrame, resolved: str) -> gpd.GeoDataFrame:
        gdf = gdf.copy()
        gdf.attrs["geology_source"] = resolved
        return gdf

    if source == "auto":
        preferred = preferred_fine_source(site_id, site_name)
        if preferred == "sgmc":
            return _tag(
                normalize_units_lithology(
                    load_sgmc_units(site_id, bbox_wgs84, root, refresh=refresh)
                ),
                "sgmc",
            )
        fine = load_fine_units(
            site_id,
            bbox_wgs84,
            preferred,
            root,
            refresh=refresh,
            geology_path=geology_path,
            site_name=site_name,
        )
        if len(fine) > 0:
            # Preserve resolved tag (e.g. nbmg_24k, sgmc_fallback_from_nbmg).
            return fine
        return _tag(
            normalize_units_lithology(
                load_sgmc_units(site_id, bbox_wgs84, root, refresh=refresh)
            ),
            "sgmc",
        )

    if source == "sgmc":
        return _tag(
            normalize_units_lithology(
                load_sgmc_units(site_id, bbox_wgs84, root, refresh=refresh)
            ),
            "sgmc",
        )

    if source == "azgs_local":
        if geology_path is None:
            raise ValueError("azgs_local requires geology_path= to a GeoPackage or GeoJSON")
        return _tag(
            _load_or_fetch(
                site_id,
                source,
                bbox_wgs84,
                root,
                refresh=refresh,
                fetch_fn=lambda: load_azgs_local_units(Path(geology_path), bbox_wgs84),
            ),
            "azgs_local",
        )

    if source == "nbmg_100k":
        fine = _load_or_fetch(
            site_id,
            source,
            bbox_wgs84,
            root,
            refresh=refresh,
            fetch_fn=lambda: fetch_nbmg_units(bbox_wgs84),
        )
        if len(fine) == 0:
            # Empty fine source must not become an all-unmapped "rescue".
            return _tag(
                normalize_units_lithology(
                    load_sgmc_units(site_id, bbox_wgs84, root, refresh=refresh)
                ),
                "sgmc_fallback_from_nbmg",
            )
        scale = fine.attrs.get("nbmg_scale", "unknown")
        return _tag(fine, f"nbmg_{scale}")

    if source == "tx_250k":
        fine = _load_or_fetch(
            site_id,
            source,
            bbox_wgs84,
            root,
            refresh=refresh,
            fetch_fn=lambda: fetch_tx_250k_units(bbox_wgs84),
        )
        if len(fine) == 0:
            return _tag(
                normalize_units_lithology(
                    load_sgmc_units(site_id, bbox_wgs84, root, refresh=refresh)
                ),
                "sgmc_fallback_from_tx_250k",
            )
        return _tag(fine, "tx_250k")

    raise ValueError(f"Unhandled geology source: {source}")


def unit_key_column(units: gpd.GeoDataFrame) -> str | None:
    """Best attribute for map-unit-level strata (finer than lith_class)."""
    if len(units) == 0:
        return None
    if "STATE" in units.columns and "SGMC_LABEL" in units.columns:
        return "sgmc_unit"
    for col in ("identifier", "rockunit_cd", "UNIT_NAME", "unit_name", "name", "SGMC_LABEL"):
        if col in units.columns:
            return col
    return None


def assign_unit_keys(units: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Add a ``unit`` column for map-unit stratified tests when possible."""
    out = units.copy()
    key = unit_key_column(out)
    if key is None:
        return out
    if key == "sgmc_unit":
        out["unit"] = out["STATE"].astype(str) + ":" + out["SGMC_LABEL"].astype(str)
    else:
        out["unit"] = out[key].astype(str)
    return out
