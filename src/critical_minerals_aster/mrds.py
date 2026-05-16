"""MRDS point filtering, CRS-safe joins to anomaly polygons, and Earth MRI classification."""

from __future__ import annotations

from typing import Any, Tuple

import geopandas as gpd
import pandas as pd

from critical_minerals_aster.config import BBox

# ---------------------------------------------------------------------------
# Earth MRI critical-mineral classification
# ---------------------------------------------------------------------------
# Based on USGS OFR 2020-1042 "Systems-Deposits-Commodities-Critical Minerals
# Table for the Earth Mapping Resources Initiative" and the 50-mineral critical
# minerals list from the 2022 Executive Order 14017.
#
# Each tuple is (category_label, [lowercase_keyword, ...]).
# The list is checked in priority order; the first matching category wins.
# Anything not matched falls to "Non-Critical".
# ---------------------------------------------------------------------------

_EARTH_MRI_ORDERED: list[tuple[str, list[str]]] = [
    (
        "Energy",
        [
            "uranium",
            "thorium",
            "coal",
            "oil shale",
            "petroleum",
            "natural gas",
            "helium",
        ],
    ),
    (
        "REE",
        [
            "rare earth",
            # abbreviation – checked as substring; safe in MRDS vocabulary
            "ree",
            "cerium",
            "lanthanum",
            "neodymium",
            "praseodymium",
            "samarium",
            "europium",
            "gadolinium",
            "terbium",
            "dysprosium",
            "holmium",
            "erbium",
            "thulium",
            "ytterbium",
            "lutetium",
            "yttrium",
            "monazite",
            "xenotime",
        ],
    ),
    (
        "Battery Metals",
        [
            "lithium",
            "cobalt",
            "nickel",
            "manganese",
            "graphite",
        ],
    ),
    (
        "PGM",
        [
            "platinum",
            "palladium",
            "rhodium",
            "iridium",
            "osmium",
            "ruthenium",
        ],
    ),
    (
        "Base Metals",
        [
            "copper",
            "zinc",
            "lead",
            "molybdenum",
        ],
    ),
    (
        "Specialty/High-Tech",
        [
            "tungsten",
            "tin",
            "cassiterite",
            "antimony",
            "vanadium",
            "bismuth",
            "tellurium",
            "indium",
            "gallium",
            "germanium",
            "beryllium",
            "niobium",
            "columbium",
            "tantalum",
            "titanium",
            "zirconium",
            "hafnium",
            "scandium",
            "rhenium",
            "selenium",
        ],
    ),
    (
        "Gold/Silver",
        [
            "gold",
            "silver",
        ],
    ),
    (
        "Industrial",
        [
            # barite/barium → critical per 2022 EO
            "barite",
            "barium",
            # fluorspar → critical
            "fluorspar",
            "fluorite",
            "fluorine-fluorite",
            # magnesium → critical
            "magnesium",
            "magnesite",
            # potash → critical
            "potash",
            "potassium",
            # aluminum/bauxite → critical
            "aluminum",
            "bauxite",
            # chromium → critical
            "chromium",
            "chromite",
            # phosphate (fertilizer strategic mineral)
            "phosphat",
            "phosphorus",
        ],
    ),
]

_NON_CRITICAL = "Non-Critical"

_CRITICAL_CATEGORIES: frozenset[str] = frozenset(
    cat for cat, _ in _EARTH_MRI_ORDERED
)


def _classify_earth_mri(commod: Any) -> str:
    """Map a single commod1 string to an Earth MRI category label."""
    c = str(commod).lower()
    for category, keywords in _EARTH_MRI_ORDERED:
        if any(kw in c for kw in keywords):
            return category
    return _NON_CRITICAL


def reclassify_mrds_earth_mri(
    df: pd.DataFrame, commod_col: str = "commod1"
) -> pd.DataFrame:
    """Add earth_mri_category column mapping MRDS commod1 to Earth MRI critical mineral groups.

    Classification follows USGS OFR 2020-1042 (Earth Mapping Resources
    Initiative) and the 2022 Executive Order 50-mineral critical minerals list.

    Priority order (first match wins):
        Energy > REE > Battery Metals > PGM > Base Metals >
        Specialty/High-Tech > Gold/Silver > Industrial > Non-Critical

    Parameters
    ----------
    df:
        DataFrame that includes a commodity column (default ``commod1``).
    commod_col:
        Column name containing the primary commodity string.

    Returns
    -------
    Copy of *df* with a new ``earth_mri_category`` column.
    """
    out = df.copy()
    out["earth_mri_category"] = out[commod_col].apply(_classify_earth_mri)
    return out


def is_critical_mineral(earth_mri_category: str) -> bool:
    """Return True for Earth MRI categories that are critical minerals (not Non-Critical).

    Categories considered critical:
        Energy, REE, Battery Metals, PGM, Base Metals, Specialty/High-Tech,
        Gold/Silver, Industrial.
    """
    return earth_mri_category in _CRITICAL_CATEGORIES


def filter_mrds_bbox(mrds: pd.DataFrame, bbox: BBox) -> pd.DataFrame:
    """Subset MRDS rows whose longitude/latitude fall inside bbox (WGS84)."""
    lon0, lat0, lon1, lat1 = bbox
    return mrds.loc[
        (mrds["longitude"] >= lon0)
        & (mrds["longitude"] <= lon1)
        & (mrds["latitude"] >= lat0)
        & (mrds["latitude"] <= lat1)
    ].copy()


def mrds_to_points_gdf(mrds_local: pd.DataFrame, target_crs) -> gpd.GeoDataFrame:
    """Build point GeoDataFrame from MRDS lat/lon; reproject to target CRS (e.g. zones.crs)."""
    gdf = gpd.GeoDataFrame(
        mrds_local,
        geometry=gpd.points_from_xy(mrds_local.longitude, mrds_local.latitude),
        crs="EPSG:4326",
    )
    return gdf.to_crs(target_crs)


def spatial_join_deposits_zones(
    deposits: gpd.GeoDataFrame,
    zones: gpd.GeoDataFrame,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Join deposits to zones with predicate 'within'. Returns joined, hits, misses."""
    joined = gpd.sjoin(deposits, zones, how="left", predicate="within")
    hits = joined[joined["index_right"].notna()]
    misses = joined[joined["index_right"].isna()]
    return joined, hits, misses
