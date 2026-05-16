"""MRDS point filtering, CRS-safe joins to anomaly polygons, and Earth MRI classification."""

from __future__ import annotations

from typing import Any, Tuple

import geopandas as gpd
import pandas as pd

from critical_minerals_aster.config import BBox

# ---------------------------------------------------------------------------
# Earth MRI critical-mineral classification  (commodity-level)
# ---------------------------------------------------------------------------
# Based on USGS OFR 2020-1042 "Systems-Deposits-Commodities-Critical Minerals
# Table for the Earth Mapping Resources Initiative" and the 50-mineral 2022
# Final List of Critical Minerals (Federal Register v.87 no.37).
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
        # Lithium brine / evaporite / pegmatite — subsurface fluid or
        # sediment-hosted; no surface hydrothermal alteration footprint,
        # so TIR alteration mapping is NOT expected to detect these.
        "Battery Metals – Li/Brine",
        [
            "lithium",
            "cesium",
            "rubidium",
        ],
    ),
    (
        # Cobalt, nickel, manganese, graphite — often hosted in porphyry,
        # skarn, mafic, or hydrothermal systems whose alteration halos
        # ARE mappable by TIR.  Treated as a separate category so hit rates
        # don't get diluted by the Li/brine null result.
        "Battery Metals – Co/Ni",
        [
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
            # scandium – on 2022 list; byproduct across many systems
            "scandium",
            "rhenium",
            "selenium",
            # arsenic – on 2022 list; arsenide/Carlin systems
            "arsenic",
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
            # barite/barium → critical per 2022 list
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

# ---------------------------------------------------------------------------
# Earth MRI mineral-system classification  (geology-level)
# ---------------------------------------------------------------------------
# The 24 mineral systems defined in OFR 2020-1042 (v1.1, May 2021).
# Earth MRI focus areas (> 800 nationally; Dicken et al. 2022) are delineated
# by mineral system type, making this the natural unit for spatial analysis.
#
# Classification uses a three-field cascade for each MRDS row:
#   1. ``model``    – USGS deposit-model code string (most specific)
#   2. ``dep_type`` – deposit type vocabulary
#   3. ``commod1``  – primary commodity (broadest fallback)
# Rules are checked in priority order; the first match wins.
# Unmatched rows receive "Unknown System".
#
# Each entry: (system_name, model_kws, dep_type_kws, commod_kws)
# All keyword lists are lowercase substrings.
# ---------------------------------------------------------------------------

_SYSTEM_UNKNOWN = "Unknown System"

# fmt: off
_MINERAL_SYSTEM_RULES: list[tuple[str, list[str], list[str], list[str]]] = [
    # ── Placer ──────────────────────────────────────────────────────────────
    # Riverine/marine/eluvial concentration of resistate minerals
    (
        "Placer",
        ["placer", "shoreline placer", "alluvial placer", "beach placer",
         "stream placer"],
        ["placer", "stream placer", "beach placer", "alluvium"],
        [],  # commod fallback ambiguous – placer Au looks same as orogenic Au
    ),
    # ── Mafic Magmatic ──────────────────────────────────────────────────────
    # LIPs, mantle plumes; Ni-Cu-PGE sulfide + chromite + Fe-Ti oxide
    (
        "Mafic Magmatic",
        ["podiform chromite", "ni-cu-pge", "nickel-copper-pge",
         "cr-pge", "alaskan-type", "ilmenite", "anorthosite",
         "stillwater", "synorogenic-synvolcanic ni-cu",
         "lateritic ni", "ni-cu sulfide",
         "volcanic-hosted magnetite"],
        ["podiform", "podiform chromite"],
        ["chromite", "chromium"],
    ),
    # ── Chemical Weathering ─────────────────────────────────────────────────
    # Tropical laterite/bauxite; Ni-Co laterite; ion-adsorption REE
    (
        "Chemical Weathering",
        ["bauxite", "laterite", "ni-co laterite", "nickel laterite",
         "ion adsorption", "regolith", "supergene",
         "lateritic ni", "sedimentary kaolin"],
        ["residual", "secondary enrichment", "laterite", "supergene"],
        [],
    ),
    # ── Magmatic REE ────────────────────────────────────────────────────────
    # Carbonatite / peralkaline syenite / alkaline intrusions
    (
        "Magmatic REE",
        ["carbonatite", "thorium-rare-earth", "rare-earth vein",
         "peralkaline", "alkaline intrusion", "syenite ree",
         "alkalic ree", "nb-ree"],
        ["carbonatite", "alkalic", "peralkaline"],
        ["rare earth", "ree", "niobium", "columbium"],
    ),
    # ── Porphyry Cu-Mo-Au ───────────────────────────────────────────────────
    # Calc-alkaline arcs; broad spectrum porphyry → skarn → epithermal
    (
        "Porphyry Cu-Mo-Au",
        ["porphyry cu", "porphyry mo", "skarn cu", "skarn zn-pb",
         "skarn fe", "skarn w", "w skarn", "epithermal quartz-alunite",
         "high sulfidation", "lithocap", "alkalic porphyry",
         "alkaline au-te", "epithermal vein, comstock",
         "epithermal vein, sado", "epithermal vein",
         "polymetallic replacement", "replacement mn",
         "epithermal vein"],
        ["porphyry", "skarn", "contact metasomatic", "contact metamorphic",
         "stockwork", "replacement"],
        ["molybdenum"],
    ),
    # ── Climax-type ─────────────────────────────────────────────────────────
    # Continental rift; A-type topaz rhyolite; Mo-W-Sn + NYF pegmatites
    (
        "Climax-type",
        ["climax", "porphyry cu-mo", "topaz rhyolite", "nyf pegmatite",
         "volcanogenic u", "volcanogenic beryllium", "rhyolite tin",
         "greisen mo", "greisen sn"],
        [],
        [],
    ),
    # ── Porphyry Sn (granite-related) ───────────────────────────────────────
    # S-type peraluminous granites; LCT pegmatites → greisen → Sn porphyry
    (
        "Porphyry Sn",
        ["lct pegmatite", "be-li pegmatite", "sn-polymetallic",
         "porphyry sn", "granite-related sn", "alluvial placer sn",
         "mica pegmatite", "beryl pegmatite"],
        ["pegmatite"],
        ["tin", "cassiterite", "cesium", "tantalum", "niobium",
         "columbium", "beryllium"],
    ),
    # ── Reduced Intrusion-Related ────────────────────────────────────────────
    # Calc-alkaline arcs assimilating carbonaceous rocks; W-Au-Ag-Te-Bi
    (
        "Reduced Intrusion-Related",
        ["w veins", "w vein", "reduced intrusion", "tungsten vein",
         "greisen-v tungsten", "intrusion-related gold"],
        [],
        ["tungsten"],
    ),
    # ── IOA-IOCG ────────────────────────────────────────────────────────────
    # Subduction/rift; iron oxide-apatite + iron oxide-copper-gold
    (
        "IOA-IOCG",
        ["iocg", "ioa", "iron oxide copper gold", "iron oxide-apatite",
         "albitite uranium"],
        [],
        [],
    ),
    # ── Carlin-type ─────────────────────────────────────────────────────────
    # Continental arc; meteoric convection in carbonaceous sediments
    (
        "Carlin-type",
        ["carlin", "sediment-hosted au", "carbonate-hosted au",
         "disseminated au", "sedimentary au", "arsenic-thallium"],
        [],
        [],
    ),
    # ── Meteoric Convection ─────────────────────────────────────────────────
    # Mantle-plume / extensional volcanics; low-sulfidation epithermal
    (
        "Meteoric Convection",
        ["hot-spring au", "hot-spring hg", "hot-spring",
         "low sulfidation", "low-sulfidation",
         "epithermal mn", "fumarolic"],
        [],
        ["mercury"],
    ),
    # ── Orogenic ────────────────────────────────────────────────────────────
    # Metamorphic dewatering of sulfidic sequences; Au-quartz veins
    (
        "Orogenic",
        ["low-sulfide au-quartz", "low sulfide au", "orogenic au",
         "simple sb", "lode gold", "mesothermal"],
        ["metamorphic", "hydrothermal vein", "vein"],
        ["gold", "antimony"],
    ),
    # ── Coeur d'Alene-type ──────────────────────────────────────────────────
    # Metamorphic dewatering of oxidized siliciclastics; Ag-Pb-Zn
    (
        "Coeur d'Alene-type",
        ["coeur d'alene", "coeur dalene", "polymetallic vein",
         "polymetallic veins"],
        ["replacement vein"],
        [],
    ),
    # ── Volcanogenic Seafloor ───────────────────────────────────────────────
    # Spreading centers / back-arc; VMS / VHMS; Mn-Fe crusts
    (
        "Volcanogenic Seafloor",
        ["massive sulfide, kuroko", "massive sulfide, besshi",
         "massive sulfide, cyprus", "massive sulfide",
         "volcanogenic mn", "seafloor mn",
         "algoma fe", "stratabound exhalative",
         "sedimentary mn"],
        ["stratabound exhalative", "chemical sediment, marine",
         "bedded, chemical sediment"],
        [],
    ),
    # ── Basin Brine Path ────────────────────────────────────────────────────
    # Marine evaporite basins; MVT Zn-Pb; sedex; Cu±Co; unconformity U
    (
        "Basin Brine Path",
        ["mississippi valley", "mvt", "sediment-hosted cu",
         "sedimentary exhalative zn-pb", "sedex",
         "unconformity u", "breccia pipe u", "strontium replacement",
         "native cu", "basaltic cu", "volcanic redbed cu",
         "bedded barite", "bedded celestite", "bedded strontium"],
        ["stratabound", "stratiform", "manto",
         "sediment-hosted copper"],
        ["zinc", "lead"],
    ),
    # ── Marine Chemocline ───────────────────────────────────────────────────
    # Basin-brine discharge; black shales; phosphate; Mn-Fe oxides
    (
        "Marine Chemocline",
        ["phosphate, warm current", "phosphate, upwelling",
         "black shale", "metalliferous shale", "superior fe",
         "iron-manganese", "epithermal mn", "volcanogenic mn, cuban",
         "replacement mn"],
        ["chemical sediment, marine"],
        ["phosphat", "phosphorus", "manganese"],
    ),
    # ── Metamorphic ─────────────────────────────────────────────────────────
    # Recrystallization; graphite; magnesite; gneiss REE
    (
        "Metamorphic",
        ["graphite (coal", "graphite (amorphous", "graphite (flake",
         "metamorphic graphite", "crystalline graphite",
         "gneiss ree", "gneiss uranium"],
        ["metamorphic"],
        ["graphite"],
    ),
    # ── Marine Evaporite ────────────────────────────────────────────────────
    # Restricted epicontinental basins; potash; Mg; dissolution brine
    (
        "Marine Evaporite",
        ["marine evaporite", "bedded gypsum", "salt dome",
         "potash evaporite", "sabkha", "sedimentary magnesite",
         "salt-dome sulfur", "diapiric salt"],
        ["salt dome", "evaporite"],
        ["gypsum", "potash", "potassium"],
    ),
    # ── Lacustrine Evaporite ────────────────────────────────────────────────
    # Closed basins; Li brine/clay; borates; potash; evaporite sequence
    (
        "Lacustrine Evaporite",
        ["lacustrine borate", "lacustrine diatomite", "lacustrine gypsum",
         "lithium clay", "lithium-boron zeolite",
         "residual brine", "salar"],
        ["lacustrine"],
        ["lithium", "borate", "borax"],
    ),
    # ── Hybrid Magmatic REE / Basin Brine Path ──────────────────────────────
    # CO2/HF-bearing magmatic volatiles mixing with basinal brines; fluorspar
    (
        "Hybrid Magmatic REE",
        ["fluorspar deposit", "fluorite vein", "fluorspar vein",
         "illinois-kentucky fluorspar"],
        [],
        ["fluorspar", "fluorite"],
    ),
    # ── Meteoric Recharge ───────────────────────────────────────────────────
    # Oxidized meteoric groundwater in sandstone/carbonate aquifers; U-V
    (
        "Meteoric Recharge",
        ["sandstone u", "carbonate u", "calcrete u", "granite u",
         "tabular u", "roll-front u", "cryptocrystalline magnesite"],
        [],
        ["uranium", "vanadium"],
    ),
    # ── Arsenide ────────────────────────────────────────────────────────────
    # Continental rifts; oxidized basement brines; five-element veins
    (
        "Arsenide",
        ["five element", "five-element", "arsenide", "cobalt arsenide",
         "silver arsenide"],
        [],
        [],
    ),
    # ── Petroleum ───────────────────────────────────────────────────────────
    # Source rock + reservoir + seal; V/Ni in porphyrins; helium in gas
    (
        "Petroleum",
        ["petroleum", "natural gas", "oil shale", "tar sand",
         "helium gas"],
        [],
        ["petroleum", "natural gas", "helium", "oil shale"],
    ),
]
# fmt: on


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
    Initiative) and the 2022 Final List of Critical Minerals (50 minerals).

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
        Energy, REE, Battery Metals – Li/Brine, Battery Metals – Co/Ni,
        PGM, Base Metals, Specialty/High-Tech, Gold/Silver, Industrial.
    """
    return earth_mri_category in _CRITICAL_CATEGORIES


# ---------------------------------------------------------------------------
# Mineral-system classification  (geology-level)
# ---------------------------------------------------------------------------


def _classify_mineral_system(
    model: Any,
    dep_type: Any,
    commod1: Any,
) -> str:
    """Infer one of the 24 Earth MRI mineral systems for a single MRDS row.

    Uses a three-field cascade (most → least specific):
        1. ``model``    – USGS deposit-model code string
        2. ``dep_type`` – deposit type descriptor
        3. ``commod1``  – primary commodity

    Returns the first matching system name, or ``"Unknown System"`` if no
    rule matches.
    """
    m = str(model).lower() if pd.notna(model) and str(model).strip() else ""
    d = str(dep_type).lower() if pd.notna(dep_type) and str(dep_type).strip() else ""
    c = str(commod1).lower() if pd.notna(commod1) and str(commod1).strip() else ""

    for system, model_kws, dep_kws, commod_kws in _MINERAL_SYSTEM_RULES:
        if m and any(kw in m for kw in model_kws):
            return system
        if d and dep_kws and any(kw in d for kw in dep_kws):
            return system
        if c and commod_kws and any(kw in c for kw in commod_kws):
            return system
    return _SYSTEM_UNKNOWN


def reclassify_mrds_mineral_system(
    df: pd.DataFrame,
    model_col: str = "model",
    dep_type_col: str = "dep_type",
    commod_col: str = "commod1",
) -> pd.DataFrame:
    """Add mineral_system column mapping each MRDS row to one of the 24 Earth MRI systems.

    Classification follows USGS OFR 2020-1042 v1.1 (Hofstra & Kreiner 2020)
    using a three-field cascade:  model > dep_type > commod1.  Missing fields
    are treated as empty strings and skipped.

    The 24 systems are:
        Placer, Chemical Weathering, Meteoric Recharge, Meteoric Convection,
        Lacustrine Evaporite, Marine Evaporite, Basin Brine Path,
        Marine Chemocline, Petroleum, Hybrid Magmatic REE, Arsenide,
        Volcanogenic Seafloor, Orogenic, Coeur d'Alene-type, Metamorphic,
        Porphyry Cu-Mo-Au, Alkalic Porphyry, Porphyry Sn, Reduced
        Intrusion-Related, Carlin-type, Climax-type, IOA-IOCG,
        Magmatic REE, Mafic Magmatic  (plus "Unknown System").

    Parameters
    ----------
    df:
        MRDS DataFrame; must contain at least one of the three key columns.
    model_col, dep_type_col, commod_col:
        Column names to use (defaults match standard MRDS CSV headers).

    Returns
    -------
    Copy of *df* with a new ``mineral_system`` column.
    """
    out = df.copy()

    def _safe_col(col: str) -> pd.Series:
        if col in out.columns:
            return out[col]
        return pd.Series("", index=out.index)

    model_s = _safe_col(model_col)
    dep_s = _safe_col(dep_type_col)
    commod_s = _safe_col(commod_col)

    out["mineral_system"] = [
        _classify_mineral_system(m, d, c)
        for m, d, c in zip(model_s, dep_s, commod_s)
    ]
    return out


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
