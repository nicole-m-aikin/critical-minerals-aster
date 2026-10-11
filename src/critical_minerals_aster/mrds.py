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
    # ── Skarn / Carbonate Replacement ───────────────────────────────────────
    # Magmatic-hydrothermal fluids replacing carbonate wallrock adjacent to a
    # porphyry-related intrusion. Genetically related to (and historically
    # lumped with) Porphyry Cu-Mo-Au below, but split out as its own system
    # here because skarn's carbonate-replacement mineralogy is a distinct,
    # strongly TIR-diagnostic surface expression (B13/B12 carbonate ratio)
    # that a deep porphyry stockwork alone does not produce. Checked BEFORE
    # Porphyry Cu-Mo-Au in this cascade so "skarn cu"-type model strings are
    # not absorbed by the broader porphyry rule. See docs/results.md
    # section 1.5 for the geological rationale and per-site cross-check.
    (
        "Skarn / Carbonate Replacement",
        ["skarn cu", "skarn zn-pb", "skarn fe", "skarn w", "w skarn",
         "polymetallic replacement"],
        ["skarn", "contact metasomatic", "contact metamorphic"],
        [],  # commod fallback deliberately empty — skarn commodities (Cu/Zn/Pb/Fe/W)
             # overlap heavily with porphyry and base-metal systems; forcing a
             # commodity-only guess here would be geologically unjustified.
    ),
    # ── Porphyry Cu-Mo-Au ───────────────────────────────────────────────────
    # Calc-alkaline arcs; porphyry stockwork core, distinct from its skarn
    # aureole (see Skarn / Carbonate Replacement above)
    (
        "Porphyry Cu-Mo-Au",
        ["porphyry cu", "porphyry mo", "epithermal quartz-alunite",
         "high sulfidation", "lithocap", "alkalic porphyry",
         "alkaline au-te", "epithermal vein, comstock",
         "epithermal vein, sado", "epithermal vein"],
        ["porphyry", "stockwork"],
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
    """Infer one of the 25 Earth MRI mineral systems for a single MRDS row.

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
    """Add mineral_system column mapping each MRDS row to one of the 24 OFR
    2020-1042 systems, plus a 25th "Skarn / Carbonate Replacement" category
    split out of Porphyry Cu-Mo-Au for this project (see docs/results.md — skarn's carbonate-replacement mineralogy is a distinct,
    strongly TIR-diagnostic surface expression that a bare porphyry stockwork
    does not share, and the two of them jointly explain most of the strongest
    site-level signals in the 45-site survey).

    Classification follows USGS OFR 2020-1042 v1.1 (Hofstra & Kreiner 2020)
    using a three-field cascade:  model > dep_type > commod1.  Missing fields
    are treated as empty strings and skipped.

    The 25 systems are:
        Placer, Chemical Weathering, Meteoric Recharge, Meteoric Convection,
        Lacustrine Evaporite, Marine Evaporite, Basin Brine Path,
        Marine Chemocline, Petroleum, Hybrid Magmatic REE, Arsenide,
        Volcanogenic Seafloor, Orogenic, Coeur d'Alene-type, Metamorphic,
        Skarn / Carbonate Replacement, Porphyry Cu-Mo-Au, Porphyry Sn,
        Reduced Intrusion-Related, Carlin-type, Climax-type, IOA-IOCG,
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


# ---------------------------------------------------------------------------
# Legacy-data validation (CIM 2018 §2.1.1.1: "validation" = checks for errors,
# mismatches, duplicates and mislabelling, distinct from "verification" of the
# underlying data, which lives in scripts/mrds_record_audit.py and
# docs/data_verification.md).
# ---------------------------------------------------------------------------

# USGS MRDS record-quality grade (ScienceBase 5d712a5ae4b0c4f70cfd4e54):
# A = most informative/reliable ... E = generally lacks references.
_SCORE_RANK = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}

# dev_stat == "Plant" is a processing facility (mill, smelter, plant), not a
# mineral occurrence; counting it as a "deposit" is a labelling error.
_NON_DEPOSIT_DEV_STAT = frozenset({"plant"})

_DEV_STAT_RANK = {
    "producer": 0,
    "past producer": 1,
    "prospect": 2,
    "occurrence": 3,
    "unknown": 4,
}


def validate_mrds(mrds: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Apply record-level validation checks to a raw MRDS table.

    Hard exclusions (each counted once, in this order):

    1. ``missing_coords`` — latitude or longitude null / non-numeric.
    2. ``out_of_range_coords`` — |lat| > 90 or |lon| > 180.
    3. ``non_deposit_plant`` — ``dev_stat == "Plant"`` (processing facility).
    4. ``duplicate_record`` — same normalised ``site_name`` at identical
       coordinates. The retained copy is the one with the best ``score`` grade,
       then the most advanced ``dev_stat``.

    Flags added (records kept):

    * ``qc_shared_coords`` — shares exact coordinates with another retained
      record (different name). Handled statistically by DBSCAN declustering,
      not dropped, because distinct workings are legitimately co-located.
    * ``qc_score_high`` — ``score`` in A–C.
    * ``qc_producer`` — ``dev_stat`` is Producer or Past Producer.

    The original index is preserved so downstream joins keyed on it still work.
    Coordinates are WGS84 per ``data/mrds.met`` (Horizontal_Datum_Name), so no
    datum shift is applied.
    """
    report: dict[str, int] = {"n_input": int(len(mrds))}
    df = mrds.copy()

    lat = pd.to_numeric(df["latitude"], errors="coerce")
    lon = pd.to_numeric(df["longitude"], errors="coerce")
    missing = lat.isna() | lon.isna()
    report["missing_coords"] = int(missing.sum())
    df = df.loc[~missing].copy()
    df["latitude"] = lat[~missing]
    df["longitude"] = lon[~missing]

    bad_range = (df["latitude"].abs() > 90) | (df["longitude"].abs() > 180)
    report["out_of_range_coords"] = int(bad_range.sum())
    df = df.loc[~bad_range]

    dev = df["dev_stat"].astype(str).str.strip().str.lower() if "dev_stat" in df.columns else pd.Series("", index=df.index)
    plant = dev.isin(_NON_DEPOSIT_DEV_STAT)
    report["non_deposit_plant"] = int(plant.sum())
    df = df.loc[~plant].copy()
    dev = dev.loc[df.index]

    score = df["score"].astype(str).str.strip().str.upper() if "score" in df.columns else pd.Series("", index=df.index)
    name = df["site_name"].astype(str).str.strip().str.lower() if "site_name" in df.columns else pd.Series("", index=df.index)
    order = pd.DataFrame(
        {
            "_name": name,
            "_lat": df["latitude"],
            "_lon": df["longitude"],
            "_score_rank": score.map(_SCORE_RANK).fillna(len(_SCORE_RANK)),
            "_dev_rank": dev.map(_DEV_STAT_RANK).fillna(len(_DEV_STAT_RANK)),
        },
        index=df.index,
    ).sort_values(["_score_rank", "_dev_rank"], kind="stable")
    dup = order.duplicated(["_name", "_lat", "_lon"], keep="first")
    dup_idx = order.index[dup]
    report["duplicate_record"] = int(len(dup_idx))
    df = df.drop(index=dup_idx)

    df["qc_shared_coords"] = df.duplicated(["latitude", "longitude"], keep=False)
    df["qc_score_high"] = score.loc[df.index].isin(["A", "B", "C"])
    df["qc_producer"] = dev.loc[df.index].isin(["producer", "past producer"])

    report["n_output"] = int(len(df))
    report["n_flag_shared_coords"] = int(df["qc_shared_coords"].sum())
    report["n_flag_score_high"] = int(df["qc_score_high"].sum())
    report["n_flag_producer"] = int(df["qc_producer"].sum())
    return df, report


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
    if not isinstance(deposits, gpd.GeoDataFrame):
        # Boolean filtering an empty GeoDataFrame with an object-dtype mask yields a DataFrame,
        # sometimes without the geometry column.
        if "geometry" not in deposits.columns:
            deposits = deposits.assign(geometry=None)
        deposits = gpd.GeoDataFrame(deposits, geometry="geometry", crs=zones.crs)
    joined = gpd.sjoin(deposits, zones, how="left", predicate="within")
    hits = joined[joined["index_right"].notna()]
    misses = joined[joined["index_right"].isna()]
    return joined, hits, misses
