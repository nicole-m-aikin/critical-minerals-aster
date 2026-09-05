"""Generate site YAMLs + index + pre-registered category predictions for the +58 expansion.

Consolidated from five parallel research subagents (see docs/results.md).
Idempotent: skips any site whose YAML already exists, any index entry already present,
and any results/deposit_system_categories.csv row already present.

Run once, from the repo root:  python scripts/expansion_add_sites.py

After this: nudge bboxes if needed, then
  python -m critical_minerals_aster run --site <id> --mosaic   (batched)
then the corrected significance chain (README.md, Methods).
"""
from __future__ import annotations

import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SITES_DIR = REPO / "sites"
INDEX = SITES_DIR / "index.yaml"
CATEGORIES_CSV = REPO / "results" / "deposit_system_categories.csv"

# site_id, name, bbox_wgs84, category, confidence, prediction, rationale
# prediction: one of "null", "positive", "positive-ish", "null-to-weak"
SITES: list[tuple] = [
    # --- Tier 1a: Uranium/Energy  (current n=1; pre-registered expectation: NULL) ---
    ("shirley_basin", "Shirley Basin, WY", [-106.38, 41.97, -105.92, 42.43],
     "Uranium/Energy", "high", "null",
     "Classic Wyoming sandstone roll-front U in Eocene Wind River Fm arkose; the type area "
     "where the roll-front oxidation-reduction model was refined. No intrusion-related alteration."),
    ("crooks_gap", "Crooks Gap-Green Mountain, WY", [-108.03, 42.07, -107.57, 42.53],
     "Uranium/Energy", "high", "null",
     "Roll-front and fault-controlled U in Eocene Battle Spring Fm sandstone S of Jeffrey City; "
     "a separate district from Gas Hills."),
    ("ambrosia_lake", "Ambrosia Lake, NM", [-108.10, 35.19, -107.64, 35.65],
     "Uranium/Energy", "high", "null",
     "Largest U district in the US; elongate tabular 'manto' ore in Jurassic Morrison Fm "
     "(Westwater Canyon sandstone), Grants Mineral Belt."),
    ("uravan", "Uravan Mineral Belt, CO", [-108.96, 38.09, -108.50, 38.55],
     "Uranium/Energy", "high", "null",
     "Carnotite tabular U-V bodies in the Salt Wash Member of the Morrison Fm; the type belt "
     "for Colorado Plateau tabular V-U ore (Gateway-Uravan-Slick Rock)."),
    ("henry_mountains", "Henry Mountains, UT", [-110.88, 37.62, -110.42, 38.08],
     "Uranium/Energy", "high", "null",
     "Tabular V-U in fluvial Salt Wash Member sandstone on the SE flank of the Henry Mtns "
     "(Tony M / Bullfrog deposits near Ticaboo)."),
    ("lisbon_valley", "Lisbon Valley, UT", [-109.51, 38.04, -109.05, 38.50],
     "Uranium/Energy", "high", "null",
     "Utah's largest U district; ore in the basal Moss Back Member of the Triassic Chinle Fm "
     "along the SW flank of the Lisbon Valley salt anticline."),
    ("karnes_county", "Karnes County, TX", [-98.13, 28.72, -97.67, 29.18],
     "Uranium/Energy", "high", "null",
     "Core of the South Texas roll-front belt; marginal-marine roll-front U (coffinite) in the "
     "Upper Eocene Jackson Group (Panna Maria, Hobson)."),
    ("arizona_strip", "Arizona Strip, AZ", [-112.83, 36.52, -112.37, 36.98],
     "Uranium/Energy", "medium", "null",
     "GENETICALLY DISTINCT SUB-TYPE: high-grade U in solution-collapse breccia pipes rooted in "
     "the Mississippian Redwall Limestone (Pigeon, Kanab North, Hack, Pinenut). Box on the "
     "northern Kanab Creek pipe cluster."),

    # --- Tier 1b: Mafic/Ultramafic  (current n=1; pre-registered expectation: NULL) ---
    ("duluth_partridge_river", "Duluth Complex - Partridge River, MN", [-92.23, 47.35, -91.77, 47.79],
     "Mafic/Ultramafic", "high", "null",
     "Basal-contact disseminated Cu-Ni-PGE sulfide of the 1.1 Ga Partridge River intrusion "
     "(NorthMet); humid boreal forest and wetland. Mafic protolith, no hydrothermal alteration, "
     "heavy vegetation cover."),
    ("duluth_south_kawishiwi", "Duluth Complex - South Kawishiwi, MN", [-91.93, 47.63, -91.47, 48.07],
     "Mafic/Ultramafic", "high", "null",
     "Basal magmatic Cu-Ni-PGE sulfide of the South Kawishiwi intrusion SE of Ely; humid "
     "boreal terrain."),
    ("webster_addie", "Webster-Addie, NC", [-83.41, 35.14, -82.95, 35.58],
     "Mafic/Ultramafic", "high", "null",
     "Alpine-type ultramafic ring (dunite-websterite-pyroxenite) with accessory chromite, "
     "southern Appalachian ophiolitic belt; humid forested Blue Ridge."),
    ("josephine_ophiolite", "Josephine Ophiolite, OR", [-123.91, 41.93, -123.45, 42.37],
     "Mafic/Ultramafic", "high", "null",
     "Podiform chromite in harzburgite/dunite of the Late Jurassic Josephine ophiolite mantle "
     "section; sparsely vegetated serpentine barrens, summer-dry Klamath Mtns."),
    ("new_idria", "New Idria, CA", [-120.90, 36.18, -120.44, 36.62],
     "Mafic/Ultramafic", "high", "null",
     "Coast Range ophiolite serpentinite diapir with podiform chromite (plus Cr-Hg-asbestos); "
     "large barren serpentine exposure in the semi-arid southern Diablo Range."),
    ("state_line_chromite", "PA-MD State Line Chromite, PA", [-76.36, 39.50, -75.90, 39.94],
     "Mafic/Ultramafic", "high", "null",
     "Podiform chromite in serpentinized Alpine peridotite of the Baltimore mafic complex; "
     "sparsely vegetated serpentine barrens in humid SE PA / Cecil Co., MD."),

    # --- Tier 1c: Sediment-hosted  (current n=3; pre-registered expectation: NULL) ---
    ("getchell", "Getchell Trend, NV", [-117.45, 41.02, -116.99, 41.48],
     "Sediment-hosted", "high", "null",
     "Premier Carlin-type Au trend along the Getchell/Range Front fault, E flank of the Osgood "
     "Mtns; micron gold in arsenian pyrite in Paleozoic carbonate. No surface alteration cap."),
    ("cortez", "Cortez-Pipeline, NV", [-117.01, 40.02, -116.55, 40.48],
     "Sediment-hosted", "high", "null",
     "Multi-Moz Carlin-type disseminated Au (Pipeline, Cortez Hills) in lower-plate Roberts "
     "Mtns / Wenban carbonate, Battle Mountain-Eureka trend."),
    ("alligator_ridge", "Alligator Ridge, NV", [-115.68, 39.52, -115.22, 39.98],
     "Sediment-hosted", "high", "null",
     "Shallow Carlin-type Au in Pilot Shale / Pogonip carbonate, discovered 1976 with no prior "
     "mining (a useful discovery-bias datapoint)."),
    ("long_canyon", "Long Canyon, NV", [-114.81, 40.62, -114.35, 41.08],
     "Sediment-hosted", "medium", "null",
     "Off-trend Carlin-type discovery on the E flank of the Pequop Mtns; sub-5 um gold in "
     "oxidized arsenian pyrite in Cambrian Notch Peak dolomite."),
    ("tri_state", "Tri-State District, MO-KS-OK", [-94.93, 36.77, -94.47, 37.23],
     "Sediment-hosted", "high", "null",
     "Textbook Mississippi Valley-Type Pb-Zn: sphalerite-galena in strata-bound solution "
     "breccias in Mississippian cherty limestone (Boone Fm)."),
    ("upper_mississippi_valley", "Upper Mississippi Valley, WI-IL-IA", [-90.53, 42.32, -90.07, 42.78],
     "Sediment-hosted", "high", "null",
     "The classic MVT district; stratabound Zn-Pb in Ordovician Galena/Platteville dolostone, "
     "Grant-Lafayette-Iowa-Jo Daviess Cos."),
    ("mascot_jefferson_city", "Mascot-Jefferson City, TN", [-83.83, 35.87, -83.37, 36.33],
     "Sediment-hosted", "high", "null",
     "Largest carbonate-hosted MVT Zn district in the eastern US; sphalerite in strata-bound "
     "breccia bodies in Lower Ordovician Knox Group dolostone, NE of Knoxville."),
    ("austinville_ivanhoe", "Austinville-Ivanhoe, VA", [-81.16, 36.58, -80.70, 37.04],
     "Sediment-hosted", "high", "null",
     "Strata-bound MVT Zn-Pb in the Cambrian Shady Dolomite on a carbonate-ramp platform "
     "margin, Wythe Co."),
    ("black_butte", "Black Butte, MT", [-111.13, 46.57, -110.67, 47.03],
     "Sediment-hosted", "high", "null",
     "Sedex-style Cu-(Co-Ag): laterally continuous massive-sulfide + barite horizons in "
     "dolomitic black shale of the Mesoproterozoic Newland Fm, Belt Basin."),

    # --- Tier 1d: Arid skarn/CRD replication  (pre-registered expectation: POSITIVE) ---
    ("magdalena_kelly", "Magdalena (Kelly), NM", [-107.44, 33.87, -106.98, 34.32],
     "Skarn/Carbonate Replacement", "high", "positive",
     "Zn-Pb skarn, replacement and vein orebodies in the Mississippian Kelly Limestone around "
     "Tertiary intrusions; garnet-pyroxene skarn + smithsonite in a bare desert range. "
     "Independent replicate of the Bisbee/Eureka regime outside SE Arizona."),
    ("lake_valley", "Lake Valley, NM", [-107.80, 32.50, -107.34, 32.95],
     "Skarn/Carbonate Replacement", "high", "positive",
     "Bridal Chamber-style Ag-halide + Mn replacement mantos in Mississippian Lake Valley "
     "Limestone, ~20 mi from the Hillsboro Laramide porphyry; well-exposed desert foothills, "
     "closely comparable to Tombstone."),
    ("organ_mountains", "Organ Mountains, NM", [-106.79, 32.16, -106.33, 32.61],
     "Skarn/Carbonate Replacement", "high", "positive",
     "Cu-Pb-Zn-Ag sulfide skarn and replacement bodies in Paleozoic limestone against the "
     "Organ quartz-monzonite batholith; steep, essentially bare granite-carbonate range."),
    ("victorio", "Victorio Mountains, NM", [-108.53, 32.00, -108.07, 32.45],
     "Skarn/Carbonate Replacement", "high", "positive",
     "W-Be-Mo skarn/tactite plus carbonate-hosted Pb-Zn replacement and porphyry Mo in "
     "Paleozoic dolostone/limestone; a full skarn-CRD-porphyry system in a near soil-free "
     "desert range."),
    ("san_francisco_frisco", "San Francisco (Frisco), UT", [-113.45, 38.22, -112.99, 38.67],
     "Skarn/Carbonate Replacement", "high", "positive",
     "Horn Silver fault-controlled Ag-Pb replacement pipe/manto in Paleozoic carbonate with "
     "associated Cu-W skarn (Cactus, Frisco Contact); box also spans the adjacent Star Range "
     "skarn/CRD. Well-exposed Basin & Range replicate."),
    ("shafter", "Shafter (Presidio), TX", [-104.55, 29.65, -104.09, 30.10],
     "Skarn/Carbonate Replacement", "high", "positive",
     "Bedding- and fault-controlled Ag replacement orebodies in Permian Mina Grande limestone "
     "on the flank of the Chinati caldera; a fully independent Chihuahuan-Desert replicate."),
    ("cherry_creek", "Cherry Creek, NV", [-115.13, 39.68, -114.67, 40.13],
     "Skarn/Carbonate Replacement", "medium", "positive",
     "Intrusion-related scheelite tactite (garnet-pyroxene skarn) plus Ag-Pb carbonate "
     "replacement, jasperoid and decalcification in Paleozoic carbonate, Egan Range - same "
     "mineralogy as Eureka."),
    ("providence_mountains", "Providence Mountains, CA", [-115.76, 34.72, -115.30, 35.17],
     "Skarn/Carbonate Replacement", "medium", "positive",
     "Lenticular intrusion-related Ag-Pb (galena) replacement orebodies in the Cambrian "
     "Bonanza King Formation; superb bare-rock desert exposure in the eastern Mojave."),
    ("eagle_mountain", "Eagle Mountain, CA", [-115.72, 33.65, -115.26, 34.10],
     "Skarn/Carbonate Replacement", "medium", "positive",
     "Jurassic magnetite-pyrite Fe skarn - extensive anhydrous calc-silicate skarn where "
     "quartz-monzonite sills intrude Paleozoic carbonate/quartzite; fully exposed desert "
     "range (former Kaiser Steel pit). The cleanest independent Fe-skarn test."),
    ("contact_district", "Contact District, NV", [-114.99, 41.53, -114.53, 41.98],
     "Porphyry Cu-Mo-Au", "medium", "positive",
     "Porphyry-Cu system with Cu skarn and calc-silicate hornfels where monzonite-granodiorite "
     "intrudes Paleozoic limestone; bare high-desert contact aureole. An independent analog of "
     "the Bisbee porphyry-to-replacement pairing (skarn-fringe style)."),

    # --- Tier 2a: VMS  (current n=2, both anti-correlated; pre-registered expectation: NULL) ---
    ("west_shasta", "West Shasta, CA", [-122.76, 40.45, -122.30, 40.90],
     "Volcanogenic/VMS", "high", "null",
     "Kuroko-type Cu-Zn-Au massive sulfide in the Devonian Balaklala Rhyolite, an exposed "
     "submarine felsic-arc sequence, semi-arid Klamath Mtns. Tests whether an exposed VMS "
     "system can show any TIR signal at all."),
    ("east_shasta", "East Shasta, CA", [-122.36, 40.56, -121.90, 41.00],
     "Volcanogenic/VMS", "high", "null",
     "Kuroko-type Cu-Zn(-Ag-Pb-Au) VMS in Triassic felsic volcanic rocks NE of Redding; same "
     "belt as West Shasta, well exposed."),
    ("big_mike", "Big Mike, NV", [-117.78, 40.33, -117.32, 40.77],
     "Volcanogenic/VMS", "high", "null",
     "Small Cu(-Zn) VMS pod in mafic volcanic rocks of the Havallah sequence, exposed in arid "
     "Pershing County - a minimal-vegetation end-member for the VMS test."),
    ("iron_king", "Iron King, AZ", [-112.47, 34.29, -112.01, 34.73],
     "Volcanogenic/VMS", "high", "null",
     "Proterozoic bimodal-mafic VMS: Au-Ag-Cu-Pb-Zn massive sulfide in metarhyolite/"
     "meta-andesite of the Yavapai arc, semi-arid Bradshaw Mtns."),
    ("holden", "Holden, WA", [-121.01, 47.98, -120.55, 48.42],
     "Volcanogenic/VMS", "high", "null",
     "Kuroko-type Cu-Zn-Au VMS in Cascade metavolcanics; large exposed gossan/tailings in a "
     "montane rain-shadowed setting E of the Cascade crest above Lake Chelan."),
    ("gossan_lead", "Great Gossan Lead, VA", [-80.95, 36.53, -80.49, 36.97],
     "Volcanogenic/VMS", "high", "null",
     "Besshi-type stratiform pyrrhotite-Cu-Zn massive sulfide, Blue Ridge; a humid vegetated "
     "Appalachian analog to Ducktown (VMS genesis + vegetation suppression)."),
    ("ore_knob", "Ore Knob, NC", [-81.63, 36.20, -81.17, 36.64],
     "Volcanogenic/VMS", "high", "null",
     "Recrystallized massive pyrrhotite-chalcopyrite-pyrite VMS in Carolina gneiss; humid "
     "vegetated Blue Ridge, largest NC copper producer."),

    # --- Tier 2b: Humid skarn controls  (pre-registered expectation: NULL - vegetation control) ---
    ("cornwall_pa", "Cornwall, PA", [-76.64, 40.05, -76.16, 40.50],
     "Skarn/Carbonate Replacement", "high", "null",
     "Type locality of 'Cornwall-type' Fe skarn: massive magnetite replacing Cambrian carbonate "
     "against Triassic diabase, diopside-garnet-actinolite-phlogopite calc-silicate gangue - "
     "same genetic family as Bisbee - but in humid (~110 cm/yr) mixed-deciduous-forest / "
     "farmland Piedmont. VEGETATION-SUPPRESSION CONTROL: genesis matches the best-performing "
     "category, terrain does not."),
    ("french_creek_pa", "French Creek - Grace Mine, PA", [-76.06, 39.96, -75.60, 40.41],
     "Skarn/Carbonate Replacement", "high", "null",
     "Two Cornwall-type magnetite skarns (Grace Mine, French Creek Mines at St. Peters) in "
     "humid forested Berks/Chester Co. woodland. Vegetation-suppression control."),
    ("dillsburg_pa", "Dillsburg, PA", [-77.23, 39.89, -76.77, 40.34],
     "Skarn/Carbonate Replacement", "high", "null",
     "Cornwall-type magnetite from 30+ pits where Triassic limestone is contact-metasomatized "
     "between two diabase bodies; humid, wooded/agricultural York Co. Piedmont. Control."),
    ("willsboro_lewis_ny", "Willsboro-Lewis, NY", [-73.70, 44.10, -73.24, 44.55],
     "Skarn/Carbonate Replacement", "high", "null",
     "Garnet-clinopyroxene-wollastonite skarn: marble replaced by infiltration metasomatism "
     "along the margin of the Marcy anorthosite massif; dense Adirondack northern-hardwood/"
     "conifer forest, heavy snowpack. Control."),
    ("tilly_foster_ny", "Tilly Foster, NY", [-73.85, 41.18, -73.39, 41.63],
     "Skarn/Carbonate Replacement", "medium", "null",
     "Magnesian (dolomite-replacement) Fe skarn - magnetite with serpentine, chondrodite, "
     "brucite, clinochlore - in the forested Hudson Highlands, humid continental, thick "
     "vegetation and glacial cover. Control."),
    ("snoqualmie_wa", "Snoqualmie, WA", [-121.66, 47.20, -121.18, 47.65],
     "Skarn/Carbonate Replacement", "medium", "null",
     "Cu-Fe skarn (Guye iron-skarn, andradite-diopside skarn) where the Oligo-Miocene "
     "Snoqualmie batholith intrudes Triassic metasedimentary rock at Snoqualmie Pass; wet "
     "west-central Cascade crest, dense Douglas-fir/hemlock forest, >250 cm/yr precip. Control."),

    # --- Tier 3a: Epithermal, exposed cap  (pre-registered expectation: POSITIVE-ISH) ---
    ("summitville", "Summitville, CO", [-106.823, 37.211, -106.373, 37.661],
     "Epithermal", "high", "positive-ish",
     "Textbook volcanic-dome-hosted high-sulfidation Au; concentrically zoned vuggy silica -> "
     "quartz-alunite -> quartz-kaolinite -> argillic -> propylitic alteration fully exposed in "
     "the South Mountain dome. Exposed-cap test of the Phase 6 below-chance epithermal result."),
    ("round_mountain", "Round Mountain, NV", [-117.295, 38.485, -116.845, 38.935],
     "Epithermal", "high", "positive-ish",
     "Giant low-sulfidation quartz-adularia hot-springs Au deposit in Round Mountain Tuff; "
     "propylitic-to-advanced-argillic alteration halo well exposed in and around the open pit."),
    ("republic", "Republic, WA", [-118.965, 48.425, -118.515, 48.875],
     "Epithermal", "high", "positive-ish",
     "Eocene low-sulfidation Au-Ag veins in the Republic graben (Sanpoil Volcanics), capped by "
     "preserved hot-spring sinter deposits marking a near-intact paleosurface."),
    ("jarbidge", "Jarbidge, NV", [-115.655, 41.645, -115.205, 42.095],
     "Epithermal", "high", "positive-ish",
     "Mid-Miocene low-sulfidation quartz-adularia veins in rhyolite with well-preserved "
     "bladed-calcite-replacement textures and ~3000 ft vertical vein extent - a shallow, "
     "minimally eroded system."),

    # --- Tier 3b: Alkaline/Carbonatite  (pre-registered expectation: NULL-TO-WEAK) ---
    ("wet_mountains", "Wet Mountains, CO", [-105.695, 38.045, -105.245, 38.495],
     "Alkaline/Carbonatite", "high", "null-to-weak",
     "Ediacaran-Cambrian mafic-ultramafic + nepheline-syenite complexes (McClure Mtn, Gem "
     "Park) cut by carbonatite, lamprophyre and syenite dikes plus Th-REE fracture veins; "
     "three exposed complexes within the box."),
    ("round_top", "Round Top Mountain, TX", [-105.699, 31.052, -105.249, 31.502],
     "Alkaline/Carbonatite", "high", "null-to-weak",
     "Eocene peralkaline rhyolite laccolith (~2 km across) hosting disseminated Y-fluorite "
     "HREE mineralization plus Be-Li-F; well exposed above the Diablo Plateau."),
    ("cornudas", "Cornudas Mountains, NM", [-105.745, 31.795, -105.295, 32.245],
     "Alkaline/Carbonatite", "medium", "null-to-weak",
     "Tertiary nepheline-syenite / phonolite laccolith-and-dike cluster (Wind Mountain) with "
     "REE-Nb-Zr (+/- Be-Li) in the basal syenite and alkaline dikes."),
    ("lemitar", "Lemitar Mountains, NM", [-107.175, 33.945, -106.725, 34.395],
     "Alkaline/Carbonatite", "medium", "null-to-weak",
     "~100+ Cambrian calcite-dolomite carbonatite dikes (bastnasite-fluorite-barite) intruding "
     "Proterozoic basement with fenitization, Lemitar + Chupadera Mtns; REE up to ~1.1%."),

    # --- Tier 3c: Porphyry Cu-Mo-Au  (pre-registered expectation: POSITIVE-ISH) ---
    ("safford_az", "Safford, AZ", [-109.85, 32.75, -109.39, 33.20],
     "Porphyry Cu-Mo-Au", "high", "positive-ish",
     "Laramide porphyry Cu cluster (Dos Pobres, San Juan, Lone Star) with thick supergene/"
     "exotic Cu, N-NE of Safford; active Freeport operation, distinct from Morenci."),
    ("san_manuel_az", "San Manuel-Kalamazoo, AZ", [-110.86, 32.38, -110.40, 32.83],
     "Porphyry Cu-Mo-Au", "high", "positive-ish",
     "Textbook porphyry Cu-Mo in Ruin Granite, a single orebody tilted and split by "
     "mid-Tertiary normal faulting into the San Manuel and blind Kalamazoo halves; classic "
     "alteration/metal zoning studies."),
    ("tyrone_nm", "Tyrone, NM", [-108.58, 32.44, -108.12, 32.89],
     "Porphyry Cu-Mo-Au", "high", "positive-ish",
     "Laramide porphyry Cu with a well-developed supergene chalcocite blanket SW of Silver "
     "City; a separate system from Santa Rita/Chino and from hanover_fierro (~25 km NE)."),
    ("questa_nm", "Questa, NM", [-105.76, 36.48, -105.30, 36.93],
     "Porphyry Cu-Mo-Au", "high", "positive-ish",
     "Climax-type porphyry Mo deposit in the Latir volcanic field, Taos Co.; stockwork "
     "molybdenite in and above a Tertiary granite porphyry, semi-arid mountainous terrain."),
]

YAML_TEMPLATE = """\
# {name} -- {category}
# {rationale_wrapped}
# Added in the +58 expansion (docs/results.md). Pre-registered expectation
# BEFORE running the pipeline: {prediction_long}
id: {site_id}
name: "{name}"
bbox_wgs84: [{bbox}]
layout: nested

classification:
  low_pct: 70
  high_pct: 90
  strong_score_min: 3

temporal:
  start: "2010-01-01"
  end: "2023-12-31"
"""

PREDICTION_LONG = {
    "null": "NO strong-zone enrichment over the site-specific null (deposit style produces no "
            "intrusion-related surface alteration, and/or terrain suppresses the spectral signature).",
    "positive": "strong-zone enrichment ABOVE the site-specific null -- this is an independent "
                "replication test of the arid skarn/carbonate-replacement success regime outside SE Arizona.",
    "positive-ish": "modest strong-zone enrichment expected (exposed alteration present); "
                    "tests whether the category's aggregate below-chance / neutral result is driven "
                    "by buried- or eroded-cap sites in the original sample.",
    "null-to-weak": "at most weak enrichment; alkaline/carbonatite alteration footprints are small "
                    "relative to a scene and only patchily TIR-expressive.",
}


def _wrap(text: str, width: int = 92, prefix: str = "# ") -> str:
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        lines.append(cur)
    return ("\n" + prefix).join(lines)


def main() -> None:
    existing_index = INDEX.read_text().splitlines()
    index_ids = {ln.split("- ", 1)[1].strip() for ln in existing_index if ln.strip().startswith("- ")}

    with CATEGORIES_CSV.open() as f:
        cat_rows = list(csv.reader(f))
    cat_ids = {r[0] for r in cat_rows[1:] if r}

    new_index_lines, new_cat_rows, written, skipped = [], [], [], []

    for site_id, name, bbox, category, confidence, prediction, rationale in SITES:
        yaml_path = SITES_DIR / f"{site_id}.yaml"
        if yaml_path.exists():
            skipped.append(site_id)
        else:
            yaml_path.write_text(YAML_TEMPLATE.format(
                name=name, category=category,
                rationale_wrapped=_wrap(rationale),
                prediction_long=_wrap(PREDICTION_LONG[prediction]),
                site_id=site_id,
                bbox=", ".join(f"{c:g}" for c in bbox),
            ))
            written.append(site_id)

        if site_id not in index_ids:
            new_index_lines.append(f"  - {site_id}")

        if site_id not in cat_ids:
            basis = (f"[EXPANSION +58] {rationale} "
                     f"Pre-registered expectation: {prediction.upper()} -- {PREDICTION_LONG[prediction]}")
            new_cat_rows.append([site_id, category, confidence, basis])

    if new_index_lines:
        with INDEX.open("a") as f:
            if not INDEX.read_text().endswith("\n"):
                f.write("\n")
            f.write("\n".join(new_index_lines) + "\n")

    if new_cat_rows:
        with CATEGORIES_CSV.open("a", newline="") as f:
            csv.writer(f).writerows(new_cat_rows)

    print(f"YAMLs written : {len(written)}")
    print(f"YAMLs skipped (already existed): {len(skipped)}  {skipped}")
    print(f"index.yaml appended: {len(new_index_lines)}")
    print(f"deposit_system_categories.csv appended: {len(new_cat_rows)}")
    print(f"\nTotal sites in index now: {len(index_ids) + len(new_index_lines)}")


if __name__ == "__main__":
    main()
