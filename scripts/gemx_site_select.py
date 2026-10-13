"""Site-list selection and pin resolution for GEMx SWIR batch runs.

Pure helpers (no cube downloads).  Used by ``gemx_swir_pilot``,
``gemx_swir_correspondence``, and ``gemx_host_vs_alteration``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

AVIRIS5_COLLECTION = "aviris5_l2a_rfl"
AVIRISC_COLLECTION = "avirisc_l2_rfl"
SWIR_FULL_COLLECTIONS = (AVIRIS5_COLLECTION, AVIRISC_COLLECTION)
COLLECTION_SENSOR = {
    AVIRIS5_COLLECTION: "AVIRIS-5",
    AVIRISC_COLLECTION: "AVIRIS-Classic",
}

# Town / deposit pins (lon, lat).  Real GNIS / published points, not bbox edges.
SITE_PINS = {
    "goldfield": {
        "town": (-117.23563, 37.70854),
        "town_label": "Goldfield",
        "target": (-117.14341, 37.71938),
        "target_label": "district",
    },
    "marysvale": {
        "town": (-112.2302, 38.4494),
        "town_label": "Marysvale",
        "target": (-112.33257, 38.38160),
        "target_label": "Alunite Ridge",
    },
    # GNIS populated place Bodie (1658094) and GNIS cliff Bodie Bluff
    # (219578). Low-sulfidation negative control: alunite-like should
    # stay rare. Not a high-sulfidation analog.
    "bodie": {
        "town": (-119.01209, 38.21214),
        "town_label": "Bodie",
        "town_offset": (-10, -10),
        "target": (-119.00125, 38.21964),
        "target_label": "Bodie Bluff",
        "target_offset": (10, 8),
    },
}

# Al-OH correspondence verdict thresholds (documented in correspondence_verdict).
JACCARD_OVERLAP = 0.15
LIFT_OVERLAP = 1.3
LIFT_COMPLEMENTARY = 1.1
MIN_OVERLAP_CELLS = 100


def bbox_centroid(bbox_wgs84: Sequence[float]) -> tuple[float, float]:
    """(lon, lat) center of a WGS84 west/south/east/north bbox."""
    west, south, east, north = bbox_wgs84
    return ((west + east) / 2.0, (south + north) / 2.0)


def yaml_center_wgs84(site_id: str, sites_dir: Path) -> tuple[float, float] | None:
    """Optional ``center_wgs84`` / ``center`` [lon, lat] from a site YAML."""
    import yaml

    path = Path(sites_dir) / f"{site_id}.yaml"
    if not path.is_file():
        return None
    raw = yaml.safe_load(path.read_text()) or {}
    c = raw.get("center_wgs84") or raw.get("center")
    if isinstance(c, (list, tuple)) and len(c) == 2:
        return float(c[0]), float(c[1])
    return None


def resolve_pins(
    site,
    *,
    pins_table: dict | None = None,
    sites_dir: Path | None = None,
) -> dict:
    """Town/target seeds for granule selection and map annotation.

    Prefer ``SITE_PINS``; else optional YAML ``center`` / ``center_wgs84``;
    else the site bbox centroid for both town and target.
    """
    table = SITE_PINS if pins_table is None else pins_table
    if site.id in table:
        out = dict(table[site.id])
        out.setdefault("pins_source", "site_pins")
        return out

    lon, lat = bbox_centroid(site.bbox_wgs84)
    source = "bbox_centroid"
    if sites_dir is not None:
        center = yaml_center_wgs84(site.id, Path(sites_dir))
        if center is not None:
            lon, lat = center
            source = "yaml_center"

    label = site.name.split(",")[0].strip() if getattr(site, "name", None) else site.id
    return {
        "town": (lon, lat),
        "town_label": label,
        "target": (lon, lat),
        "target_label": "center",
        "pins_source": source,
    }


def primary_fdr_sites(repo_root: Path) -> list[str]:
    """Phase-4 ``sig_fdr_binomial`` site ids (currently 19)."""
    from critical_minerals_aster.significance import significant_site_sets

    fdr, _ = significant_site_sets(repo_root)
    return list(fdr)


def swir_covered_sites(
    repo_root: Path,
    *,
    coverage_csv: Path | None = None,
) -> list[str]:
    """Sites with ``has_swir`` in ``results/gemx_coverage_by_site.csv``."""
    import pandas as pd

    path = coverage_csv or (Path(repo_root) / "results" / "gemx_coverage_by_site.csv")
    df = pd.read_csv(path)
    mask = df["has_swir"].astype(bool) if "has_swir" in df.columns else df["swir_coverage_frac"] > 0
    return sorted(df.loc[mask, "site_id"].astype(str))


def aviris5_coverage_by_site(
    repo_root: Path,
    *,
    long_csv: Path | None = None,
) -> dict[str, float]:
    """``site_id -> coverage_frac`` for AVIRIS-5 (``aviris5_l2a_rfl``)."""
    import pandas as pd

    path = long_csv or (Path(repo_root) / "results" / "gemx_coverage.csv")
    df = pd.read_csv(path)
    sub = df[df["collection"] == AVIRIS5_COLLECTION]
    out: dict[str, float] = {}
    for row in sub.itertuples(index=False):
        frac = float(row.coverage_frac)
        if frac > 0:
            out[str(row.site_id)] = frac
    return out



def select_census_collection(
    site_id: str,
    long_rows,
    *,
    good_coverage: float = 0.75,
) -> tuple[str, float, str] | None:
    """Pick one full-SWIR collection for *site_id*.

    Rule: among collections with coverage >= *good_coverage* (usable/good),
    prefer AVIRIS-5; otherwise take greatest coverage, with AVIRIS-5
    tie-break. Returns ``(collection_key, coverage_frac, sensor)`` or ``None``.
    """
    rows = [
        r for r in long_rows
        if str(r.get("site_id")) == str(site_id)
        and str(r.get("collection")) in SWIR_FULL_COLLECTIONS
        and float(r.get("coverage_frac") or 0.0) > 0.0
    ]
    if not rows:
        return None
    usable = [r for r in rows if float(r["coverage_frac"]) >= good_coverage]
    pool = usable if usable else rows
    # Prefer AVIRIS-5 inside the usable/good pool.
    av5 = [r for r in pool if str(r["collection"]) == AVIRIS5_COLLECTION]
    if av5:
        chosen = max(av5, key=lambda r: float(r["coverage_frac"]))
    else:
        best_frac = max(float(r["coverage_frac"]) for r in pool)
        tied = [r for r in pool if float(r["coverage_frac"]) == best_frac]
        chosen = next(
            (r for r in tied if str(r["collection"]) == AVIRIS5_COLLECTION),
            tied[0],
        )
    key = str(chosen["collection"])
    return key, float(chosen["coverage_frac"]), COLLECTION_SENSOR[key]


def census_collection_by_site(
    repo_root: Path,
    *,
    long_csv: Path | None = None,
) -> dict[str, tuple[str, float, str]]:
    """``site_id -> (collection_key, coverage_frac, sensor)`` for the 109-site census."""
    import pandas as pd

    path = long_csv or (Path(repo_root) / "results" / "gemx_coverage.csv")
    df = pd.read_csv(path)
    rows = df.to_dict(orient="records")
    out: dict[str, tuple[str, float, str]] = {}
    for site_id in sorted({str(r["site_id"]) for r in rows}):
        picked = select_census_collection(site_id, rows)
        if picked is not None:
            out[site_id] = picked
    return out


def swir_coverage_by_site(
    repo_root: Path,
    *,
    long_csv: Path | None = None,
    collection: str | None = None,
) -> dict[str, float]:
    """``site_id -> coverage_frac`` for one collection or the census pick."""
    if collection == AVIRIS5_COLLECTION:
        return aviris5_coverage_by_site(repo_root, long_csv=long_csv)
    picked = census_collection_by_site(repo_root, long_csv=long_csv)
    if collection is None:
        return {sid: frac for sid, (_k, frac, _s) in picked.items()}
    return {
        sid: frac
        for sid, (key, frac, _s) in picked.items()
        if key == collection
    }


def select_sites(
    *,
    repo_root: Path,
    site: str | None = None,
    sites: Sequence[str] | None = None,
    primary_fdr: bool = False,
    all_covered: bool = False,
    default_primary_fdr: bool = True,
) -> list[str]:
    """Resolve a site list from mutually exclusive CLI modes.

    Default (no flags) is ``--primary-fdr`` when ``default_primary_fdr`` is True.
    """
    n_modes = sum(
        [
            site is not None,
            bool(sites),
            primary_fdr,
            all_covered,
        ]
    )
    if n_modes > 1:
        raise ValueError(
            "choose only one of --site, --sites, --primary-fdr, --all-covered"
        )
    if site is not None:
        return [site]
    if sites:
        return list(sites)
    if all_covered:
        return swir_covered_sites(repo_root)
    if primary_fdr or (n_modes == 0 and default_primary_fdr):
        return primary_fdr_sites(repo_root)
    raise ValueError("no site selection mode given")


def correspondence_verdict(
    *,
    status: str = "ok",
    swir_coverage: float | None = None,
    jaccard_aloh: float | None = None,
    lift_aloh: float | None = None,
    n_overlap: int | None = None,
    n_tir_strong: int | None = None,
    min_overlap: int = MIN_OVERLAP_CELLS,
) -> str:
    """Classify TIR strong-zone vs SWIR Al-OH spatial relationship.

    Heuristic (Al-OH is primary):

    - ``no_coverage`` — no AVIRIS-5 / SWIR coverage for the site
    - ``insufficient`` — missing metrics, tiny overlap, or no TIR-strong cells;
      also the mixed band (Jaccard < 0.15 and 1.1 ≤ lift < 1.3) treated as
      inconclusive / insufficient for a firm complementary|overlapping call
    - ``overlapping`` — Jaccard ≥ 0.15 **or** lift ≥ 1.3 (TIR zones coincide
      with clay alteration)
    - ``complementary`` — Jaccard < 0.15 **and** lift < 1.1 (TIR not SWIR-clay)
    - ``error`` — caller already marked a hard failure
    - ``mixed`` — reserved if caller wants the inconclusive band named
      explicitly; this function returns ``insufficient`` for that band so the
      CSV verdict stays in
      {complementary, overlapping, insufficient, no_coverage, error}
    """
    if status == "error":
        return "error"
    if status == "no_coverage" or (swir_coverage is not None and swir_coverage <= 0):
        return "no_coverage"
    if status in ("insufficient", "no_tiffs", "no_tir"):
        return "insufficient"
    if n_overlap is not None and n_overlap < min_overlap:
        return "insufficient"
    if n_tir_strong is not None and n_tir_strong <= 0:
        return "insufficient"
    if jaccard_aloh is None or lift_aloh is None:
        return "insufficient"
    try:
        j = float(jaccard_aloh)
        lift = float(lift_aloh)
    except (TypeError, ValueError):
        return "insufficient"
    if j != j or lift != lift:  # NaN
        return "insufficient"
    if j >= JACCARD_OVERLAP or lift >= LIFT_OVERLAP:
        return "overlapping"
    if j < JACCARD_OVERLAP and lift < LIFT_COMPLEMENTARY:
        return "complementary"
    # Mixed / inconclusive band → insufficient for a firm call
    return "insufficient"


def add_site_select_args(parser) -> None:
    """Attach ``--site`` / ``--sites`` / ``--primary-fdr`` / ``--all-covered``."""
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--site", type=str, default=None, help="Single site id")
    g.add_argument(
        "--sites",
        nargs="+",
        default=None,
        help="Explicit site id list",
    )
    g.add_argument(
        "--primary-fdr",
        action="store_true",
        help="Phase-4 sig_fdr_binomial sites (default for batch)",
    )
    g.add_argument(
        "--all-covered",
        action="store_true",
        help="All sites with has_swir in gemx_coverage_by_site.csv",
    )


def sites_from_namespace(ns, repo_root: Path) -> list[str]:
    """``select_sites`` from an argparse namespace."""
    return select_sites(
        repo_root=repo_root,
        site=getattr(ns, "site", None),
        sites=getattr(ns, "sites", None),
        primary_fdr=bool(getattr(ns, "primary_fdr", False)),
        all_covered=bool(getattr(ns, "all_covered", False)),
        default_primary_fdr=True,
    )


def dedupe_preserve(ids: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for s in ids:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out
