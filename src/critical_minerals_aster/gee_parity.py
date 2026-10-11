"""Local rasterio vs Google Earth Engine parity for the ASTER TIR workflow.

This is a reproducibility / comparison module, not a rewrite of the 109-site
pipeline.  For a small representative panel it runs the same six steps on
both sides:

    1. footprint / scene selection
    2. TIR B10–B14 retrieval
    3. silica / carbonate / mafic ratio calculation
    4. clip to ``site.bbox_wgs84``  (percentiles *after* the clip)
    5. 70th / 90th percentile thresholds
    6. 3-class maps and combined score ≥ 3

and then quantifies agreement.

Why this exists
---------------
The local pipeline already found that skipping the bbox clip moves
percentile thresholds enough to reverse site-level conclusions.  A
cloud-native port has more places to silently diverge: GEE hosts
``ASTER/AST_L1T_003`` while EarthData serves LP DAAC v004; ``reduceRegion``
defaults to EPSG:4326 (resampling the native UTM 90 m grid); the
percentile reducer is histogram-based; TIR unit-conversion coefficients
differ by band, so DN ratios and radiance ratios are not the same map.

The scientific question is therefore implementation-level:

    Does a cloud-native EO implementation reproduce the local raster
    pipeline, and which of those choices move the classified map?

Live Earth Engine is optional.  Tests and the offline benchmark never
import ``ee``.  The GEE recipe (Python graph + Code Editor JavaScript)
is always generated.  When ``earthengine-api`` is installed and
authenticated, ``run_gee_workflow`` pulls cloud-native rasters and the
same agreement metrics apply.

Resume-facing claim this module is built to support
---------------------------------------------------
Ported representative ASTER workflows to Google Earth Engine and
benchmarked cloud-native outputs against a local rasterio pipeline,
testing agreement across scene selection, masking, reprojection and
percentile-based classification.
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Literal, Mapping

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import reproject, transform_bounds

from critical_minerals_aster.classification import (
    classify_percentiles,
    combined_score,
)
from critical_minerals_aster.config import BBox, ClassificationParams, SiteConfig
from critical_minerals_aster.paths import SitePaths
from critical_minerals_aster.spectral import (
    alteration_ratios,
    clip_bands_to_bbox,
    load_ratio_mosaic,
    load_tir_bands_10_14,
)

# GEE hosts the v003 L1T collection.  The local pipeline downloads LP DAAC
# v004 granules via earthaccess.  Same instrument, different archive cut.
GEE_ASTER_L1T_COLLECTION = "ASTER/AST_L1T_003"
TIR_BANDS = ("B10", "B11", "B12", "B13", "B14")

# ASTER L1T TIR unit conversion coefficients (W m^-2 sr^-1 µm^-1 / DN).
# Radiance = (DN − 1) × UCC.  Coefficients differ by band, so a DN ratio
# is not a radiance ratio.  Local pipeline uses raw DN (0 = nodata).
TIR_UCC: dict[str, float] = {
    "B10": 0.006822,
    "B11": 0.006780,
    "B12": 0.006590,
    "B13": 0.005693,
    "B14": 0.005225,
}

# Default panel: 8 sites that already have local TIR, spanning the
# regimes the paper actually talks about.  Not a random convenience sample.
DEFAULT_PARITY_SITES: tuple[str, ...] = (
    "goldfield",  # HS epithermal; pinned granule; SWIR pilot
    "bisbee",  # robust arid skarn
    "christmas",  # robust skarn replicate
    "tombstone",  # robust CRD analog
    "mcdermitt",  # original / first site
    "eureka",  # FDR then clustering failure
    "mountain_pass",  # reversed by the site-specific null
    "marysvale",  # HS + alunite; SWIR complementary
)

SITE_ROLES: dict[str, str] = {
    "goldfield": "HS epithermal",
    "bisbee": "robust arid skarn",
    "christmas": "robust skarn",
    "tombstone": "robust CRD",
    "mcdermitt": "original site",
    "eureka": "clustering failure",
    "mountain_pass": "null-model reversal",
    "marysvale": "HS + alunite",
}

_GRANULE_ACQ_RE = re.compile(
    r"AST_L1T_(\d{3})(\d{2})(\d{2})(\d{4})(\d{2})(\d{2})(\d{2})"
)

ProbeName = Literal[
    "clip_order",
    "reproject_4326",
    "histogram_percentile",
    "radiance_ucc",
    "gee_live",
]


# ---------------------------------------------------------------------------
# Granule ID ↔ GEE image ID
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AsterAcquisition:
    """Acquisition clock parsed from an LP DAAC ``AST_L1T_*`` granule id."""

    version: str
    acquired: datetime
    granule_id: str

    @property
    def gee_image_id(self) -> str:
        stamp = self.acquired.strftime("%Y%m%dT%H%M%S")
        return f"{GEE_ASTER_L1T_COLLECTION}/{stamp}"


def parse_aster_l1t_granule_id(granule_id: str) -> AsterAcquisition | None:
    """Parse ``AST_L1T_00402072010054840_…`` → 2010-02-07 05:48:40 UTC.

    Local ids encode ``VVV + MMDDYYYY + HHMMSS``.  Returns None for mosaic
    ids (``goldfield_mosaic``) or anything that is not an L1T granule.
    """
    match = _GRANULE_ACQ_RE.search(granule_id)
    if match is None:
        return None
    version, month, day, year, hour, minute, second = match.groups()
    try:
        acquired = datetime(
            int(year),
            int(month),
            int(day),
            int(hour),
            int(minute),
            int(second),
            tzinfo=timezone.utc,
        )
    except ValueError:
        return None
    return AsterAcquisition(
        version=version, acquired=acquired, granule_id=granule_id
    )


def gee_image_id_from_granule(granule_id: str) -> str | None:
    parsed = parse_aster_l1t_granule_id(granule_id)
    return None if parsed is None else parsed.gee_image_id


# ---------------------------------------------------------------------------
# Shared recipe (the thing both sides are supposed to execute)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GeeRecipe:
    """Implementation-agnostic description of one site's TIR workflow.

    ``crs=None`` means "use the image's native projection" — the local
    pipeline's behaviour.  Omitting CRS in GEE ``reduceRegion`` silently
    reprojects to EPSG:4326; that is a discrepancy source, not a default
    we want to copy.
    """

    site_id: str
    bbox_wgs84: BBox
    start: str
    end: str
    low_pct: float = 70.0
    high_pct: float = 90.0
    strong_score_min: int = 3
    collection: str = GEE_ASTER_L1T_COLLECTION
    bands: tuple[str, ...] = TIR_BANDS
    scale_m: float = 90.0
    crs: str | None = None
    apply_ucc: bool = False
    pin_image_id: str | None = None
    mosaic: bool = False
    coverage_min: float = 0.30

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["bbox_wgs84"] = list(self.bbox_wgs84)
        d["bands"] = list(self.bands)
        return d


def recipe_for_site(site: SiteConfig, *, mosaic: bool = False) -> GeeRecipe:
    """Build the GEE recipe that mirrors this site's YAML + classification."""
    cp = site.classification or ClassificationParams()
    pin = gee_image_id_from_granule(site.granule_id) if site.granule_id else None
    return GeeRecipe(
        site_id=site.id,
        bbox_wgs84=site.bbox_wgs84,
        start=site.temporal_start,
        end=site.temporal_end,
        low_pct=cp.low_pct,
        high_pct=cp.high_pct,
        strong_score_min=cp.strong_score_min,
        pin_image_id=pin,
        mosaic=mosaic,
    )


def gee_recipe_javascript(recipe: GeeRecipe) -> str:
    """Earth Engine Code Editor JavaScript for one site.

    Generated so the port is inspectable without ``earthengine-api``.
    Comments flag the discrepancy sources the benchmark measures.
    """
    lon0, lat0, lon1, lat1 = recipe.bbox_wgs84
    bbox = f"[{lon0}, {lat0}, {lon1}, {lat1}]"
    pin_block = (
        f"var image = ee.Image('{recipe.pin_image_id}');\n"
        if recipe.pin_image_id
        else (
            "var image = col\n"
            "  .map(function(img) {\n"
            "    var inter = img.geometry().intersection(region, 1);\n"
            "    return img.set('parity_coverage', inter.area(1).divide(region.area(1)));\n"
            "  })\n"
            f"  .filter(ee.Filter.gt('parity_coverage', {recipe.coverage_min}))\n"
            + (
                "  .median();\n"
                if recipe.mosaic
                else "  .sort('parity_coverage', false).first();\n"
            )
        )
    )
    crs_line = (
        f"  crs: '{recipe.crs}',\n"
        if recipe.crs
        else "  // crs omitted → GEE defaults to EPSG:4326 (resampling discrepancy)\n"
    )
    ucc_block = ""
    if recipe.apply_ucc:
        terms = ",\n".join(
            f"    {b}: image.select('{b}').subtract(1).multiply({TIR_UCC[b]})"
            for b in recipe.bands
            if b in TIR_UCC
        )
        ucc_block = (
            "// Radiance = (DN - 1) * UCC.  Off by default: local pipeline uses DN.\n"
            "image = ee.Image.cat({\n"
            f"{terms}\n"
            "});\n"
        )
    return f"""// ASTER TIR parity recipe — {recipe.site_id}
// Collection {recipe.collection} (GEE v003). Local pipeline uses LP DAAC v004.
// Clip the study bbox *before* percentiles. Skipping that clip is a known
// local-pipeline bug class; the GEE analogue is reduceRegion without region.
var region = ee.Geometry.Rectangle({bbox}, 'EPSG:4326', false);
var col = ee.ImageCollection('{recipe.collection}')
  .filterBounds(region)
  .filterDate('{recipe.start}', '{recipe.end}')
  .select({list(recipe.bands)});
{pin_block}image = image.select({list(recipe.bands)}).clip(region);
{ucc_block}image = image.updateMask(image.select('B10').gt(0));
var silica    = image.select('B13').divide(image.select('B14')).rename('silica');
var carbonate = image.select('B13').divide(image.select('B12')).rename('carbonate');
var mafic     = image.select('B12').divide(image.select('B13')).rename('mafic');
var ratios = silica.addBands(carbonate).addBands(mafic);
var pct = ratios.reduceRegion({{
  reducer: ee.Reducer.percentile([{recipe.low_pct}, {recipe.high_pct}]),
  geometry: region,
  scale: {recipe.scale_m},
{crs_line}  maxPixels: 1e9
}});
"""


# ---------------------------------------------------------------------------
# Agreement metrics (work on any pair of aligned arrays)
# ---------------------------------------------------------------------------


def cohen_kappa(a: np.ndarray, b: np.ndarray, labels: Iterable[int] = (0, 1, 2)) -> float:
    """Cohen's κ on two integer class maps of the same shape."""
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    aa = np.asarray(a).ravel()
    bb = np.asarray(b).ravel()
    n = aa.size
    if n == 0:
        return float("nan")
    po = float(np.mean(aa == bb))
    pe = 0.0
    for lab in labels:
        pe += float(np.mean(aa == lab) * np.mean(bb == lab))
    if pe >= 1.0:
        return 1.0 if po >= 1.0 else 0.0
    return float((po - pe) / (1.0 - pe))


def jaccard_binary(mask_a: np.ndarray, mask_b: np.ndarray) -> float:
    """Jaccard index of two boolean masks.  Both empty → 1.0 (agree)."""
    if mask_a.shape != mask_b.shape:
        raise ValueError(f"shape mismatch: {mask_a.shape} vs {mask_b.shape}")
    inter = int(np.logical_and(mask_a, mask_b).sum())
    union = int(np.logical_or(mask_a, mask_b).sum())
    if union == 0:
        return 1.0
    return float(inter / union)


def overall_accuracy(a: np.ndarray, b: np.ndarray) -> float:
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    n = a.size
    if n == 0:
        return float("nan")
    return float(np.mean(np.asarray(a) == np.asarray(b)))


@dataclass
class AgreementStats:
    """Pixel-wise and threshold-wise agreement between two TIR products."""

    n_overlap: int
    overall_accuracy: float
    kappa: float
    jaccard_strong: float
    score_mae: float
    strong_frac_ref: float
    strong_frac_other: float
    strong_frac_delta: float
    p70_delta: dict[str, float] = field(default_factory=dict)
    p90_delta: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "n_overlap": self.n_overlap,
            "overall_accuracy": self.overall_accuracy,
            "kappa": self.kappa,
            "jaccard_strong": self.jaccard_strong,
            "score_mae": self.score_mae,
            "strong_frac_ref": self.strong_frac_ref,
            "strong_frac_other": self.strong_frac_other,
            "strong_frac_delta": self.strong_frac_delta,
        }
        for name, val in self.p70_delta.items():
            d[f"p70_delta_{name}"] = val
        for name, val in self.p90_delta.items():
            d[f"p90_delta_{name}"] = val
        return d


def compare_classified(
    ref_score: np.ndarray,
    other_score: np.ndarray,
    valid: np.ndarray,
    *,
    min_score: int = 3,
    ref_thresholds: Mapping[str, tuple[float, float]] | None = None,
    other_thresholds: Mapping[str, tuple[float, float]] | None = None,
) -> AgreementStats:
    """Agreement on the overlap of two combined-score maps.

    ``valid`` is the pixel-wise AND of both sides' finite masks, already
    aligned to the same grid.  Strong-zone Jaccard uses score ≥ ``min_score``.
    """
    if ref_score.shape != other_score.shape or ref_score.shape != valid.shape:
        raise ValueError("ref, other, and valid must share a shape")
    mask = np.asarray(valid, dtype=bool)
    n = int(mask.sum())
    if n == 0:
        return AgreementStats(
            n_overlap=0,
            overall_accuracy=float("nan"),
            kappa=float("nan"),
            jaccard_strong=float("nan"),
            score_mae=float("nan"),
            strong_frac_ref=float("nan"),
            strong_frac_other=float("nan"),
            strong_frac_delta=float("nan"),
        )
    r = np.asarray(ref_score)[mask]
    o = np.asarray(other_score)[mask]
    strong_r = r >= min_score
    strong_o = o >= min_score
    p70_delta: dict[str, float] = {}
    p90_delta: dict[str, float] = {}
    if ref_thresholds and other_thresholds:
        for key in ref_thresholds:
            if key not in other_thresholds:
                continue
            p70_delta[key] = float(other_thresholds[key][0] - ref_thresholds[key][0])
            p90_delta[key] = float(other_thresholds[key][1] - ref_thresholds[key][1])
    return AgreementStats(
        n_overlap=n,
        overall_accuracy=overall_accuracy(r, o),
        kappa=cohen_kappa(r, o, labels=range(7)),
        jaccard_strong=jaccard_binary(strong_r, strong_o),
        score_mae=float(np.mean(np.abs(r.astype(float) - o.astype(float)))),
        strong_frac_ref=float(np.mean(strong_r)),
        strong_frac_other=float(np.mean(strong_o)),
        strong_frac_delta=float(np.mean(strong_o) - np.mean(strong_r)),
        p70_delta=p70_delta,
        p90_delta=p90_delta,
    )


# ---------------------------------------------------------------------------
# Local workflow (rasterio path — same math as pipeline.run_classification)
# ---------------------------------------------------------------------------


@dataclass
class TIRWorkflowProduct:
    """One side of a parity comparison (local, GEE, or a discrepancy probe)."""

    site_id: str
    source: str
    granule_id: str | None
    gee_image_id: str | None
    crs: str
    silica: np.ndarray
    carbonate: np.ndarray
    mafic: np.ndarray
    silica_cls: np.ndarray
    carbonate_cls: np.ndarray
    mafic_cls: np.ndarray
    combined: np.ndarray
    valid: np.ndarray
    thresholds: dict[str, tuple[float, float]]
    transform: rasterio.Affine | None = None
    raw_bands: dict[str, np.ndarray] | None = None
    unclipped_ratios: dict[str, np.ndarray] | None = None
    notes: str = ""

    @property
    def strong_frac(self) -> float:
        n = int(self.valid.sum())
        if n == 0:
            return float("nan")
        return float(np.mean(self.combined[self.valid] >= 3))

    @property
    def n_valid(self) -> int:
        return int(self.valid.sum())


def _classify_ratios(
    silica: np.ndarray,
    carbonate: np.ndarray,
    mafic: np.ndarray,
    low_pct: float,
    high_pct: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, tuple[float, float]]]:
    silica_cls, s_lo, s_hi = classify_percentiles(silica, low_pct, high_pct)
    carb_cls, c_lo, c_hi = classify_percentiles(carbonate, low_pct, high_pct)
    mafic_cls, m_lo, m_hi = classify_percentiles(mafic, low_pct, high_pct)
    combined = combined_score(silica_cls, carb_cls, mafic_cls)
    thresholds = {
        "silica": (s_lo, s_hi),
        "carbonate": (c_lo, c_hi),
        "mafic": (m_lo, m_hi),
    }
    return silica_cls, carb_cls, mafic_cls, combined, thresholds


def _raw_granule_prefix(aster_dir: Path, site: SiteConfig) -> str | None:
    if site.granule_id:
        if (aster_dir / f"{site.granule_id}_TIR_B10.tif").is_file():
            return site.granule_id
    candidates = sorted(
        p for p in aster_dir.glob("AST_L1T_*_TIR_B10.tif") if "mosaic" not in p.name
    )
    if candidates:
        return candidates[0].name.replace("_TIR_B10.tif", "")
    return None


def run_local_workflow(site: SiteConfig, paths: SitePaths) -> TIRWorkflowProduct:
    """Execute the local TIR recipe: load → clip → ratios → percentiles → score.

    Recomputes classification from ratios (or raw bands) so thresholds are
    available for comparison.  Does not call MRDS or write figures.
    """
    from critical_minerals_aster.pipeline import resolve_granule_id

    cp = site.classification or ClassificationParams()
    aster_dir = paths.aster_dir
    granule_id = resolve_granule_id(site, paths)
    raw_prefix = _raw_granule_prefix(aster_dir, site)

    raw_bands: dict[str, np.ndarray] | None = None
    unclipped: dict[str, np.ndarray] | None = None
    notes: list[str] = []

    if raw_prefix is not None:
        b10, _, b12, b13, b14, _, transform, crs = load_tir_bands_10_14(
            aster_dir, raw_prefix
        )
        unclipped = {}
        u_silica, u_carb, u_mafic = alteration_ratios(b12, b13, b14)
        unclipped["silica"] = u_silica
        unclipped["carbonate"] = u_carb
        unclipped["mafic"] = u_mafic
        (b10, b12, b13, b14), transform = clip_bands_to_bbox(
            [b10, b12, b13, b14], transform, crs, site.bbox_wgs84
        )
        silica, carbonate, mafic = alteration_ratios(b12, b13, b14)
        raw_bands = {"B10": b10, "B12": b12, "B13": b13, "B14": b14}
        notes.append(f"raw granule {raw_prefix}")
        granule_id = raw_prefix
    else:
        silica, carbonate, mafic, b10, transform, crs = load_ratio_mosaic(
            aster_dir, granule_id
        )
        unclipped = {
            "silica": silica.copy(),
            "carbonate": carbonate.copy(),
            "mafic": mafic.copy(),
        }
        (b10, silica, carbonate, mafic), transform = clip_bands_to_bbox(
            [b10, silica, carbonate, mafic], transform, crs, site.bbox_wgs84
        )
        notes.append(f"ratio mosaic {granule_id}")

    silica_cls, carb_cls, mafic_cls, combined, thresholds = _classify_ratios(
        silica, carbonate, mafic, cp.low_pct, cp.high_pct
    )
    valid = np.isfinite(silica) & np.isfinite(carbonate) & np.isfinite(mafic)
    return TIRWorkflowProduct(
        site_id=site.id,
        source="local",
        granule_id=granule_id,
        gee_image_id=gee_image_id_from_granule(granule_id),
        crs=str(crs) if crs is not None else "",
        silica=silica,
        carbonate=carbonate,
        mafic=mafic,
        silica_cls=silica_cls,
        carbonate_cls=carb_cls,
        mafic_cls=mafic_cls,
        combined=combined,
        valid=valid,
        thresholds=thresholds,
        transform=transform,
        raw_bands=raw_bands,
        unclipped_ratios=unclipped,
        notes="; ".join(notes),
    )


# ---------------------------------------------------------------------------
# Discrepancy probes (offline stand-ins for GEE implementation choices)
# ---------------------------------------------------------------------------


def histogram_percentile(
    values: np.ndarray,
    pct: float,
    n_bins: int = 10_000,
) -> float:
    """Histogram-bin percentile, analogue of GEE ``Reducer.percentile``.

    ``np.nanpercentile`` is exact on the sample.  GEE's reducer bins first.
    On smooth TIR ratios the two stay close; on peaked distributions the
    bin edges move the 70th / 90th enough to flip edge pixels.
    """
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan")
    vmin = float(finite.min())
    vmax = float(finite.max())
    if vmin == vmax:
        return vmin
    counts, edges = np.histogram(finite, bins=n_bins, range=(vmin, vmax))
    total = int(counts.sum())
    if total == 0:
        return vmin
    cdf = np.cumsum(counts) / total
    return float(np.interp(pct / 100.0, cdf, edges[1:]))


def classify_histogram_percentiles(
    ratio: np.ndarray,
    low_pct: float = 70.0,
    high_pct: float = 90.0,
    n_bins: int = 10_000,
) -> tuple[np.ndarray, float, float]:
    """Same contract as ``classify_percentiles``, histogram thresholds."""
    low = histogram_percentile(ratio, low_pct, n_bins=n_bins)
    high = histogram_percentile(ratio, high_pct, n_bins=n_bins)
    classes = np.zeros_like(ratio, dtype=np.uint8)
    classes[ratio >= low] = 1
    classes[ratio >= high] = 2
    classes[np.isnan(ratio)] = 0
    return classes, low, high


def dn_to_radiance(dn: np.ndarray, band: str) -> np.ndarray:
    """Apply the ASTER L1T TIR UCC.  DN 0 stays nodata."""
    if band not in TIR_UCC:
        raise KeyError(f"no UCC for {band}")
    out = (dn.astype(float) - 1.0) * TIR_UCC[band]
    out[~np.isfinite(dn) | (dn == 0)] = np.nan
    return out


def _product_from_ratios(
    site_id: str,
    source: str,
    silica: np.ndarray,
    carbonate: np.ndarray,
    mafic: np.ndarray,
    *,
    low_pct: float,
    high_pct: float,
    crs: str,
    granule_id: str | None,
    gee_image_id: str | None,
    transform: rasterio.Affine | None = None,
    notes: str = "",
    histogram: bool = False,
) -> TIRWorkflowProduct:
    if histogram:
        silica_cls, s_lo, s_hi = classify_histogram_percentiles(silica, low_pct, high_pct)
        carb_cls, c_lo, c_hi = classify_histogram_percentiles(carbonate, low_pct, high_pct)
        mafic_cls, m_lo, m_hi = classify_histogram_percentiles(mafic, low_pct, high_pct)
        combined = combined_score(silica_cls, carb_cls, mafic_cls)
        thresholds = {
            "silica": (s_lo, s_hi),
            "carbonate": (c_lo, c_hi),
            "mafic": (m_lo, m_hi),
        }
    else:
        silica_cls, carb_cls, mafic_cls, combined, thresholds = _classify_ratios(
            silica, carbonate, mafic, low_pct, high_pct
        )
    valid = np.isfinite(silica) & np.isfinite(carbonate) & np.isfinite(mafic)
    return TIRWorkflowProduct(
        site_id=site_id,
        source=source,
        granule_id=granule_id,
        gee_image_id=gee_image_id,
        crs=crs,
        silica=silica,
        carbonate=carbonate,
        mafic=mafic,
        silica_cls=silica_cls,
        carbonate_cls=carb_cls,
        mafic_cls=mafic_cls,
        combined=combined,
        valid=valid,
        thresholds=thresholds,
        transform=transform,
        notes=notes,
    )


def probe_clip_order(
    local: TIRWorkflowProduct,
    site: SiteConfig,
) -> TIRWorkflowProduct | None:
    """Apply *unclipped* percentiles to the clipped grid.

    This is the diagnostic-script bug the project already found: compute
    70/90 over the whole scene, then classify the study box.  The GEE
    analogue is ``reduceRegion`` on the unclipped image footprint.
    """
    if local.unclipped_ratios is None:
        return None
    cp = site.classification or ClassificationParams()
    un = local.unclipped_ratios
    # Thresholds from the unclipped scene…
    _, s_lo, s_hi = classify_percentiles(un["silica"], cp.low_pct, cp.high_pct)
    _, c_lo, c_hi = classify_percentiles(un["carbonate"], cp.low_pct, cp.high_pct)
    _, m_lo, m_hi = classify_percentiles(un["mafic"], cp.low_pct, cp.high_pct)
    # …applied to the clipped pixels (same grid as ``local``).
    silica_cls = np.zeros_like(local.silica, dtype=np.uint8)
    silica_cls[local.silica >= s_lo] = 1
    silica_cls[local.silica >= s_hi] = 2
    silica_cls[np.isnan(local.silica)] = 0
    carb_cls = np.zeros_like(local.carbonate, dtype=np.uint8)
    carb_cls[local.carbonate >= c_lo] = 1
    carb_cls[local.carbonate >= c_hi] = 2
    carb_cls[np.isnan(local.carbonate)] = 0
    mafic_cls = np.zeros_like(local.mafic, dtype=np.uint8)
    mafic_cls[local.mafic >= m_lo] = 1
    mafic_cls[local.mafic >= m_hi] = 2
    mafic_cls[np.isnan(local.mafic)] = 0
    combined = combined_score(silica_cls, carb_cls, mafic_cls)
    return TIRWorkflowProduct(
        site_id=local.site_id,
        source="probe:clip_order",
        granule_id=local.granule_id,
        gee_image_id=local.gee_image_id,
        crs=local.crs,
        silica=local.silica,
        carbonate=local.carbonate,
        mafic=local.mafic,
        silica_cls=silica_cls,
        carbonate_cls=carb_cls,
        mafic_cls=mafic_cls,
        combined=combined,
        valid=local.valid,
        thresholds={
            "silica": (s_lo, s_hi),
            "carbonate": (c_lo, c_hi),
            "mafic": (m_lo, m_hi),
        },
        transform=local.transform,
        notes="percentiles from unclipped scene, applied to bbox-clipped pixels",
    )


def probe_histogram_percentile(
    local: TIRWorkflowProduct,
    site: SiteConfig,
) -> TIRWorkflowProduct:
    """Replace ``nanpercentile`` with a GEE-like histogram reducer."""
    cp = site.classification or ClassificationParams()
    return _product_from_ratios(
        local.site_id,
        "probe:histogram_percentile",
        local.silica,
        local.carbonate,
        local.mafic,
        low_pct=cp.low_pct,
        high_pct=cp.high_pct,
        crs=local.crs,
        granule_id=local.granule_id,
        gee_image_id=local.gee_image_id,
        transform=local.transform,
        notes="histogram-bin percentiles (GEE Reducer.percentile analogue)",
        histogram=True,
    )


def probe_radiance_ucc(
    local: TIRWorkflowProduct,
    site: SiteConfig,
) -> TIRWorkflowProduct | None:
    """Recompute ratios from radiance instead of DN.

    GEE catalog docs describe L1T TIR as DN; some EE recipes convert with
    the published UCCs first.  Because UCC_B13 ≠ UCC_B14, silica = B13/B14
    changes, and so can the classified map.
    """
    if local.raw_bands is None:
        return None
    required = {"B12", "B13", "B14"}
    if not required.issubset(local.raw_bands):
        return None
    cp = site.classification or ClassificationParams()
    b12 = dn_to_radiance(local.raw_bands["B12"], "B12")
    b13 = dn_to_radiance(local.raw_bands["B13"], "B13")
    b14 = dn_to_radiance(local.raw_bands["B14"], "B14")
    silica, carbonate, mafic = alteration_ratios(b12, b13, b14)
    return _product_from_ratios(
        local.site_id,
        "probe:radiance_ucc",
        silica,
        carbonate,
        mafic,
        low_pct=cp.low_pct,
        high_pct=cp.high_pct,
        crs=local.crs,
        granule_id=local.granule_id,
        gee_image_id=local.gee_image_id,
        transform=local.transform,
        notes="ratios from (DN-1)*UCC radiance; local pipeline uses DN",
    )


def _reproject_array(
    arr: np.ndarray,
    src_transform: rasterio.Affine,
    src_crs: CRS,
    dst_crs: CRS,
    *,
    resampling: Resampling = Resampling.nearest,
    dst_resolution: float | None = None,
) -> tuple[np.ndarray, rasterio.Affine]:
    """Warp one band to ``dst_crs``.  NaNs stay NaN (nodata = nan)."""
    src = np.asarray(arr, dtype=np.float64)
    left, bottom, right, top = transform_bounds(
        src_crs,
        dst_crs,
        src_transform.c,
        src_transform.f + src.shape[0] * src_transform.e,
        src_transform.c + src.shape[1] * src_transform.a,
        src_transform.f,
    )
    if dst_resolution is None:
        # ~90 m in degrees at mid-latitude of the source grid.
        mid_lat = (bottom + top) / 2.0
        metres = 90.0
        dst_resolution = metres / (111_320.0 * max(0.2, math.cos(math.radians(mid_lat))))
    width = max(1, int(math.ceil((right - left) / dst_resolution)))
    height = max(1, int(math.ceil((top - bottom) / dst_resolution)))
    dst_transform = rasterio.Affine(dst_resolution, 0.0, left, 0.0, -dst_resolution, top)
    dest = np.full((height, width), np.nan, dtype=np.float64)
    reproject(
        source=src,
        destination=dest,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        src_nodata=np.nan,
        dst_nodata=np.nan,
        resampling=resampling,
    )
    return dest, dst_transform


def probe_reproject_4326(
    local: TIRWorkflowProduct,
    site: SiteConfig,
) -> TIRWorkflowProduct | None:
    """Reproject ratios to EPSG:4326, classify, warp classes back.

    GEE ``reduceRegion`` without an explicit ``crs`` samples in WGS84.
    That is a different grid than the local UTM 90 m raster, even at a
    nominal 90 m scale.  Classes are warped back with nearest-neighbour
    so Jaccard is on the local grid.
    """
    if local.transform is None or not local.crs:
        return None
    try:
        src_crs = CRS.from_user_input(local.crs)
    except Exception:
        return None
    if src_crs.to_epsg() == 4326:
        return None
    dst_crs = CRS.from_epsg(4326)
    cp = site.classification or ClassificationParams()

    sil_w, dst_transform = _reproject_array(
        local.silica, local.transform, src_crs, dst_crs, resampling=Resampling.bilinear
    )
    carb_w, _ = _reproject_array(
        local.carbonate, local.transform, src_crs, dst_crs, resampling=Resampling.bilinear
    )
    maf_w, _ = _reproject_array(
        local.mafic, local.transform, src_crs, dst_crs, resampling=Resampling.bilinear
    )
    warped = _product_from_ratios(
        local.site_id,
        "probe:reproject_4326",
        sil_w,
        carb_w,
        maf_w,
        low_pct=cp.low_pct,
        high_pct=cp.high_pct,
        crs="EPSG:4326",
        granule_id=local.granule_id,
        gee_image_id=local.gee_image_id,
        transform=dst_transform,
        notes="bilinear warp to EPSG:4326 ~90 m, then classify (GEE default CRS)",
    )
    # Warp the combined score back to the local grid for pixel-wise metrics.
    back = np.full(local.combined.shape, 0, dtype=np.float64)
    reproject(
        source=warped.combined.astype(np.float64),
        destination=back,
        src_transform=dst_transform,
        src_crs=dst_crs,
        dst_transform=local.transform,
        dst_crs=src_crs,
        src_nodata=np.nan,
        dst_nodata=np.nan,
        resampling=Resampling.nearest,
    )
    back[~np.isfinite(back)] = 0
    warped.combined = back.astype(np.uint8)
    warped.valid = local.valid
    warped.silica = local.silica
    warped.carbonate = local.carbonate
    warped.mafic = local.mafic
    warped.transform = local.transform
    warped.crs = local.crs
    return warped


# ---------------------------------------------------------------------------
# Live Earth Engine (optional; never imported by tests)
# ---------------------------------------------------------------------------


def initialize_earth_engine(project: str | None = None) -> Any:
    """``ee.Initialize`` with a clear error if the extra is missing."""
    try:
        import ee
    except ImportError as exc:
        raise RuntimeError(
            "earthengine-api is not installed. "
            "pip install 'critical-minerals-aster[gee]' "
            "(or: pip install earthengine-api)"
        ) from exc
    try:
        if project:
            ee.Initialize(project=project)
        else:
            ee.Initialize()
    except Exception as exc:
        raise RuntimeError(
            "Earth Engine is not authenticated. "
            "Run `earthengine authenticate` and pass a Cloud project "
            "via --project / initialize_earth_engine(project=...)."
        ) from exc
    return ee


def build_gee_image(ee: Any, recipe: GeeRecipe) -> tuple[Any, Any, Any]:
    """Construct (region, tir_image, ratios) without computing anything.

    Mirrors ``gee_recipe_javascript``.  Callers trigger computation with
    ``reduceRegion`` / ``getDownloadURL``.
    """
    region = ee.Geometry.Rectangle(list(recipe.bbox_wgs84), "EPSG:4326", False)
    if recipe.pin_image_id:
        image = ee.Image(recipe.pin_image_id)
    else:
        col = (
            ee.ImageCollection(recipe.collection)
            .filterBounds(region)
            .filterDate(recipe.start, recipe.end)
            .select(list(recipe.bands))
        )

        def _score(img: Any) -> Any:
            inter = img.geometry().intersection(region, 1)
            return img.set("parity_coverage", inter.area(1).divide(region.area(1)))

        col = col.map(_score).filter(ee.Filter.gt("parity_coverage", recipe.coverage_min))
        image = col.median() if recipe.mosaic else col.sort("parity_coverage", False).first()

    image = image.select(list(recipe.bands)).clip(region)
    if recipe.apply_ucc:
        parts = [
            image.select(band).subtract(1).multiply(TIR_UCC[band]).rename(band)
            for band in recipe.bands
            if band in TIR_UCC
        ]
        image = ee.Image.cat(parts)
    image = image.updateMask(image.select("B10").gt(0))
    silica = image.select("B13").divide(image.select("B14")).rename("silica")
    carbonate = image.select("B13").divide(image.select("B12")).rename("carbonate")
    mafic = image.select("B12").divide(image.select("B13")).rename("mafic")
    ratios = silica.addBands(carbonate).addBands(mafic)
    return region, image, ratios


def run_gee_workflow(
    site: SiteConfig,
    recipe: GeeRecipe | None = None,
    *,
    project: str | None = None,
    ee: Any | None = None,
) -> dict[str, Any]:
    """Execute the GEE recipe and return reduceRegion stats + image id.

    Pixel arrays are *not* pulled here (export is slow and size-capped).
    The benchmark script can call ``download_gee_ratios`` separately when
    a GeoTIFF comparison is wanted.  This function is the live half of
    the resume claim: it actually talks to Earth Engine.
    """
    if ee is None:
        ee = initialize_earth_engine(project)
    recipe = recipe or recipe_for_site(site)
    region, image, ratios = build_gee_image(ee, recipe)
    reducer = ee.Reducer.percentile([recipe.low_pct, recipe.high_pct])
    kwargs: dict[str, Any] = {
        "reducer": reducer,
        "geometry": region,
        "scale": recipe.scale_m,
        "maxPixels": 1e9,
        "bestEffort": True,
    }
    if recipe.crs:
        kwargs["crs"] = recipe.crs
    stats = ratios.reduceRegion(**kwargs)
    image_id = image.get("system:id")
    coverage = image.get("parity_coverage")
    info = ee.Dictionary(
        {"percentiles": stats, "image_id": image_id, "coverage": coverage}
    ).getInfo()
    return {
        "site_id": site.id,
        "source": "gee_live",
        "recipe": recipe.to_dict(),
        "gee": info,
        "javascript": gee_recipe_javascript(recipe),
    }


# ---------------------------------------------------------------------------
# Site-level comparison table
# ---------------------------------------------------------------------------


@dataclass
class ParityRow:
    site_id: str
    site_name: str
    role: str
    probe: str
    n_overlap: int
    overall_accuracy: float
    kappa: float
    jaccard_strong: float
    score_mae: float
    strong_frac_local: float
    strong_frac_other: float
    strong_frac_delta: float
    p70_delta_silica: float | None
    p90_delta_silica: float | None
    local_granule_id: str | None
    gee_image_id: str | None
    local_crs: str
    notes: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def parity_row(
    site: SiteConfig,
    local: TIRWorkflowProduct,
    other: TIRWorkflowProduct,
    stats: AgreementStats,
) -> ParityRow:
    return ParityRow(
        site_id=site.id,
        site_name=site.name,
        role=SITE_ROLES.get(site.id, ""),
        probe=other.source,
        n_overlap=stats.n_overlap,
        overall_accuracy=stats.overall_accuracy,
        kappa=stats.kappa,
        jaccard_strong=stats.jaccard_strong,
        score_mae=stats.score_mae,
        strong_frac_local=stats.strong_frac_ref,
        strong_frac_other=stats.strong_frac_other,
        strong_frac_delta=stats.strong_frac_delta,
        p70_delta_silica=stats.p70_delta.get("silica"),
        p90_delta_silica=stats.p90_delta.get("silica"),
        local_granule_id=local.granule_id,
        gee_image_id=local.gee_image_id or other.gee_image_id,
        local_crs=local.crs,
        notes=other.notes,
    )


def run_site_probes(
    site: SiteConfig,
    paths: SitePaths,
) -> tuple[TIRWorkflowProduct, list[ParityRow]]:
    """Local TIR product + every offline discrepancy probe that applies."""
    local = run_local_workflow(site, paths)
    rows: list[ParityRow] = []
    probes: list[TIRWorkflowProduct | None] = [
        probe_clip_order(local, site),
        probe_histogram_percentile(local, site),
        probe_radiance_ucc(local, site),
        probe_reproject_4326(local, site),
    ]
    min_score = (site.classification or ClassificationParams()).strong_score_min
    for other in probes:
        if other is None:
            continue
        stats = compare_classified(
            local.combined,
            other.combined,
            local.valid & other.valid,
            min_score=min_score,
            ref_thresholds=local.thresholds,
            other_thresholds=other.thresholds,
        )
        rows.append(parity_row(site, local, other, stats))
    return local, rows


def summarize_probe_effect(rows: Iterable[ParityRow]) -> dict[str, dict[str, float]]:
    """Mean |Δ strong fraction| and mean Jaccard per probe, across sites."""
    by_probe: dict[str, list[ParityRow]] = {}
    for row in rows:
        by_probe.setdefault(row.probe, []).append(row)
    out: dict[str, dict[str, float]] = {}
    for probe, group in by_probe.items():
        jacc = [r.jaccard_strong for r in group if np.isfinite(r.jaccard_strong)]
        delta = [abs(r.strong_frac_delta) for r in group if np.isfinite(r.strong_frac_delta)]
        kappa = [r.kappa for r in group if np.isfinite(r.kappa)]
        out[probe] = {
            "n_sites": float(len(group)),
            "mean_jaccard_strong": float(np.mean(jacc)) if jacc else float("nan"),
            "mean_abs_strong_frac_delta": float(np.mean(delta)) if delta else float("nan"),
            "mean_kappa": float(np.mean(kappa)) if kappa else float("nan"),
        }
    return out
