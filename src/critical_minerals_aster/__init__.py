"""Critical minerals ASTER alteration pipeline (library package)."""

from critical_minerals_aster.classification import (
    classify_percentiles,
    combined_score,
    vectorize_strong_zones,
)
from critical_minerals_aster.config import (
    ClassificationParams,
    SiteConfig,
    StructureLayer,
    list_site_ids,
    load_site_by_id,
    load_site_config,
    search_bbox,
)
from critical_minerals_aster.metrics import (
    compute_site_summary,
    read_mrds_national,
    simplify_commodity,
)
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    is_critical_mineral,
    mrds_to_points_gdf,
    reclassify_mrds_earth_mri,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.paths import SitePaths, site_paths_for
from critical_minerals_aster.pipeline import (
    run_batch,
    run_site,
    save_band_ratio_figure,
    save_composite_figure,
)
from critical_minerals_aster.spectral import (
    alteration_ratios,
    band_ratio,
    clip_bands_to_bbox,
    extract_granule_id,
    load_tir_band,
    load_tir_bands_10_14,
    score_granule,
    select_granule,
)
from critical_minerals_aster.structure import (
    annotate_deposits_with_structure,
    load_structure_layers,
    nearest_structure_distance_m,
)
from critical_minerals_aster.synthesis import load_site_summaries, write_national_summary

__all__ = [
    "ClassificationParams",
    "SiteConfig",
    "SitePaths",
    "StructureLayer",
    "alteration_ratios",
    "annotate_deposits_with_structure",
    "band_ratio",
    "clip_bands_to_bbox",
    "classify_percentiles",
    "combined_score",
    "compute_site_summary",
    "extract_granule_id",
    "filter_mrds_bbox",
    "is_critical_mineral",
    "list_site_ids",
    "load_site_by_id",
    "load_site_config",
    "load_site_summaries",
    "load_structure_layers",
    "load_tir_band",
    "load_tir_bands_10_14",
    "mrds_to_points_gdf",
    "reclassify_mrds_earth_mri",
    "nearest_structure_distance_m",
    "read_mrds_national",
    "run_batch",
    "run_site",
    "save_band_ratio_figure",
    "save_composite_figure",
    "score_granule",
    "search_bbox",
    "select_granule",
    "simplify_commodity",
    "site_paths_for",
    "spatial_join_deposits_zones",
    "vectorize_strong_zones",
    "write_national_summary",
]
