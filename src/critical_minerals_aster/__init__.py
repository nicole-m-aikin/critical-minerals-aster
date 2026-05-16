"""Critical minerals ASTER alteration pipeline (library package)."""

from critical_minerals_aster.classification import (
    classify_percentiles,
    combined_score,
    vectorize_strong_zones,
)
from critical_minerals_aster.config import ClassificationParams, SiteConfig, load_site_config
from critical_minerals_aster.mrds import (
    filter_mrds_bbox,
    mrds_to_points_gdf,
    spatial_join_deposits_zones,
)
from critical_minerals_aster.paths import SitePaths, site_paths_for
from critical_minerals_aster.spectral import (
    alteration_ratios,
    band_ratio,
    load_tir_band,
    load_tir_bands_10_14,
)

__all__ = [
    "ClassificationParams",
    "SiteConfig",
    "SitePaths",
    "alteration_ratios",
    "band_ratio",
    "classify_percentiles",
    "combined_score",
    "filter_mrds_bbox",
    "load_site_config",
    "load_tir_band",
    "load_tir_bands_10_14",
    "mrds_to_points_gdf",
    "site_paths_for",
    "spatial_join_deposits_zones",
    "vectorize_strong_zones",
]
