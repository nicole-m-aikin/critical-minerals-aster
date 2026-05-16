"""Filesystem paths for site data and figures."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from critical_minerals_aster.config import SiteConfig


@dataclass(frozen=True)
class SitePaths:
    """Resolved directories for ASTER rasters, vectors, figures, and results."""

    repo_root: Path
    site: SiteConfig

    @property
    def aster_dir(self) -> Path:
        if self.site.layout == "nested":
            return self.repo_root / "data" / "sites" / self.site.id / "aster"
        return self.repo_root / "data" / "aster"

    @property
    def vectors_dir(self) -> Path:
        if self.site.layout == "nested":
            return self.repo_root / "data" / "sites" / self.site.id / "vectors"
        return self.repo_root / "data" / "vectors"

    @property
    def figures_dir(self) -> Path:
        if self.site.layout == "nested":
            return self.repo_root / "figures" / "sites" / self.site.id
        return self.repo_root / "figures"

    @property
    def mrds_dir(self) -> Path:
        return self.repo_root / "data" / "mrds"

    @property
    def mrds_csv(self) -> Path:
        legacy = self.repo_root / "data" / "mrds.csv"
        if legacy.is_file():
            return legacy
        return self.mrds_dir / "mrds.csv"

    @property
    def results_dir(self) -> Path:
        return self.repo_root / "results"

    @property
    def site_summary_csv(self) -> Path:
        return self.results_dir / f"{self.site.id}_summary.csv"

    @property
    def site_provenance_json(self) -> Path:
        return self.results_dir / f"{self.site.id}_provenance.json"

    @property
    def dem_dir(self) -> Path:
        return self.repo_root / "data" / "dem" / self.site.id

    @property
    def strong_zones_geojson(self) -> Path:
        return self.vectors_dir / "strong_anomaly_zones.geojson"


def site_paths_for(site: SiteConfig, repo_root: Path) -> SitePaths:
    return SitePaths(repo_root=repo_root.resolve(), site=site)
