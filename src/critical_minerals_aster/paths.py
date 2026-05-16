"""Filesystem paths for site data and figures."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from critical_minerals_aster.config import SiteConfig


@dataclass(frozen=True)
class SitePaths:
    """Resolved directories for ASTER rasters, vectors, and figures."""

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


def site_paths_for(site: SiteConfig, repo_root: Path) -> SitePaths:
    return SitePaths(repo_root=repo_root.resolve(), site=site)
