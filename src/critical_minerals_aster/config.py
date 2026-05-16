"""Configuration for one study site (YAML-driven)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Tuple, Union

import yaml

BBox = Tuple[float, float, float, float]


@dataclass
class ClassificationParams:
    low_pct: float = 70.0
    high_pct: float = 90.0
    strong_score_min: int = 3


@dataclass
class StructureLayer:
    path: str
    type: Literal["faults", "contacts", "folds"] = "faults"
    buffer_m: float = 500.0


@dataclass
class SiteConfig:
    id: str
    name: str
    bbox_wgs84: BBox
    granule_id: str | None = None
    layout: Literal["flat", "nested"] = "flat"
    buffer_deg: float = 0.0
    classification: ClassificationParams | None = None
    temporal_start: str = "2010-01-01"
    temporal_end: str = "2023-12-31"
    structure_layers: list[StructureLayer] = field(default_factory=list)


def _as_classification(obj: Any) -> ClassificationParams:
    if obj is None:
        return ClassificationParams()
    if isinstance(obj, ClassificationParams):
        return obj
    if isinstance(obj, dict):
        return ClassificationParams(**obj)
    raise TypeError(f"Invalid classification config: {type(obj)}")


def _as_structure_layers(raw: Any) -> list[StructureLayer]:
    if not raw:
        return []
    layers: list[StructureLayer] = []
    for item in raw:
        if isinstance(item, StructureLayer):
            layers.append(item)
        elif isinstance(item, dict):
            layers.append(StructureLayer(**item))
        else:
            raise TypeError(f"Invalid structure layer entry: {item!r}")
    return layers


def search_bbox(site: SiteConfig) -> BBox:
    """WGS84 bbox expanded by buffer_deg for granule search."""
    if site.buffer_deg <= 0:
        return site.bbox_wgs84
    lon0, lat0, lon1, lat1 = site.bbox_wgs84
    b = site.buffer_deg
    return (lon0 - b, lat0 - b, lon1 + b, lat1 + b)


def load_site_config(path: Union[str, Path]) -> SiteConfig:
    path = Path(path)
    raw = yaml.safe_load(path.read_text())

    temporal = raw.get("temporal") or {}
    granule_id = raw.get("granule_id")
    if granule_id in (None, "", "null"):
        granule_id = None

    return SiteConfig(
        id=raw["id"],
        name=raw["name"],
        bbox_wgs84=tuple(raw["bbox_wgs84"]),
        granule_id=granule_id,
        layout=raw.get("layout", "flat"),
        buffer_deg=float(raw.get("buffer_deg", 0.0)),
        classification=_as_classification(raw.get("classification")),
        temporal_start=temporal.get("start", "2010-01-01"),
        temporal_end=temporal.get("end", "2023-12-31"),
        structure_layers=_as_structure_layers(raw.get("structure_layers")),
    )


def list_site_ids(sites_dir: Union[str, Path]) -> list[str]:
    """Site ids from sites/index.yaml or all *.yaml except index."""
    sites_dir = Path(sites_dir)
    index_path = sites_dir / "index.yaml"
    if index_path.is_file():
        raw = yaml.safe_load(index_path.read_text()) or {}
        return list(raw.get("sites", []))
    return sorted(
        p.stem for p in sites_dir.glob("*.yaml") if p.stem != "index"
    )


def load_site_by_id(site_id: str, sites_dir: Union[str, Path]) -> SiteConfig:
    return load_site_config(Path(sites_dir) / f"{site_id}.yaml")
