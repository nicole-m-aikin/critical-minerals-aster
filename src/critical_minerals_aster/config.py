"""Configuration for one study site (YAML-driven)."""

from __future__ import annotations

from dataclasses import dataclass
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
class SiteConfig:
    id: str
    name: str
    bbox_wgs84: BBox
    granule_id: str
    layout: Literal["flat", "nested"] = "flat"
    classification: ClassificationParams | None = None
    temporal_start: str = "2010-01-01"
    temporal_end: str = "2023-12-31"


def _as_classification(obj: Any) -> ClassificationParams:
    if obj is None:
        return ClassificationParams()
    if isinstance(obj, ClassificationParams):
        return obj
    if isinstance(obj, dict):
        return ClassificationParams(**obj)
    raise TypeError(f"Invalid classification config: {type(obj)}")


def load_site_config(path: Union[str, Path]) -> SiteConfig:
    path = Path(path)
    raw = yaml.safe_load(path.read_text())

    temporal = raw.get("temporal") or {}
    return SiteConfig(
        id=raw["id"],
        name=raw["name"],
        bbox_wgs84=tuple(raw["bbox_wgs84"]),
        granule_id=raw["granule_id"],
        layout=raw.get("layout", "flat"),
        classification=_as_classification(raw.get("classification")),
        temporal_start=temporal.get("start", "2010-01-01"),
        temporal_end=temporal.get("end", "2023-12-31"),
    )
