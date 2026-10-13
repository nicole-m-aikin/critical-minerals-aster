"""Early-exploration prospectivity ML (sensor-agnostic feature contract).

Ranks ground for triage — “worth a closer look” — using thermal host-rock
cues plus hyperspectral clay/alteration cues under site- and era-blocked
validation. Not a drill-target or ore detector.
"""

from critical_minerals_aster.prospectivity.adapters import (
    HyperspectralScene,
    ThermalScene,
    load_aster_thermal_scene,
    load_aviris_hyperspectral_scene,
)
from critical_minerals_aster.prospectivity.split import (
    EraBlockSplit,
    era_block_split,
    site_group_kfold,
)

__all__ = [
    "HyperspectralScene",
    "ThermalScene",
    "load_aster_thermal_scene",
    "load_aviris_hyperspectral_scene",
    "EraBlockSplit",
    "era_block_split",
    "site_group_kfold",
]
