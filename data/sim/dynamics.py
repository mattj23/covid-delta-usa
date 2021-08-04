"""
    Classes and tools for handling system dynamics
"""

from dataclasses import dataclass, asdict
from typing import List, Dict


@dataclass
class DiscreteFunction:
    offset: int
    values: List[float]


@dataclass
class VariantProperties:
    incubation: List[float]
    infectivity: DiscreteFunction
    vax_immunity: DiscreteFunction
    natural_immunity: DiscreteFunction


@dataclass
class WorldProperties:
    alpha: VariantProperties
    delta: VariantProperties


def prepare_variant_properties(v: VariantProperties) -> Dict:
    return {
        "incubation": v.incubation,
        "infectivity": asdict(v.infectivity),
        "vax_immunity": asdict(v.vax_immunity),
        "natural_immunity": asdict(v.natural_immunity)
    }


def prepare_world_properties(w: WorldProperties) -> Dict:
    return {
        "alpha": prepare_variant_properties(w.alpha),
        "delta": prepare_variant_properties(w.delta)
    }


def default_world_properties() -> WorldProperties:
