"""
    Classes and tools for handling system dynamics
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import List, Dict
import numpy as np


@dataclass
class DiscreteFunction:
    offset: int
    values: List[float]

    def scale_by(self, factor: float) -> DiscreteFunction:
        return DiscreteFunction(self.offset, [y * factor for y in self.values])

    def mean_filter(self, window_size: int) -> DiscreteFunction:
        values = np.array(self.values)
        filtered = []
        for i in range(len(values)):
            start = i - int(window_size/2)
            end = start + window_size
            start = max(0, start)
            end = min(len(values), end)
            sliced: np.array = values[start:end]
            filtered.append(sliced.mean())
        return DiscreteFunction(self.offset, filtered)


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

# def default_world_properties() -> WorldProperties:
#     alpha
