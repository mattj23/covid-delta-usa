"""
    Classes and tools for handling system dynamics
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import List, Dict, Tuple, Optional

import numpy
import numpy as np
from matplotlib.axes._axes import Axes
from scipy.interpolate import interp1d


@dataclass
class DiscreteFunction:
    offset: int
    values: List[float]

    def _to_arrays(self) -> Tuple[np.array, np.array]:
        x, y = [np.array(p) for p in zip(*enumerate(self.values))]
        return x - self.offset, y

    def plot(self, ax: Axes, color: Optional[str], linewidth: Optional[int], label: Optional[str]):
        x, y = self._to_arrays()
        kwargs = {}
        if linewidth is not None:
            kwargs["linewidth"] = linewidth
        if color is not None:
            kwargs["color"] = color
        if label is not None:
            kwargs["label"] = label
        ax.plot(x, y, **kwargs)

    def scale_y(self, factor: float) -> DiscreteFunction:
        return DiscreteFunction(self.offset, [y * factor for y in self.values])

    def scale_x(self, factor_array: List[float]) -> DiscreteFunction:
        assert len(factor_array) == len(self.values)
        x, y = self._to_arrays()
        x_ = np.multiply(x, factor_array)
        f = interp1d(x_, y)
        x__ = np.arange(min(x_), max(x_))
        offset = round(self.offset * factor_array[self.offset])
        return DiscreteFunction(offset, list(f(x__)))

    def scale_x_positive(self, starting_x: int, scale_factor, kind="linear") -> DiscreteFunction:
        x, y = self._to_arrays()
        working = 0
        while working < len(x):
            original = x[working]
            if original >= starting_x:
                span = original - starting_x
                x[working] = original + scale_factor * span
            working += 1

        f = interp1d(x, y, kind=kind)
        x__ = np.arange(min(x), max(x))
        return DiscreteFunction(self.offset, list(f(x__)))

    def normalize_area(self) -> DiscreteFunction:
        area = sum(self.values)
        return self.scale_y(1/area)

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

