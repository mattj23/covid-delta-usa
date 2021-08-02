import json
from typing import Dict, List, Optional, Callable, Tuple

import numpy

import settings
import subprocess
import time
from datetime import date as Date
from dataclasses import dataclass

from sim.program_input import ProgramInput


@dataclass
class ValueDistribution:
    mean: List[float]
    upper: List[float]
    lower: List[float]


@dataclass
class PlottableSteps:
    dates: List[Date]
    infections: ValueDistribution
    new_infections: ValueDistribution
    vaccinated: ValueDistribution
    vaccine_saves: ValueDistribution


@dataclass
class StepResult:
    date: Date
    infections: int
    new_infections: int
    positive_tests: int
    vaccinated: int
    vaccine_saves: int


@dataclass
class SimulationResult:
    run_time: float
    results: Dict[str, List[List[StepResult]]]

    def get_plottable(self, state: str, start: Optional[Date] = None, end: Optional[Date] = None) -> PlottableSteps:
        state_data = self.results[state]

        dates = [step.date for step in state_data[0]]
        if start is not None:
            dates = [d for d in dates if d >= start]
        if end is not None:
            dates = [d for d in dates if d <= end]

        return PlottableSteps(dates=dates,
                              infections=self._get_high_mean_lo(state, lambda step: step.infections, start, end),
                              new_infections=self._get_high_mean_lo(state, lambda step: step.new_infections, start, end),
                              vaccinated=self._get_high_mean_lo(state, lambda step: step.vaccinated, start, end),
                              vaccine_saves=self._get_high_mean_lo(state, lambda step: step.vaccine_saves, start, end)
                              )

    def _get_high_mean_lo(self, state: str,
                          extractor: Callable,
                          start: Optional[Date] = None,
                          end: Optional[Date] = None,
                          ) -> ValueDistribution:
        state_data = self.results[state]
        distributions = [[] for _ in state_data[0]]
        for run in state_data:
            for i, step in enumerate(run):
                if start is not None and step.date < start:
                    continue
                if end is not None and step.date > end:
                    continue
                distributions[i].append(extractor(step))

        mean_ = []
        upper_ = []
        lower_ = []
        for dist in distributions:
            d = numpy.array(dist)
            mean_.append(numpy.mean(d))
            upper_.append(numpy.quantile(d, 0.95))
            lower_.append(numpy.quantile(d, 0.05))

        return ValueDistribution(mean_, upper_, lower_)


class Simulator:
    def __init__(self, input_data: ProgramInput, input_file=None):
        self.input_data = input_data
        self.input_file = settings.default_input_file if input_file is None else input_file

    def run(self) -> SimulationResult:
        self.input_data.write(self.input_file)
        start_time = time.time()
        command = [settings.binary_path, self.input_file]
        process = subprocess.Popen(command)
        process.communicate()
        end_time = time.time()

        result = self._load_results()
        return SimulationResult(end_time - start_time, result)

    def _load_results(self):
        with open(self.input_data.output_file, "r") as handle:
            raw_data = json.load(handle)
        results = {}
        for row in raw_data:
            state_name = row["name"]
            if state_name not in results:
                results[state_name] = []

            run_result = []
            for d in row["results"]:
                run_result.append(StepResult(
                    date=Date(d["year"], d["month"], d["day"]),
                    infections=d["infections"],
                    new_infections=d["new_infections"],
                    positive_tests=d["positive_tests"],
                    vaccinated=d["vaccinated"],
                    vaccine_saves=d["vaccine_saves"]
                ))

            results[state_name].append(sorted(run_result, key=lambda x: x.date))

        return results
