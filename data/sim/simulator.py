from __future__ import annotations

import json
import hashlib
import os
import pickle
from typing import Dict, List, Optional, Callable, Tuple, Union, Any

import numpy

import settings
import subprocess
import time
from datetime import date as Date
from datetime import timedelta as TimeDelta
from dataclasses import dataclass

from sim.program_input import ProgramInput, ProgramMode, ProgramOptions

_reference_date = Date(2019, 1, 1)


def from_integer_date(d: int) -> Date:
    return _reference_date + TimeDelta(days=d)


def _stream_process(process):
    go = process.poll() is None
    for line in process.stdout:
        print(line)
    return go


@dataclass
class ValueDistribution:
    mean: List[float]
    upper: List[float]
    lower: List[float]

    def __truediv__(self, other: ValueDistribution) -> ValueDistribution:
        t = ValueDistribution([], [], [])
        assert len(self.mean) == len(other.mean)
        for i in range(len(self.mean)):
            t.mean.append(self.mean[i] / other.mean[i])
            t.upper.append(self.upper[i] / other.upper[i])
            t.lower.append(self.lower[i] / other.lower[i])
        return t

    def __mul__(self, other: Union[float, ValueDistribution]) -> ValueDistribution:
        if isinstance(other, float):
            return ValueDistribution(numpy.array(self.mean) * other,
                                     numpy.array(self.upper) * other,
                                     numpy.array(self.lower) * other)

        if isinstance(other, ValueDistribution):
            t = ValueDistribution([], [], [])
            assert len(self.mean) == len(other.mean)
            for i in range(len(self.mean)):
                t.mean.append(self.mean[i] * other.mean[i])
                t.upper.append(self.upper[i] * other.upper[i])
                t.lower.append(self.lower[i] * other.lower[i])
            return t


@dataclass
class PlottableSteps:
    dates: List[Date]
    total_infections: ValueDistribution
    new_infections: ValueDistribution
    total_vaccinated: ValueDistribution
    never_infected: ValueDistribution
    total_delta_infections: ValueDistribution
    new_delta_infections: ValueDistribution
    total_alpha_infections: ValueDistribution
    new_alpha_infections: ValueDistribution
    reinfections: ValueDistribution
    new_reinfections: ValueDistribution
    vaccine_saves: ValueDistribution
    new_vaccine_saves: ValueDistribution
    natural_saves: ValueDistribution
    new_natural_saves: ValueDistribution
    vaccinated_infections: ValueDistribution
    virus_carriers: ValueDistribution
    population_infectiousness: ValueDistribution


@dataclass
class StepResult:
    date: Date
    total_infections: int
    total_vaccinated: int
    never_infected: int
    total_delta_infections: int
    total_alpha_infections: int
    reinfections: int
    vaccinated_infections: int
    vaccine_saves: int
    natural_saves: int
    virus_carriers: int
    population_infectiousness: float = 0
    new_infections: int = 0
    new_delta_infections: int = 0
    new_alpha_infections: int = 0
    new_reinfections: int = 0
    new_natural_saves: int = 0
    new_vaccine_saves: int = 0

    def set_difference_values(self, other: StepResult):
        self.new_infections = self.total_infections - other.total_infections
        self.new_delta_infections = self.total_delta_infections - other.total_delta_infections
        self.new_alpha_infections = self.total_alpha_infections - other.total_alpha_infections
        self.new_reinfections = self.reinfections - other.reinfections
        self.new_natural_saves = self.natural_saves - other.natural_saves
        self.new_vaccine_saves = self.vaccine_saves - other.vaccine_saves


@dataclass
class ContactSearchResult:
    days: List[Date]
    probabilities: List[float]
    stdevs: List[float]


@dataclass
class SimulationResult:
    run_time: float
    results: Union[Dict[str, List[List[StepResult]]], List[float]]

    def get_plottable(self, state: str, start: Optional[Date] = None, end: Optional[Date] = None) -> PlottableSteps:
        state_data = self.results[state]

        dates = [step.date for step in state_data[0]]
        if start is not None:
            dates = [d for d in dates if d >= start]
        if end is not None:
            dates = [d for d in dates if d <= end]

        return PlottableSteps(
            dates=dates,
            total_infections=self._get_high_mean_lo(state, lambda step: step.total_infections, start, end),
            new_infections=self._get_high_mean_lo(state, lambda step: step.new_infections, start, end),
            total_vaccinated=self._get_high_mean_lo(state, lambda step: step.total_vaccinated, start, end),
            never_infected=self._get_high_mean_lo(state, lambda step: step.never_infected, start, end),
            total_delta_infections=self._get_high_mean_lo(state, lambda step: step.total_delta_infections, start, end),
            total_alpha_infections=self._get_high_mean_lo(state, lambda step: step.total_alpha_infections, start, end),
            vaccine_saves=self._get_high_mean_lo(state, lambda step: step.vaccine_saves, start, end),
            natural_saves=self._get_high_mean_lo(state, lambda step: step.natural_saves, start, end),
            new_vaccine_saves=self._get_high_mean_lo(state, lambda step: step.new_vaccine_saves, start, end),
            new_natural_saves=self._get_high_mean_lo(state, lambda step: step.new_natural_saves, start, end),
            vaccinated_infections=self._get_high_mean_lo(state, lambda step: step.vaccinated_infections, start, end),
            reinfections=self._get_high_mean_lo(state, lambda step: step.reinfections, start, end),
            new_delta_infections=self._get_high_mean_lo(state, lambda step: step.new_delta_infections, start, end),
            new_alpha_infections=self._get_high_mean_lo(state, lambda step: step.new_alpha_infections, start, end),
            virus_carriers=self._get_high_mean_lo(state, lambda step: step.virus_carriers, start, end),
            population_infectiousness=self._get_high_mean_lo(state, lambda step: step.population_infectiousness, start,
                                                             end),
            new_reinfections=self._get_high_mean_lo(state, lambda step: step.new_reinfections, start, end)
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
            if not dist:
                continue
            d = numpy.array(dist)
            mean_.append(numpy.mean(d))
            if max(d) - min(d) == 0:
                upper_.append(d[0])
                lower_.append(d[0])
            else:
                upper_.append(numpy.quantile(d, 0.95))
                lower_.append(numpy.quantile(d, 0.05))

        return ValueDistribution(mean_, upper_, lower_)


class Simulator:
    def __init__(self, input_data: ProgramInput, input_file=None, no_cache_results=False):
        self.input_data = input_data
        self.input_file = settings.default_input_file if input_file is None else input_file
        self.no_cache_results = no_cache_results

        self._input_text: Optional[str] = None
        self._cache_path: Optional[str] = None

    def _set_cache_info(self):
        if self._input_text is not None and self._cache_path is not None:
            return

        self._input_text = self.input_data.to_str()
        digest = hashlib.sha1(self._input_text.encode()).hexdigest()
        self._cache_path = os.path.join(settings.cache_folder, "results", f"{digest}.pickle")

    def _clear_cache_info(self):
        self._input_text = None
        self._cache_path = None

    def _check_cache(self) -> Optional[Any]:
        self._set_cache_info()

        if os.path.exists(self._cache_path):
            with open(self._cache_path, "rb") as handle:
                return pickle.load(handle)
        return None

    def _write_input_text(self):
        self._set_cache_info()
        with open(self.input_file, "w") as handle:
            handle.write(self._input_text)

    def _cache_results(self, result_obj):
        assert self._cache_path is not None
        base_dir, _ = os.path.split(self._cache_path)
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
        with open(self._cache_path, "wb") as handle:
            pickle.dump(result_obj, handle)

    def run(self, full_history: bool = False, expensive_stats: bool = False, use_cache=False) -> SimulationResult:
        self.input_data.options = ProgramOptions(full_history, expensive_stats, ProgramMode.Simulate)

        # Check for cached results
        if use_cache:
            final_result = self._check_cache()
            if final_result is not None:
                self._clear_cache_info()
                print("Retrieved results from cache")
                return SimulationResult(0, final_result)

        self._write_input_text()

        start_time = time.time()
        command = [settings.binary_path, self.input_file]
        process = subprocess.Popen(command)
        process.communicate()
        end_time = time.time()

        result = self._load_results()

        if not self.no_cache_results:
            self._cache_results(result)

        self._clear_cache_info()
        return SimulationResult(end_time - start_time, result)

    def find_contact_prob(self, use_cache=False) -> ContactSearchResult:
        self.input_data.options = ProgramOptions(False, False, ProgramMode.FindContactProb)

        if use_cache:
            final_result = self._check_cache()
            if final_result is not None:
                final_result: ContactSearchResult
                self._clear_cache_info()
                print("Retrieved results from cache")
                return final_result

        self._write_input_text()

        command = [settings.binary_path, self.input_file]
        with subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) as process:
            for line in process.stdout:
                print(line.decode('utf-8').strip("\n"))

        result: Dict = self._load_results()
        days = [from_integer_date(d) for d in result['days']]
        result['days'] = days

        search_result = ContactSearchResult(**result)

        if not self.no_cache_results:
            self._cache_results(search_result)

        return search_result

    def _load_simulation_results(self, raw_data) -> Dict[str, List[List[StepResult]]]:
        results = {}
        for row in raw_data:
            state_name = row["name"]
            if state_name not in results:
                results[state_name] = []

            run_result = []
            for d in row["results"]:
                d_ = dict(d)
                del d_["day"]
                d_["date"] = from_integer_date(d["day"])
                run_result.append(StepResult(**d_))

            run_result.sort(key=lambda x: x.date)
            for i in range(len(run_result) - 1):
                run_result[i + 1].set_difference_values(run_result[i])

            # add all but the first value, which is before the start day
            results[state_name].append(run_result[1:])

        return results

    def _load_results(self):
        with open(self.input_data.output_file, "r") as handle:
            raw_data = json.load(handle)

        if self.input_data.options.mode == ProgramMode.Simulate:
            return self._load_simulation_results(raw_data)

        if self.input_data.options.mode == ProgramMode.FindContactProb:
            return raw_data
