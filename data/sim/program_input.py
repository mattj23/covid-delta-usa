"""
    Simulation tools
"""

import json
from enum import IntEnum
from datetime import date as Date
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional

from history.estimates import StateEstimates
from history.states import StateHistory, StateInfo
from history.vaccine import StateVaccineHistory
from sim.dynamics import WorldProperties, prepare_world_properties


class ProgramMode(IntEnum):
    Simulate = 1
    FindContactProb = 2


@dataclass
class ProgramOptions:
    full_history: bool
    expensive_stats: bool
    mode: ProgramMode


@dataclass
class ProgramInput:
    """
    Settings for running the simulator program, contains all configuration and setup
    """
    output_file: str
    state: str
    start_day: Date
    end_day: Date
    contact_prob: float
    world_properties: WorldProperties
    state_info: Dict[str, StateInfo]
    contact_day_interval: int = 1
    options: Optional[ProgramOptions] = None
    population_scale: Optional[int] = 10
    run_count: Optional[int] = 1
    test_history: Optional[Dict[str, StateHistory]] = None
    infected_history: Optional[Dict[str, StateEstimates]] = None
    vax_history: Optional[Dict[str, StateVaccineHistory]] = None
    variant_history: Optional[Dict[str, List[Dict]]] = None

    def to_str(self) -> str:
        if self.options is None:
            self.options = ProgramOptions(False, False, ProgramMode.Simulate)
        output = {
            "output_file": self.output_file,
            "state": self.state,
            "world_properties": prepare_world_properties(self.world_properties),
            "start_day": self.start_day.strftime("%Y-%m-%d"),
            "end_day": self.end_day.strftime("%Y-%m-%d"),
            "contact_probability": self.contact_prob,
            "population_scale": self.population_scale,
            "contact_day_interval": self.contact_day_interval,
            "run_count": self.run_count,
            "options": asdict(self.options),
            "infected_history": _prep_estimates(self.infected_history),
            "vax_history": _prep_vaccines(self.vax_history),
            "test_history": _prep_history(self.test_history),
            "state_info": _prep_state_info(self.state_info),
            "variant_history": _prep_variant_history(self.variant_history)
        }
        return json.dumps(output, indent=2)


def _prep_variant_history(variant_data: Optional[Dict[str, List[Dict]]]) -> Dict[str, List]:
    all_prepared = {}
    if variant_data is None:
        return all_prepared

    for state, data in variant_data.items():
        prepared = []
        for row in data:
            prepared.append({"date": row["date"].strftime("%Y-%m-%d"), "variants": row["variants"]})
        all_prepared[state] = prepared

    return all_prepared


def _prep_state_info(data: Dict[str, StateInfo]) -> Dict:
    prepared = {}
    if data is None:
        return prepared

    for state, info in data.items():
        prepared[state] = {
            "population": info.population,
            "adjacent": info.adjacent
        }
    return prepared


def _prep_vaccines(data: Optional[Dict[str, StateVaccineHistory]]) -> Dict:
    prepared = {}
    if data is None:
        return prepared

    for state, vax in data.items():
        prepared[state] = {}
        for d, e in vax.records.items():
            prepared[state][d.strftime("%Y-%m-%d")] = {
                "total_completed_vax": e.completed
            }
    return prepared


def _prep_history(data: Optional[Dict[str, StateHistory]]) -> Dict:
    prepared = {}
    if data is None:
        return prepared

    for state, hist in data.items():
        prepared[state] = {}
        for d, e in hist.records.items():
            prepared[state][d.strftime("%Y-%m-%d")] = {
                "total_known_cases": e.total_cases
            }
    return prepared


def _prep_estimates(data: Optional[Dict[str, StateEstimates]]) -> Dict:
    prepared = {}
    if data is None:
        return prepared

    for state, est in data.items():
        prepared[state] = {}
        for d, e in est.estimates.items():
            prepared[state][d.strftime("%Y-%m-%d")] = {
                "total_infections": e.total_infections,
                "total_cases": e.total_cases
            }
    return prepared
