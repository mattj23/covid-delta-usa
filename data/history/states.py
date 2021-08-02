import json
import os
import pickle
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass

from datetime import date as Date
from datetime import datetime as DateTime
from states import state_abbrevs, adjacency
import settings

_state_info_cache = os.path.join(settings.cache_folder, "state_info_cache.pickle")
_state_history_cache = os.path.join(settings.cache_folder, "state_history_cache.pickle")


@dataclass
class StateInfo:
    name: str
    population: int
    adjacent: List[str]


@dataclass
class DailyRecord:
    total_cases: int
    new_cases: int


@dataclass
class PlottableHistory:
    dates: List[Date]
    total_known_cases: List
    new_known_cases: List


@dataclass
class StateHistory:
    name: str
    records: Dict[Date, DailyRecord]

    def get_in_range(self, start: Date, end: Date) -> Dict[Date, DailyRecord]:
        return {d: v for d, v in self.records.items() if start >= d >= end}

    def get_plottable(self, start: Optional[Date] = None, end: Optional[Date] = None) -> PlottableHistory:
        dates = sorted(self.records.keys())
        records = []
        for date in dates:
            if start is not None and date < start:
                continue
            if end is not None and date > end:
                continue
            r = self.records[date]
            records.append((date, r.total_cases, r.new_cases))
        # For clarity
        dates, total_cases, new_cases = zip(*records)
        return PlottableHistory(dates, total_cases, new_cases)


def load_state_info() -> Dict[str, StateInfo]:
    if not os.path.exists(_state_info_cache):
        data = _get_state_info()
        results = {}
        for k, v in data.items():
            results[k] = StateInfo(**v)

        with open(_state_info_cache, "wb") as handle:
            pickle.dump(results, handle)

    with open(_state_info_cache, "rb") as handle:
        return pickle.load(handle)


def load_state_histories() -> Dict[str, StateHistory]:
    if not os.path.exists(_state_history_cache):
        data = _get_all_histories_no_cache()
        with open(_state_history_cache, "wb") as handle:
            pickle.dump(data, handle)

    with open(_state_history_cache, "rb") as handle:
        return pickle.load(handle)


def _get_state_info() -> Dict:
    all_states = {}
    for state in state_abbrevs():
        state_file_path = os.path.join(settings.cache_folder, f"covid_act_now-{state}.json")
        with open(state_file_path, "r") as handle:
            data = json.load(handle)

        population = data['population']
        all_states[state] = {
            "name": state,
            "population": population,
            "adjacent": adjacency(state),
        }
    return all_states


def _get_all_histories_no_cache() -> Dict[str, StateHistory]:
    all_data = {}

    for state in state_abbrevs():
        all_data[state] = StateHistory(state, {})
        state_file_path = os.path.join(settings.cache_folder, f"covid_act_now-{state}.json")
        with open(state_file_path, "r") as handle:
            data = json.load(handle)

        for record in data['actualsTimeseries']:
            record: Dict[str, Any]
            date = DateTime.strptime(record['date'], "%Y-%m-%d").date()
            total_cases = record['cases']

            # Some dates seem to be missing case values, in which case we'll remove them rather than have blank
            # entries in the history
            if total_cases is None or (total_cases <= 0 and date.year >= 2021):
                continue

            all_data[state].records[date] = DailyRecord(
                total_cases=0 if record["cases"] is None else record["cases"],
                new_cases=0 if record["newCases"] is None else record["newCases"],
            )

    return all_data


