"""
    Module for working with the estimated undetected cases of Covid-19 based on the data from covidestim.org

    This module will load the data from the covid-estim.csv file produced by the fetch_data.py script. Make sure you've
    used that script to retrieve the data before trying to use this.
"""
from typing import Dict, Tuple, Optional, List

import csv
import datetime
from datetime import date as Date
from dataclasses import dataclass
import os
import settings

from states import state_to_abbrev

_estimate_file = os.path.join(settings.cache_folder, "covid-estim.csv")
_to_abbrev = state_to_abbrev()


@dataclass
class DailyEstimate:
    cases: int
    diagnoses: int
    infections: int
    total_cases: int
    total_diagnoses: int
    total_infections: int


@dataclass
class PlottableEstimates:
    dates: List[Date]
    cases: List[int]
    diagnoses: List[int]
    infections: List[int]
    total_cases: List[int]
    total_diagnoses: List[int]
    total_infections: List[int]


@dataclass
class StateEstimates:
    name: str
    estimates: Dict[Date, DailyEstimate]

    def get_in_range(self, start: Date, end: Date) -> Dict[Date, DailyEstimate]:
        return {d: v for d, v in self.estimates.items() if start >= d >= end}

    def get_plottable(self, start: Optional[Date] = None, end: Optional[Date] = None) -> PlottableEstimates:
        dates = sorted(self.estimates.keys())
        records = []
        for date in dates:
            if start is not None and date < start:
                continue
            if end is not None and date > end:
                continue
            r = self.estimates[date]
            records.append((date, r.cases, r.diagnoses, r.infections, r.total_cases,
                            r.total_diagnoses, r.total_infections))
        # For clarity
        dates, cases, diagnoses, infections, total_cases, total_diagnoses, total_infections = zip(*records)
        return PlottableEstimates(dates, cases, diagnoses, infections, total_cases, total_diagnoses, total_infections)


def load_state_estimates() -> Dict[str, StateEstimates]:
    values = {}
    with open(_estimate_file, "r") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row['state'] not in _to_abbrev:
                continue

            state = _to_abbrev[row["state"]]
            if state not in values:
                values[state] = {}

            date = datetime.datetime.strptime(row['date'], "%Y-%m-%d").date()
            estimate = DailyEstimate(cases=int(float(row["cases.fitted"])),
                                     diagnoses=int(float(row["diagnoses"])),
                                     infections=int(float(row["infections"])),
                                     total_cases=0,
                                     total_diagnoses=0,
                                     total_infections=0)
            values[state][date] = estimate

    result = {}
    for k, v in values.items():
        # Compute the cumulative totals
        result[k] = _integrate_estimates(StateEstimates(k, v))

    return result


def _integrate_estimates(state_estimates: StateEstimates) -> StateEstimates:
    all_dates = sorted(state_estimates.estimates.keys())
    total_cases = 0
    total_diagnoses = 0
    total_infections = 0
    for date in all_dates:
        total_cases += state_estimates.estimates[date].cases
        total_diagnoses += state_estimates.estimates[date].diagnoses
        total_infections += state_estimates.estimates[date].infections

        state_estimates.estimates[date].total_cases = total_cases
        state_estimates.estimates[date].total_diagnoses = total_diagnoses
        state_estimates.estimates[date].total_infections = total_infections

    return state_estimates


def main():
    estimates = load_state_estimates()
    print(estimates)


if __name__ == '__main__':
    main()
