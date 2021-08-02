"""
    Module for loading historical vaccine data
"""
import json
import os
import pickle
from datetime import date as Date
from datetime import datetime as DateTime

from states import state_abbrevs
from typing import Dict, List, Any
from dataclasses import dataclass

import settings

_vaccine_cache_file = os.path.join(settings.cache_folder, "state_vaccine_history_cache.pickle")


@dataclass
class VaccineDailyRecord:
    distributed: int
    completed: int
    initiated: int
    administered: int


@dataclass
class StateVaccineHistory:
    name: str
    records: Dict[Date, VaccineDailyRecord]


def load_vaccine_histories() -> Dict[str, StateVaccineHistory]:
    if not os.path.exists(_vaccine_cache_file):
        data = _get_vaccine_histories_no_cache()
        with open(_vaccine_cache_file, "wb") as handle:
            pickle.dump(data, handle)

    with open(_vaccine_cache_file, "rb") as handle:
        return pickle.load(handle)


def _get_vaccine_histories_no_cache() -> Dict[str, StateVaccineHistory]:
    all_data = {}

    for state in state_abbrevs():
        all_data[state] = StateVaccineHistory(state, {})
        state_file_path = os.path.join(settings.cache_folder, f"covid_act_now-{state}.json")
        with open(state_file_path, "r") as handle:
            data = json.load(handle)

        for record in data['actualsTimeseries']:
            record: Dict[str, Any]
            date = DateTime.strptime(record['date'], "%Y-%m-%d").date()

            if "vaccinationsCompleted" in record:
                if record["vaccinationsCompleted"] is None:
                    continue
                all_data[state].records[date] = VaccineDailyRecord(
                    distributed=record["vaccinesDistributed"],
                    completed=record["vaccinationsCompleted"],
                    initiated=record["vaccinationsInitiated"],
                    administered=record["vaccinesAdministered"]
                )

    return all_data


if __name__ == '__main__':
    vax_record = load_vaccine_histories()
    print(vax_record)
