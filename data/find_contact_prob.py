from typing import Optional, Dict, Tuple

import numpy as np
from datetime import date as Date

import settings
from sim.dynamics import VariantProperties, DiscreteFunction, WorldProperties
from sim import default_world_properties, ProgramInput, Simulator, ContactSearchResult
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories, StateEstimates, DailyEstimate, StateInfo)


def find_contact_prob(state: str, start_date: Date, end_date: Date, interval: int, world: WorldProperties,
                      state_info: Optional[Dict[str, StateInfo]] = None,
                      vax_history: Optional = None,
                      variant_history: Optional = None,
                      infected_history: Optional = None) -> ContactSearchResult:
    state_info = state_info if state_info is not None else load_state_info()
    scale = int(state_info[state].population / 500000)

    input_data = ProgramInput(output_file=settings.default_output_file,
                              state=state,
                              world_properties=world,
                              start_day=start_date,
                              end_day=end_date,
                              contact_day_interval=interval,
                              contact_prob=0,
                              state_info=state_info,
                              population_scale=scale,
                              vax_history=vax_history if vax_history is not None else load_vaccine_histories(),
                              variant_history=variant_history if variant_history is not None else load_variant_history(),
                              infected_history=infected_history if infected_history is not None else load_state_estimates(),
                              run_count=200)

    simulator = Simulator(input_data, settings.default_input_file)
    return simulator.find_contact_prob()


def find_days_contact_prob(state: str, date: Date, world: WorldProperties,
                           state_info: Optional[Dict[str, StateInfo]] = None,
                           vax_history: Optional = None,
                           variant_history: Optional = None,
                           infected_history: Optional = None) -> float:
    result = find_contact_prob(state, date, date, 1, world, state_info, vax_history, variant_history, infected_history)
    return result.probabilities[0]
