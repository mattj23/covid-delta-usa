from typing import Optional, Dict, Tuple

import numpy as np
from datetime import date as Date

import settings
from sim.dynamics import VariantProperties, DiscreteFunction, WorldProperties
from sim import default_world_properties, ProgramInput, Simulator
from sim.program_input import ProgramMode, ProgramOptions
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories, StateEstimates, DailyEstimate, StateInfo)

import scipy.stats as stats


def find_contact_prob(state: str, date: Date, world: WorldProperties, guess_prob: Optional[float] = None,
                      state_info: Optional[Dict[str, StateInfo]] = None,
                      vax_history: Optional = None,
                      variant_history: Optional = None,
                      infected_history: Optional = None) -> float:
    upper = 2.0
    lower = 0.25

    if guess_prob is not None:
        upper = guess_prob + 0.5
        lower = guess_prob - 0.5

    state_info = state_info if state_info is not None else load_state_info()
    scale = int(state_info[state].population / 500000)

    input_data = ProgramInput(output_file=settings.default_output_file,
                              options=ProgramOptions(full_history=False,
                                                     expensive_stats=False,
                                                     mode=ProgramMode.FindContactProb),
                              state=state,
                              world_properties=world,
                              start_day=date,
                              end_day=date,
                              contact_prob=0,
                              state_info=state_info,
                              population_scale=scale,
                              vax_history=vax_history if vax_history is not None else load_vaccine_histories(),
                              variant_history=variant_history if variant_history is not None else load_variant_history(),
                              infected_history=infected_history if infected_history is not None else load_state_estimates(),
                              run_count=200)

    results = []
    while True:
        upper_result = _get_error(upper, input_data)
        lower_result = _get_error(lower, input_data)
        results.append(upper_result)
        results.append(lower_result)

        if upper_result[0] < 0:
            upper *= 2
        elif lower_result[0] > 0:
            lower_result *= 0.5
        else:
            break

    tries = 0
    while tries < 10:
        just_above = min([r for r in results if r[0] > 0], key=lambda x: x[0])
        just_below = max([r for r in results if r[0] < 0], key=lambda x: x[0])

        span = just_above[-1] - just_below[-1]
        f = (0 - just_below[0]) / (just_above[0] - just_below[0])
        guess = f * span + just_below[-1]

        above_span = just_above[-1] - guess
        below_span = guess - just_below[-1]
        upper = guess + above_span * 0.3
        lower = guess - below_span * 0.3

        upper_result = _get_error(upper, input_data)
        lower_result = _get_error(lower, input_data)
        results.append(upper_result)
        results.append(lower_result)

        # Check if 0 error is included in both the lower result's +1 stdev and the upper
        # result's -1 stdev, and if so we're close enough to just return with an
        # interpolated guess
        if lower_result[0] + lower_result[1] > 0 and upper_result[0] - upper_result[1] < 0:
            return guess

        tries += 1
    return guess


def _get_error(p: float, input_data: ProgramInput) -> Tuple[float, float, float]:
    input_data.contact_prob = p
    simulator = Simulator(input_data, settings.default_input_file)
    result = simulator.run()
    values = np.array(result.results)
    *args, loc, scale = stats.norm.fit(values)
    return loc, scale, p


if __name__ == '__main__':
    find_contact_prob("CO", Date(2020, 12, 1), default_world_properties())
