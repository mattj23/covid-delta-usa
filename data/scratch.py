import os
import pickle
import settings
from datetime import date as Date

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes

from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories)
from sim import ProgramInput, Simulator, default_world_properties, ContactSearchResult, SimulationResult
from sim.dynamics import DiscreteFunction


def efficacy_curves():
    fig: Figure = plt.figure(figsize=(16, 10))
    ax: Axes = fig.subplots()
    ax.set_title("Vaccine Efficacy Curves")
    ax.set_xlabel("Days since Dose 1 (Assuming second dose follows)")
    ax.set_ylabel("Efficacy")

    from sim.world_defaults import (pfizer_alpha_efficacy, pfizer_delta_efficacy_lopez_uk,
                                    pfizer_delta_efficacy_israeli_moh)
    pfizer_alpha_efficacy().plot(ax, "forestgreen", 3, "Pfizer, vs Alpha (Thomas et al. 2021)")
    pfizer_delta_efficacy_lopez_uk().plot(ax, "slateblue", 3, "Pfizer, vs Delta (Thomas et al. scaled to Lopez Bernal et. al)")
    pfizer_delta_efficacy_israeli_moh().plot(ax, "dodgerblue", 3, "Pfizer, vs Delta (Israeli MOH trend spliced onto scaled Lopez Bernal)", "dotted")



    ax.legend(loc="lower right")
    fig.show()





def fl_anomolous_data():
    input_data = ProgramInput(output_file=settings.default_output_file,
                              state="FL",
                              world_properties=default_world_properties(),
                              start_day=Date(2021, 8, 7),
                              end_day=Date(2021, 12, 1),
                              state_info=load_state_info(),
                              population_scale=100,
                              contact_prob=1.75,
                              test_history=load_state_histories(),
                              vax_history=load_vaccine_histories(),
                              infected_history=load_state_estimates(),
                              variant_history=load_variant_history(),
                              run_count=50)

    with open(os.path.join(settings.cache_folder, "data_2021-08-15.pickle"), "rb") as handle:
        old_input_data = pickle.load(handle)

    plot_start = Date(2021, 6, 1)  # input_data.start_day - TimeDelta(days=5)
    # plot_start = input_data.start_day
    plot_end = input_data.end_day

    pl0 = old_input_data.infected_history[input_data.state].get_plottable(plot_start, plot_end)
    pl1 = input_data.infected_history[input_data.state].get_plottable(plot_start, plot_end)

    fig: Figure = plt.figure(figsize=(16, 10))
    ax: Axes = fig.subplots()
    ax.set_title("FL Covidestim.org Anomalous Model Output")
    ax.set_xlabel("Date")
    ax.set_ylabel("Estimated Daily Infections")
    ax.plot(pl0.dates, pl0.infections, color="lightblue", linewidth=4, label=f"Covidestim.org Data Ending {pl0.dates[-1]}")
    ax.plot(pl1.dates, pl1.infections, color="lightgreen", linewidth=4, label=f"Covidestim.org Data Ending {pl1.dates[-1]}")
    # ax.plot(pl.dates, n_.total_vaccinated.mean, color="lightgreen", linewidth=4, label="New Binary")
    ax.legend()

    fig.show()


if __name__ == '__main__':
    efficacy_curves()
