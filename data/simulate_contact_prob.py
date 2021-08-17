"""
    First simulation of the system
"""
import os
import pickle
from datetime import date as Date
from datetime import timedelta as TimeDelta

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes

import settings
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories)
from sim import ProgramInput, Simulator, default_world_properties, ContactSearchResult
from find_contact_prob import find_contact_prob


def main():
    cache_file = os.path.join(settings.cache_folder, "temp.pickle")
    state = "FL"

    infected_history = load_state_estimates()
    if not os.path.exists(cache_file):
        result = find_contact_prob(state, Date(2020, 5, 1), Date(2021, 8, 7), 3, default_world_properties())
        with open(cache_file, "wb") as handle:
            pickle.dump(result, handle)

    with open(cache_file, "rb") as handle:
        contact_results: ContactSearchResult = pickle.load(handle)

    fig: Figure = plt.figure(figsize=(16, 10))

    ax0 = fig.subplots(1)
    ax0: Axes
    ax0.set_title(f"Contact Ratio for {state}, Delta 2.5x modifier")
    ax0.set_xlabel("Date")
    ax0.set_ylabel(f"Contact ratio")
    dates = np.array(contact_results.days)
    probs = np.array(contact_results.probabilities)
    stdevs = np.array(contact_results.stdevs)
    ax0.fill_between(dates, probs + 3*stdevs, probs - 3*stdevs, facecolor="lightblue", alpha=0.5)
    ax0.plot(dates, probs, color="lightblue", linewidth=4, label="Estimated contact ratio")

    plt_i = infected_history[state].get_plottable(dates.min(), dates.max())
    ax1 = ax0.twinx()
    ax1.plot(plt_i.dates, plt_i.infections, color="deeppink", linewidth=4, label="Daily Infections")

    h0, l0 = ax0.get_legend_handles_labels()
    h1, l1 = ax1.get_legend_handles_labels()
    ax0.legend(h0+h1, l0+l1, loc="upper left")

    fig.show()


if __name__ == '__main__':
    main()
