"""
    First simulation of the system
"""
import os
import pickle
from datetime import date as Date
from datetime import timedelta as TimeDelta

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes

import settings
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories)
from sim import ProgramInput, Simulator, default_world_properties, ContactSearchResult


def main():
    input_data = ProgramInput(output_file=settings.default_output_file,
                              state="PA",
                              world_properties=default_world_properties(),
                              start_day=Date(2021, 8, 7),
                              end_day=Date(2022, 1, 1),
                              state_info=load_state_info(),
                              population_scale=50,
                              contact_prob=1.75,
                              run_count=50)

    plot_start = Date(2021, 7, 1)  # input_data.start_day - TimeDelta(days=5)
    # plot_start = input_data.start_day
    plot_end = input_data.end_day

    input_data.test_history = load_state_histories()
    input_data.vax_history = load_vaccine_histories()
    input_data.infected_history = load_state_estimates()
    input_data.variant_history = load_variant_history()


    fig: Figure = plt.figure(figsize=(12, 12))

    # ax0 = fig.subplots(1)
    ax0, ax1 = fig.subplots(2)
    ax0: Axes
    ax1: Axes
    ax0.set_title(
        f"Simulated {input_data.state} against Covidestim.org Infection Estimates")
    ax0.set_xlabel("Date")
    ax0.set_ylabel(f"Infected People (Pop={input_data.state_info[input_data.state].population / 1e6:0.1f}M)")

    plt_i = input_data.infected_history[input_data.state].get_plottable(plot_start, plot_end)

    ax0.plot(plt_i.dates, plt_i.total_infections, "deeppink", linewidth=3, label="Covidestim.org Infections")
    ax1.plot(plt_i.dates, plt_i.infections, "coral", linewidth=3, label="Covidestim.org Daily Infections")
    ax1.set_title("Daily Infections")
    ax1.set_ylabel("People")

    probabilities = [1., 1.5, 1.75, 2.0 ]
    cmap = cm.get_cmap("cool")

    for prob in probabilities:
        scaled = (prob - min(probabilities)) / (max(probabilities) - min(probabilities))
        color = cmap(scaled)
        input_data.contact_prob = prob
        simulator = Simulator(input_data, settings.default_input_file)
        result = simulator.run(use_cache=True)
        print(f"took {result.run_time:0.2f}s to run")
        plt_r = result.get_plottable(input_data.state, plot_start, plot_end)

        ax0.fill_between(plt_r.dates, plt_r.total_infections.upper, plt_r.total_infections.lower, facecolor=color,
                         alpha=0.5)

        ax0.plot(plt_r.dates, plt_r.total_infections.mean, color=color, linewidth=2,
                 label=f"Simulated Total Infections ({prob:0.2f})")

        ax1.fill_between(plt_r.dates, plt_r.new_infections.lower, plt_r.new_infections.upper, facecolor=color, alpha=0.5)
        ax1.plot(plt_r.dates, plt_r.new_infections.mean, color=color, linewidth=2, label=f"Simulated Daily Infections ({prob:0.2f})")

    # input_data.variant_history = None
    # simulator = Simulator(input_data, _input_path)
    # result = simulator.run()
    # plt_nd = result.get_plottable(input_data.state, plot_start, plot_end)
    #
    # ax0.fill_between(plt_nd.dates, plt_nd.infections.upper, plt_nd.infections.lower, facecolor="green", alpha=0.3)
    # ax0.plot(plt_nd.dates, plt_nd.infections.mean, "darkgreen", linewidth=2, label="Simulated Infections w/o Delta")

    # # ax1 = ax.twinx()
    # plt_h = input_data.test_history[input_data.state].get_plottable(plot_start, plot_end)
    # ax1.plot(plt_i.dates, plt_i.total_cases, "darkred", linewidth=3)
    # ax1.plot(plt_h.dates, plt_h.positive_tests, "red", marker="+", linewidth=0)

    ax1.legend()
    ax0.legend()
    fig.show()


if __name__ == '__main__':
    main()
