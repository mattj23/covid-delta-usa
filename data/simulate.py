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


def main():
    cache_file = os.path.join(settings.cache_folder, "temp.pickle")
    with open(cache_file, "rb") as handle:
        contact_results: ContactSearchResult = pickle.load(handle)

    input_data = ProgramInput(output_file=settings.default_output_file,
                              state="FL",
                              world_properties=default_world_properties(),
                              start_day=Date(2021, 8, 7),
                              end_day=Date(2021, 12, 1),
                              state_info=load_state_info(),
                              population_scale=100,
                              contact_prob=1,
                              run_count=200)

    plot_start = Date(2020, 6, 1)  # input_data.start_day - TimeDelta(days=5)
    # plot_start = input_data.start_day
    plot_end = input_data.end_day

    input_data.test_history = load_state_histories()
    input_data.vax_history = load_vaccine_histories()
    input_data.infected_history = load_state_estimates()
    input_data.variant_history = load_variant_history()

    fig: Figure = plt.figure(figsize=(16, 10))

    # ax0 = fig.subplots(1)
    ax0 = fig.subplots(1)
    ax0: Axes
    ax0.set_title(f"Contact Ratio for {input_data.state}, Delta 2.5x modifier")
    ax0.set_xlabel("Date")
    ax0.set_ylabel(f"Contact ratio")
    dates = np.array(contact_results.days)
    probs = np.array(contact_results.probabilities)
    stdevs = np.array(contact_results.stdevs)
    ax0.fill_between(dates, probs + 3*stdevs, probs - 3*stdevs, facecolor="lightblue", alpha=0.5)
    ax0.plot(dates, probs, color="lightblue", linewidth=4, label="Estimated contact ratio")

    plt_i = input_data.infected_history[input_data.state].get_plottable(dates.min(), dates.max())
    ax1 = ax0.twinx()
    ax1.plot(plt_i.dates, plt_i.infections, color="deeppink", linewidth=4, label="Daily Infections")

    h0, l0 = ax0.get_legend_handles_labels()
    h1, l1 = ax1.get_legend_handles_labels()
    ax0.legend(h0+h1, l0+l1, loc="upper left")

    fig.show()
    return


    # Vaccination counterfactual of 2x actual
    # state_pop = input_data.state_info[input_data.state].population
    # for date, record in input_data.vax_history[input_data.state].records.items():
    #     record.completed_vax = min(state_pop, record.completed_vax * 2)

    simulator = Simulator(input_data, settings.default_input_file)
    result = simulator.find_contact_prob()


    return
    result = simulator.run()

    print(f"took {result.run_time:0.2f}s to run")

    fig: Figure = plt.figure(figsize=(12, 12))

    # ax0 = fig.subplots(1)
    ax0, ax1 = fig.subplots(2)
    ax0: Axes
    ax1: Axes
    ax0.set_title(
        f"Simulated {input_data.state} against Covidestim.org Infection Estimates ({input_data.contact_prob:0.2f} contact probability)")
    ax0.set_xlabel("Date")
    ax0.set_ylabel(f"Infected People (Pop={input_data.state_info[input_data.state].population / 1e6:0.1f}M)")

    plt_r = result.get_plottable(input_data.state, plot_start, plot_end)
    plt_i = input_data.infected_history[input_data.state].get_plottable(plot_start, plot_end)

    ax0.fill_between(plt_r.dates, plt_r.total_infections.upper, plt_r.total_infections.lower, facecolor="orange",
                     alpha=0.5)
    ax0.plot(plt_i.dates, plt_i.total_infections, "deeppink", linewidth=3, label="Covidestim.org Infections")
    ax0.plot(plt_r.dates, plt_r.total_infections.mean, "darkorange", linewidth=2,
             label="Simulated Total Infections w/ Delta")

    ax1.fill_between(plt_r.dates, plt_r.new_infections.lower, plt_r.new_infections.upper, facecolor="gold", alpha=0.5)
    ax1.plot(plt_i.dates, plt_i.infections, "coral", linewidth=3, label="Covidestim.org Daily Infections")
    ax1.plot(plt_r.dates, plt_r.new_infections.mean, "goldenrod", linewidth=2, label="Simulated Daily Infections")
    ax1.set_title("Daily Infections")
    ax1.set_ylabel("People")
    ax1.legend()

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

    ax0.legend()
    fig.show()


if __name__ == '__main__':
    main()
