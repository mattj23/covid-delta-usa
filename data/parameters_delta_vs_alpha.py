"""
    This script uses the simulator to try to estimate the relative difference in infectivity between the
    alpha/ancestral strains and the delta variant.
"""
from typing import List, Dict, Optional

import numpy
import numpy as np
from datetime import date as Date
from datetime import timedelta as TimeDelta

import settings
from states import state_abbrevs
from sim.dynamics import VariantProperties, DiscreteFunction, WorldProperties
from history.estimates import StateEstimates, DailyEstimate
from history.states import StateInfo
from sim import default_world_properties, ProgramInput, Simulator
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories)

import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes


def main():
    # Create some world properties
    properties = default_world_properties()

    # Accumulate the whole USA into one population
    state_info = load_state_info()
    estimates = load_state_estimates()

    infected_history = StateEstimates("USA", {})
    all_dates = set()
    for state in state_abbrevs():
        for k in estimates[state].estimates.keys():
            all_dates.add(k)

    working_date = min(all_dates)
    latest_date = max(all_dates)
    while working_date <= latest_date:
        daily = DailyEstimate(0, 0, 0, 0, 0, 0)
        for state in state_abbrevs():
            if working_date not in estimates[state].estimates:
                continue

            daily.cases += estimates[state].estimates[working_date].cases
            daily.diagnoses += estimates[state].estimates[working_date].diagnoses
            daily.infections += estimates[state].estimates[working_date].infections
            daily.total_cases += estimates[state].estimates[working_date].total_cases
            daily.total_diagnoses += estimates[state].estimates[working_date].total_diagnoses
            daily.total_infections += estimates[state].estimates[working_date].total_infections

        infected_history.estimates[working_date] = daily
        working_date += TimeDelta(days=1)

    usa_info = StateInfo("USA", 0, [])
    for state in state_abbrevs():
        usa_info.population += state_info[state].population

    input_data = ProgramInput(
        output_file=settings.default_output_file,
        state="USA",
        world_properties=properties,
        start_day=Date(2020, 11, 1),
        end_day=Date(2021, 8, 1),
        contact_prob=2.16,
        state_info={"USA": usa_info},
        population_scale=1000,
        vax_history=load_vaccine_histories(),
        variant_history=load_variant_history(),
        infected_history={"USA": infected_history},
        run_count=5
    )

    simulator = Simulator(input_data, settings.default_input_file)
    result = simulator.run()

    plot_start = Date(2020, 6, 1)  # input_data.start_day - TimeDelta(days=5)
    # plot_start = input_data.start_day
    plot_end = input_data.end_day

    print(f"took {result.run_time:0.2f}s to run")

    fig: Figure = plt.figure(figsize=(12, 12))

    ax0, ax1 = fig.subplots(2)
    ax0: Axes
    ax1: Axes
    ax0.set_title(
        f"Simulated {input_data.state} against Covidestim.org Infection Estimates ({input_data.contact_prob:0.2f} contact probability)")
    ax0.set_xlabel("Date")
    ax0.set_ylabel(f"Infected People (Pop={input_data.state_info[input_data.state].population / 1e6:0.1f}M)")

    plt_r = result.get_plottable(input_data.state, plot_start, plot_end)
    plt_i = input_data.infected_history[input_data.state].get_plottable(plot_start, plot_end)

    ax0.fill_between(plt_r.dates, plt_r.total_infections.upper, plt_r.total_infections.lower, facecolor="orange", alpha=0.5)
    ax0.plot(plt_i.dates, plt_i.total_infections, "deeppink", linewidth=3, label="Covidestim.org Infections")
    ax0.plot(plt_r.dates, plt_r.total_infections.mean, "darkorange", linewidth=2, label="Simulated Total Infections w/ Delta")

    ax1.fill_between(plt_r.dates, plt_r.new_infections.lower, plt_r.new_infections.upper, facecolor="gold", alpha=0.5)
    ax1.plot(plt_i.dates, plt_i.infections, "coral", linewidth=3, label="Covidestim.org Daily Infections")
    ax1.plot(plt_r.dates, plt_r.new_infections.mean, "goldenrod", linewidth=2, label="Simulated Daily Infections")
    ax1.set_title("Daily Infections")
    ax1.set_ylabel("People")
    ax1.legend()
    fig.show()


def plot_world_properties(properties: WorldProperties):
    count = 4
    fig: Figure = plt.figure(figsize=(12, 8 * count))
    axes = fig.subplots(count)
    fig.subplots_adjust(top=0.95, bottom=0.05)

    # Plotting of the incubation curves
    # =======================================================================================
    def plot_incubation(ax: Axes):
        ax.set_title("Incubation Distribution, Alpha/Pre-Alpha vs Delta")
        ax.set_xlabel("Days from infection event")
        ax.set_ylabel("Probability of symptom onset")
        width = 0.4

        def _plot_bars(v: VariantProperties, shift, color, label):
            x_, y_ = zip(*enumerate(v.incubation))
            x_ = np.array(x_[1:])
            y_ = np.diff(np.array(y_), axis=0)
            ax.bar(x_ + shift, y_, width=width, facecolor=color, label=label)
            # f = scipy.interpolate.interp1d(x_, y_, kind="quadratic")
            # x_l = numpy.linspace(min(x_), max(x_), 100)
            # ax0.plot(x_l, f(x_l), color=color, linewidth=3, linestyle="--")

        _plot_bars(properties.alpha, -width / 2, "orange", "Alpha/Pre-Alpha")
        _plot_bars(properties.delta, width / 2, "red", "Delta")

        ax.legend()

    # Plotting of the infectivity curves
    # =======================================================================================
    def _plot_discrete_function(ax: Axes, f: DiscreteFunction, color, width, label):
        x_, y_ = [np.array(v) for v in zip(*enumerate(f.values))]
        x_ = x_ - f.offset
        ax.plot(x_, y_, color=color, linewidth=width, label=label)

    def plot_infectivity(ax: Axes):
        ax.set_title("Relative Infectivity, Alpha/Pre-Alpha vs Delta")
        ax.set_xlabel("Days from symptom onset")
        ax.set_ylabel("Relative infectiousness")

        _plot_discrete_function(ax, properties.alpha.infectivity, "orange", 2, "Alpha/Pre-Alpha")
        _plot_discrete_function(ax, properties.delta.infectivity, "red", 2, "Delta")
        ax.legend()

    # Plotting of the vaccine immunity curves
    # =======================================================================================
    def plot_vaccine_immunity(ax: Axes):
        ax.set_title("Vaccine Efficacy, Alpha/Pre-Alpha vs Delta")
        ax.set_xlabel("Days from first shot (assuming 2nd is at day 21)")
        ax.set_ylabel("Efficacy (Reduction in risk ratio)")

        _plot_discrete_function(ax, properties.alpha.vax_immunity, "orange", 2, "Alpha/Pre-Alpha")
        _plot_discrete_function(ax, properties.delta.vax_immunity, "red", 2, "Delta")
        ax.legend()

    # Plotting of the natural immunity curves
    # =======================================================================================
    def plot_natural_immunity(ax: Axes):
        ax.set_title("Efficacy of Natural Immunity, Alpha/Pre-Alpha vs Delta")
        ax.set_xlabel("Days from infection")
        ax.set_ylabel("Efficacy (Reduction in risk ratio)")
        ax.set_ylim(0, 1.2)

        _plot_discrete_function(ax, properties.alpha.natural_immunity, "orange", 2, "Alpha/Pre-Alpha")
        _plot_discrete_function(ax, properties.delta.natural_immunity, "red", 2, "Delta")
        ax.legend()

    # Create Plots
    plot_incubation(axes[0])
    plot_infectivity(axes[1])
    plot_natural_immunity(axes[2])
    plot_vaccine_immunity(axes[3])

    fig.show()


if __name__ == '__main__':
    main()
