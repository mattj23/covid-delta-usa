"""
    This script uses the simulator to try to estimate the relative difference in infectivity between the
    alpha/ancestral strains and the delta variant.

    bin dates end on 5/8, 5/22, 6/5, 6/19, 7/3, 7/17, and 7/31
"""
from typing import List, Dict, Optional

import numpy as np
from datetime import date as Date
from datetime import timedelta as TimeDelta

import settings
from sim.dynamics import VariantProperties, DiscreteFunction, WorldProperties
from sim import default_world_properties, ProgramInput, Simulator
from sim.program_input import ProgramMode, ProgramOptions
from sim.world_defaults import custom_infectivity_curve
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories, StateEstimates, DailyEstimate, StateInfo)

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes
import scipy.stats as stats

from find_contact_prob import find_contact_prob


def main():
    state = "TN"
    start_date = Date(2021, 6, 20)  # Corresponds with a collection bin
    end_date = Date(2021, 8, 1)
    plot_start = Date(2021, 5, 15)  # input_data.start_day - TimeDelta(days=5)
    # plot_start = input_data.start_day
    plot_end = end_date

    state_info = load_state_info()
    infected_history = load_state_estimates()

    # Create some world properties
    properties = default_world_properties()

    fig: Figure = plt.figure(figsize=(12, 6 * 4))
    fig.subplots_adjust(top=0.95, bottom=0.05)
    ax0, ax1, ax2 = fig.subplots(3)

    ax0: Axes
    ax0.set_title(f"Simulated {state} against Covidestim.org Infection Estimates")
    ax0.set_xlabel("Date")
    ax0.set_ylabel(f"Infected People (Pop={state_info[state].population / 1e6:0.1f}M)")

    plt_i = infected_history[state].get_plottable(plot_start, plot_end)
    ax0.plot(plt_i.dates, plt_i.total_infections, "deeppink", linewidth=5, label="Covidestim.org Infections")

    ax1: Axes
    ax1.set_title("Daily Infections")
    ax1.set_ylabel("People")
    ax1.plot(plt_i.dates, plt_i.infections, "deeppink", linewidth=5, label="Covidestim.org Daily Infections")

    ax2: Axes
    ax2.set_title("Proportion of Delta variant in region, 2-wk new case bins")
    ax2.set_xlabel("End date of two week proportion bin")
    ax2.set_ylabel("Proportion of B.1.617.2 in new cases")
    ax2.set_ylim(-.05, 1.1)

    all_variant_history = load_variant_history()
    variant_history = all_variant_history[state]
    x_, y_ = zip(*sorted((v["date"], v["variants"]["delta"]) for v in variant_history))
    ax2.plot(x_, y_, "deeppink", linewidth=4, label="CDC Variant Tracking", marker="o")

    test_cases = np.arange(2.0, 4.0, 0.25)
    test_colors = (test_cases - min(test_cases)) / (max(test_cases) - min(test_cases))
    colors = cm.get_cmap("jet")(test_colors)

    for d_scale, color in zip(test_cases, colors):
        delta_scale = d_scale
        properties.delta.infectivity = properties.alpha.infectivity.scale_y(delta_scale)
        state_info = load_state_info()
        vax_history = load_vaccine_histories()
        variant_history = load_variant_history()
        infected_history = load_state_estimates()

        contact_prob = find_contact_prob(state, start_date, properties, state_info=state_info,
                                         vax_history=vax_history, variant_history=variant_history,
                                         infected_history=infected_history)

        input_data = ProgramInput(output_file=settings.default_output_file,
                                  options=ProgramOptions(full_history=True,
                                                         expensive_stats=False,
                                                         mode=ProgramMode.Simulate),
                                  state=state,
                                  world_properties=properties,
                                  start_day=start_date,
                                  end_day=end_date,
                                  contact_prob=contact_prob,
                                  state_info=state_info,
                                  population_scale=25,
                                  vax_history=vax_history,
                                  variant_history=variant_history,
                                  infected_history=infected_history,
                                  run_count=100)
        simulator = Simulator(input_data, settings.default_input_file)
        result = simulator.run()
        print(f"took {result.run_time:0.2f}s to run")

        plt_r = result.get_plottable(input_data.state, plot_start, plot_end)
        for date, value in zip(plt_r.dates, plt_r.new_infections.mean):
            if value == 0:
                print(f"Broken date: {date}")

        # Cumulative infections
        # ax0.fill_between(plt_r.dates, plt_r.total_infections.upper, plt_r.total_infections.lower, facecolor=color,
        #                  alpha=0.5)
        ax0.plot(plt_r.dates, plt_r.total_infections.mean, color=color, linewidth=2, label=f"{delta_scale:0.2f}x, {contact_prob:0.2f} contact")

        # Daily infections
        ax1.plot(plt_r.dates, plt_r.new_infections.mean, color=color, linewidth=2, label=f"{delta_scale:0.2f}x, {contact_prob:0.2f} contact")
        # ax1.fill_between(plt_r.dates, plt_r.new_infections.lower, plt_r.new_infections.upper, facecolor=color,
        #                  alpha=0.5)
        # ax1.plot(plt_r.dates, plt_r.new_delta_infections.mean, "cyan", linewidth=2, label="Daily delta infections")
        # ax1.fill_between(plt_r.dates, plt_r.new_delta_infections.mean, 0, facecolor="cyan", alpha=0.5)

        # Variant proportions
        new_delta = zip(plt_r.dates, plt_r.new_delta_infections.mean)
        new_all = zip(plt_r.dates, plt_r.new_infections.mean)
        bins = {k['date']: {"delta": 0, "total": 0} for k in variant_history[state]}
        bin_dates = sorted(bins.keys())
        for sim_date, cases in new_delta:
            for bin_date in bin_dates:
                if sim_date <= bin_date:
                    bins[bin_date]["delta"] += cases
                    break
        for sim_date, cases in new_all:
            for bin_date in bin_dates:
                if sim_date <= bin_date:
                    bins[bin_date]["total"] += cases
                    break

        ratio_x, ratio_y = zip(*sorted((k, v["delta"] / v["total"]) for k, v in bins.items() if v["total"]))
        ax2.plot(ratio_x, ratio_y, color=color, linewidth=2, label=f"Simulated {delta_scale:0.2f}x Infectivity")
        # ratio = plt_r.new_delta_infections / plt_r.new_infections
        # ax2.plot(plt_r.dates, ratio.mean, "mediumpurple", label="New infections delta ratio")

    ax2.plot([start_date, start_date], [0, 1], "lightcoral", linestyle="--")
    ax0.legend()
    ax1.legend()
    ax2.legend()
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
