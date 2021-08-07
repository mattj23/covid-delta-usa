"""
    This script uses the simulator to try to estimate the relative difference in infectivity between the
    alpha/ancestral strains and the delta variant.
"""
from typing import List, Dict, Optional

import numpy as np
from datetime import date as Date
from datetime import timedelta as TimeDelta

import settings
from sim.dynamics import VariantProperties, DiscreteFunction, WorldProperties
from sim import default_world_properties, ProgramInput, Simulator
from sim.world_defaults import custom_infectivity_curve
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories, StateEstimates, DailyEstimate, StateInfo)

import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes


def main():
    state = "WA"
    start_date = Date(2021, 6, 20) # Corresponds with a collection bin
    end_date = Date(2021, 7, 17)
    plot_start = Date(2021, 5, 15)  # input_data.start_day - TimeDelta(days=5)
    # plot_start = input_data.start_day
    plot_end = end_date

    state_info = load_state_info()
    infected_history = load_state_estimates()

    # Create some world properties
    properties = default_world_properties()

    fig: Figure = plt.figure(figsize=(12, 6*4))
    fig.subplots_adjust(top=0.95, bottom=0.05)
    ax0, ax1, ax2, ax3 = fig.subplots(4)

    ax0: Axes
    ax0.set_title(f"Simulated {state} against Covidestim.org Infection Estimates")
    ax0.set_xlabel("Date")
    ax0.set_ylabel(f"Infected People (Pop={state_info[state].population / 1e6:0.1f}M)")

    plt_i = infected_history[state].get_plottable(plot_start, plot_end)
    ax0.plot(plt_i.dates, plt_i.total_infections, "deeppink", linewidth=3, label="Covidestim.org Infections")

    ax1: Axes
    ax1.set_title("Daily Infections")
    ax1.set_ylabel("People")
    ax1.plot(plt_i.dates, plt_i.infections, "coral", linewidth=3, label="Covidestim.org Daily Infections")

    ax2: Axes
    ax2.set_title("Proportion of Delta variant in region, 2-wk new case bins")
    ax2.set_xlabel("End date of two week proportion bin")
    ax2.set_ylabel("Proportion of B.1.617.2 in new cases")
    ax2.set_ylim(-.05, 1.1)

    all_variant_history = load_variant_history()
    variant_history = all_variant_history[state]
    x_, y_ = zip(*sorted((v["date"], v["variants"]["delta"]) for v in variant_history))
    ax2.plot(x_, y_, "deeppink", linewidth=4, label="CDC Variant Tracking", marker="o")

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

    def plot_vaccine_immunity(ax: Axes):
        ax.set_title("Vaccine Efficacy, Alpha/Pre-Alpha vs Delta")
        ax.set_xlabel("Days from first shot (assuming 2nd is at day 21)")
        ax.set_ylabel("Efficacy (Reduction in risk ratio)")

        _plot_discrete_function(ax, properties.alpha.vax_immunity, "orange", 2, "Alpha/Pre-Alpha")
        _plot_discrete_function(ax, properties.delta.vax_immunity, "red", 2, "Delta")
        ax.legend()


    bundles = [
        dict(factor=2.3),
               ]

    # properties.delta.vax_immunity = properties.alpha.vax_immunity.scale_y(0.8)
    plot_vaccine_immunity(ax3)

    for bundle in bundles:
        delta_scale = bundle['factor']
        properties.delta.infectivity = properties.alpha.infectivity.scale_y(delta_scale)
        input_data = ProgramInput(
            output_file=settings.default_output_file,
            state=state,
            world_properties=properties,
            start_day=start_date,
            end_day=end_date,
            contact_prob=1.8,
            state_info=load_state_info(),
            population_scale=50,
            vax_history=load_vaccine_histories(),
            variant_history=load_variant_history(),
            infected_history=load_state_estimates(),
            run_count=20
        )

        simulator = Simulator(input_data, settings.default_input_file)
        result = simulator.run()

        print(f"took {result.run_time:0.2f}s to run")

        label = f""
        color = "darkorange"

        plt_r = result.get_plottable(input_data.state, plot_start, plot_end)

        new_delta = zip(plt_r.dates, plt_r.new_delta_infections.mean)
        new_all = zip(plt_r.dates, plt_r.new_infections.mean)
        bins = {k['date']: {"delta": 0, "total": 0} for k in variant_history if k['date'] > start_date and k['date'] <= end_date}
        dates = sorted(bins.keys())
        for date, cases in new_delta:
            for d in dates:
                if date <= d:
                    bins[d]["delta"] += cases
        for date, cases in new_all:
            for d in dates:
                if date <= d:
                    bins[d]["total"] += cases

        ratio_x, ratio_y = zip(*sorted((k, v["delta"] / v["total"]) for k, v in bins.items()))
        ax2.plot(ratio_x, ratio_y, "teal", linewidth=2, label=f"Simulated {delta_scale:0.2f}x Infectivity", marker="x")




        ax0.fill_between(plt_r.dates, plt_r.total_infections.upper, plt_r.total_infections.lower, facecolor=color, alpha=0.5)
        ax0.plot(plt_r.dates, plt_r.total_infections.mean, color, linewidth=2, label=label)

        ax1.plot(plt_r.dates, plt_r.new_infections.mean, color, linewidth=2, label=label)
        ax1.fill_between(plt_r.dates, plt_r.new_infections.lower, plt_r.new_infections.upper, facecolor=color, alpha=0.5)

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
