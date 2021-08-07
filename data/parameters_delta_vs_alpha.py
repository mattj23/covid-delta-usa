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
    state = "IL"
    start_date = Date(2020, 9, 1)
    end_date = Date(2020, 10, 15)
    plot_start = Date(2020, 8, 15)  # input_data.start_day - TimeDelta(days=5)
    # plot_start = input_data.start_day
    plot_end = end_date

    state_info = load_state_info()
    infected_history = load_state_estimates()

    # Create some world properties
    properties = default_world_properties()

    fig: Figure = plt.figure(figsize=(12, 6*3))
    fig.subplots_adjust(top=0.95, bottom=0.05)
    ax0, ax1, ax2 = fig.subplots(3)

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
    ax2.set_title("Infection Curve")
    properties.alpha.infectivity.plot(ax2, color="darkorange", linewidth=3, label="Ashcroft et al.")
    ax2.set_xlabel("Days from symptom onset")
    ax2.set_ylabel("Relative infectivity")


    # Screw with the infection curve
    shape, rate, shift = 97.18750, 3.71875, 25.6250
    # c1 = custom_infectivity_curve(97, 1.71875, 50) # really lagged
    c1 = custom_infectivity_curve(97, 2.71875, 33)
    c1.plot(ax2, "lightblue", 3, "Lagged curve")

    bundles = [
        dict(curve=properties.alpha.infectivity, label="Ashcroft et al", color="darkorange", contact=1.8),
        # dict(curve=properties.alpha.infectivity, label="Ashcroft et al", color="lightgreen", contact=2.2),
        # dict(curve=c1, label="Lagged curve", color="lightblue", contact=1.15)
               ]

    for bundle in bundles:
        properties.alpha.infectivity = bundle["curve"]
        input_data = ProgramInput(
            output_file=settings.default_output_file,
            state=state,
            world_properties=properties,
            start_day=start_date,
            end_day=end_date,
            contact_prob=bundle['contact'],
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

        label = f"{bundle['label']}, p={bundle['contact']}"

        plt_r = result.get_plottable(input_data.state, plot_start, plot_end)
        ax0.fill_between(plt_r.dates, plt_r.total_infections.upper, plt_r.total_infections.lower, facecolor=bundle['color'], alpha=0.5)
        ax0.plot(plt_r.dates, plt_r.total_infections.mean, bundle['color'], linewidth=2, label=label)
        ax1.plot(plt_r.dates, plt_r.new_infections.mean, bundle['color'], linewidth=2, label=label)
        ax1.fill_between(plt_r.dates, plt_r.new_infections.lower, plt_r.new_infections.upper, facecolor=bundle['color'], alpha=0.5)

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
