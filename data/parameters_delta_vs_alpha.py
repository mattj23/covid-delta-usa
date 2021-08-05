"""
    This script uses the simulator to try to estimate the relative difference in infectivity between the
    alpha/ancestral strains and the delta variant.
"""
from typing import List

import numpy
import numpy as np

import scipy.interpolate
import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes
from sim.dynamics import VariantProperties, DiscreteFunction, WorldProperties
from sim import default_world_properties, ProgramInput, Simulator
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories)


def main():
    # Create some world properties
    properties = default_world_properties()

    # Accumulate the whole USA into one population
    estimates = load_state_estimates()


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
