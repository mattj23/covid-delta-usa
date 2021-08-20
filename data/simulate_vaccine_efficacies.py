"""
    First simulation of the system
"""
from datetime import date as Date

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes

import settings
from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories)
from find_contact_prob import find_days_contact_prob
from sim import ProgramInput, Simulator, default_world_properties, ContactSearchResult
from sim.world_defaults import (pfizer_delta_efficacy_lopez_uk, pfizer_delta_efficacy_israeli_moh,
                                pfizer_alpha_efficacy)


def main():
    state = "PA"
    input_data = ProgramInput(output_file=settings.default_output_file,
                              state=state,
                              world_properties=default_world_properties(),
                              start_day=Date(2021, 8, 7),
                              end_day=Date(2022, 1, 1),
                              state_info=load_state_info(),
                              population_scale=50,
                              contact_prob=1.75,
                              run_count=50)

    plot_start = Date(2021, 5, 1)  # input_data.start_day - TimeDelta(days=5)
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
    ax0.set_ylabel(f"Total Infected People (Pop={input_data.state_info[input_data.state].population / 1e6:0.1f}M)")

    plt_i = input_data.infected_history[input_data.state].get_plottable(plot_start, plot_end)

    ax0.plot(plt_i.dates, plt_i.total_infections, "deeppink", linewidth=3, label="Covidestim.org Infections")
    ax1.plot(plt_i.dates, plt_i.infections, "coral", linewidth=3, label="Covidestim.org Daily Infections")
    ax1.set_title("Daily Infections")
    ax1.set_ylabel("People Infected Per Day")

    def simulate(efficacy, label, color):
        input_data.world_properties.delta.vax_immunity = efficacy
        contact_prob = find_days_contact_prob(state, input_data.start_day,
                                              input_data.world_properties,
                                              input_data.state_info,
                                              input_data.vax_history,
                                              input_data.variant_history,
                                              input_data.infected_history)
        input_data.contact_prob = contact_prob

        simulator = Simulator(input_data, settings.default_input_file)
        result = simulator.run(use_cache=True)
        print(f"took {result.run_time:0.2f}s to run")
        plt_r = result.get_plottable(input_data.state, plot_start, plot_end)

        ax0.fill_between(plt_r.dates, plt_r.total_infections.upper, plt_r.total_infections.lower, facecolor=color,
                         alpha=0.5)

        ax0.plot(plt_r.dates, plt_r.total_infections.mean, color=color, linewidth=2,
                 label=f"Simulated Total ({label}, cp={contact_prob:0.2f})")

        ax1.fill_between(plt_r.dates, plt_r.new_infections.lower, plt_r.new_infections.upper, facecolor=color,
                         alpha=0.5)
        ax1.plot(plt_r.dates, plt_r.new_infections.mean, color=color, linewidth=2,
                 label=f"Simulated Daily ({label}, cp={contact_prob:0.2f})")

    # simulate(pfizer_alpha_efficacy(), "Same as Alpha", "seagreen")
    simulate(pfizer_delta_efficacy_lopez_uk(), "Lopez Bernal et al.", "gold")
    simulate(pfizer_delta_efficacy_israeli_moh(), "Israeli MOH", "lightcoral")

    ax1.legend(loc="upper left")
    ax0.legend(loc="upper left")
    fig.show()


if __name__ == '__main__':
    main()
