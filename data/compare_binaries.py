import os
import pickle
import settings

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import Figure
from matplotlib.axes._axes import Axes

from history import (load_state_estimates, load_state_histories, load_state_info, load_variant_history,
                     load_vaccine_histories)
from sim import ProgramInput, Simulator, default_world_properties, ContactSearchResult, SimulationResult


def main():

    with open(os.path.join(settings.cache_folder, "new_bin.pickle"), "rb") as handle:
        new_results: SimulationResult = pickle.load(handle)

    with open(os.path.join(settings.cache_folder, "old_bin.pickle"), "rb") as handle:
        old_results: SimulationResult = pickle.load(handle)

    n_ = new_results.get_plottable("FL")
    o_ = old_results.get_plottable("FL")

    fig: Figure = plt.figure(figsize=(16, 10))
    ax: Axes = fig.subplots()
    ax.plot(o_.dates, o_.total_vaccinated.mean, color="lightblue", linewidth=4, label="Old Binary")
    ax.plot(n_.dates, n_.total_vaccinated.mean, color="lightgreen", linewidth=4, label="New Binary")
    ax.legend()

    fig.show()




if __name__ == '__main__':
    main()
