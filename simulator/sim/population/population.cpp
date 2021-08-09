#include "population.hpp"
#include <cmath>

sim::Population::Population(int unscaled_size, int scale) {
    scale_ = scale;
    long scaled_population = static_cast<long>(std::round(static_cast<double>(unscaled_size) / scale));
    for (long i = 0; i < scaled_population; ++i) {
        people.emplace_back();
    }
}

void sim::Population::Reset() {
    for (auto &p : people) {
        p.Reset();
    }

    today = 0;
    vaccine_saves = 0;
    natural_saves = 0;
    total_infections = 0;
    total_vaccinated = 0;
    never_infected = 0;
    total_delta_infections = 0;
    total_alpha_infections = 0;
    reinfections =0;
    vaccinated_infections = 0;
}
