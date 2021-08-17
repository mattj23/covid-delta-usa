#include "population.hpp"
#include "../timer.hpp"
#include <cmath>

sim::Population::Population(int unscaled_size, int scale) {
    scale_ = scale;
    long scaled_population = static_cast<long>(std::round(static_cast<double>(unscaled_size) / scale));
    for (long i = 0; i < scaled_population; ++i) {
        people.emplace_back();
        infectious.emplace_back(i);
    }
}

void sim::Population::Reset() {
    today = 0;
    vaccine_saves = 0;
    natural_saves = 0;
    total_infections = 0;
    total_vaccinated = 0;
    never_infected = static_cast<int>(people.size());
    total_delta_infections = 0;
    total_alpha_infections = 0;
    reinfections = 0;
    vaccinated_infections = 0;

    infectious_ptr_ = 0;
    infectious_indices.clear();
    infectious.clear();
    unvaxxed_indices.clear();

    for (auto &p : people)
        p.Reset();

    for (size_t i = 0; i < people.size(); ++i)
        unvaxxed_indices.push_back(i);
}

void sim::Population::CopyFrom(const Population &other) {
    if (people.size() != other.people.size()) {
        // Throw exception
    }

    today = other.today;
    vaccine_saves = other.vaccine_saves;
    natural_saves = other.natural_saves;
    total_infections = other.total_infections;
    total_vaccinated = other.total_vaccinated;
    never_infected = other.never_infected;
    total_delta_infections = other.total_delta_infections;
    total_alpha_infections = other.total_alpha_infections;
    reinfections = other.reinfections;
    vaccinated_infections = other.vaccinated_infections;
    infectious_ptr_ = other.infectious_ptr_;

    // These are stl container assignments, which should efficiently make deep copies
    // of the container and its contents
    people = other.people;
    infectious = other.infectious;
    infectious_indices = other.infectious_indices;
    unvaxxed_indices = other.unvaxxed_indices;
}
