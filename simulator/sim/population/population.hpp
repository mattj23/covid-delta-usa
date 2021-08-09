#pragma once

#include <vector>
#include <unordered_set>
#include "person.hpp"

namespace sim {

    /** @class Population
     *
     * @brief Is a data-only representation of a population of individuals at a given time.
     */
    class Population {
    public:
        /** @brief Creates a population from an unscaled number of individuals and a scale factor
         *
         * @param unscaled_size the number of people in the population *before* scaling
         * @param scale an integer that defines how many people in the real population are represented by each
         * simulated individual, also can be thought of as the model being built to a 1:scale scale
         */
        Population(int unscaled_size, int scale);

        void Reset();
        void CopyFrom(const Population& other);

        [[nodiscard]] inline int Scale() const { return scale_; }

        std::vector<Person> people;
        std::unordered_set<size_t> infectious_indices;
        std::vector<size_t> unvaxxed_indices;

        inline int TotalInfections() const { return total_infections * scale_; }
        inline int TotalVaccinated() const { return total_vaccinated * scale_; }
        inline int NeverInfected() const { return never_infected * scale_; }
        inline int TotalDeltaInfections() const { return total_delta_infections * scale_; }
        inline int TotalAlphaInfections() const { return total_alpha_infections * scale_; }
        inline int Reinfections() const { return reinfections * scale_; }
        inline int VaccineSaves() const { return vaccine_saves * scale_; }
        inline int NaturalSaves() const { return natural_saves * scale_; }

        int today{};
        int vaccine_saves{};
        int natural_saves{};

        int total_infections{};
        int total_vaccinated{};
        int never_infected{};
        int total_delta_infections{};
        int total_alpha_infections{};
        int reinfections{};
        int vaccinated_infections{};
    private:
        int scale_{};

    };
}