#pragma once
#include <random>
#include "data.hpp"

namespace sim {
    class VariantProbabilities {
    public:
        explicit VariantProbabilities(const data::VariantProperties& variant);

        [[nodiscard]] double GetInfectivity(int days_from_symptoms) const;
        int GetRandomIncubation(std::mt19937_64 &mt);
        double GetVaxImmunity(int days_from_vax) const;
        double GetNaturalImmunity(int days_from_infection) const;

    private:
        std::discrete_distribution<int> incubation_;
        data::VariantProperties properties_;
    };
}