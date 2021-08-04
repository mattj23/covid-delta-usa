#pragma once
#include <random>
#include "data.hpp"
#include "person.hpp"

namespace sim {
    class VariantProbabilities {
    public:
        explicit VariantProbabilities(const data::VariantProperties& variant);

        [[nodiscard]] double GetInfectivity(int days_from_symptoms) const;
        [[nodiscard]] double GetVaxImmunity(int days_from_vax) const;
        [[nodiscard]] double GetNaturalImmunity(int days_from_infection) const;

        int GetRandomIncubation(std::mt19937_64 &mt);

        [[nodiscard]] bool IsPersonVaxImmune(const Person& person, int today) const;
        [[nodiscard]] bool IsPersonNatImmune(const Person& person, int today) const;
    private:
        std::discrete_distribution<int> incubation_;
        data::VariantProperties properties_;
    };
}