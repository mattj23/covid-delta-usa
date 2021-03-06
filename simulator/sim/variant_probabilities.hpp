#pragma once
#include <random>
#include "data.hpp"
#include "population/person.hpp"

namespace sim {
    class VariantProbabilities {
    public:
        VariantProbabilities(const data::VariantProperties& variant_properties, Variant variant);

        [[nodiscard]] double GetInfectivity(int days_from_symptoms) const;
        [[nodiscard]] double GetVaxImmunity(int days_from_vax) const;
        [[nodiscard]] double GetNaturalImmunity(int days_from_infection) const;

        [[nodiscard]] int GetRandomIncubation(std::mt19937_64 &mt) const;

        [[nodiscard]] bool IsPersonVaxImmune(const Person& person, int today) const;
        [[nodiscard]] bool IsPersonNatImmune(const Person& person, int today) const;

        [[nodiscard]] Variant GetVariant() const { return variant_; }
    private:
        Variant variant_;
        std::vector<double> incubation_;
        data::VariantProperties properties_;
    };

    using VariantDictionary = std::unordered_map<Variant, std::unique_ptr<VariantProbabilities>>;
}