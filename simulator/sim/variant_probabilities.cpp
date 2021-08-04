#include "variant_probabilities.hpp"

sim::VariantProbabilities::VariantProbabilities(const sim::data::VariantProperties &variant)
    : incubation_(variant.incubation.begin(), variant.incubation.end()), properties_(variant) {}

double sim::VariantProbabilities::GetInfectivity(int days_from_symptoms) const {
    return properties_.infectivity(days_from_symptoms);
}

int sim::VariantProbabilities::GetRandomIncubation(std::mt19937_64 &mt) {
    return incubation_(mt);
}

double sim::VariantProbabilities::GetVaxImmunity(int days_from_vax) const {
    return properties_.vax_immunity(days_from_vax);
}

double sim::VariantProbabilities::GetNaturalImmunity(int days_from_infection) const {
    return properties_.natural_immunity(days_from_infection);
}


