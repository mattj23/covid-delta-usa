#include "variant_probabilities.hpp"

sim::VariantProbabilities::VariantProbabilities(const data::VariantProperties& variant_properties, Variant variant)
        : variant_(variant),
        incubation_(variant_properties.incubation.begin(), variant_properties.incubation.end()),
        properties_(variant_properties) {}

double sim::VariantProbabilities::GetInfectivity(int days_from_symptoms) const {
    return properties_.infectivity(days_from_symptoms);
}

int sim::VariantProbabilities::GetRandomIncubation(std::mt19937_64 &mt) const {
    std::uniform_real_distribution<double> dist(0, 1);
    auto value = dist(mt);
    for (int i = 0; i < incubation_.size(); ++i) {
        if (value <= incubation_[i])
            return i;
    }
    return (int)incubation_.size();
}

double sim::VariantProbabilities::GetVaxImmunity(int days_from_vax) const {
    return properties_.vax_immunity(days_from_vax);
}

double sim::VariantProbabilities::GetNaturalImmunity(int days_from_infection) const {
    return properties_.natural_immunity(days_from_infection);
}

bool sim::VariantProbabilities::IsPersonNatImmune(const sim::Person &person, int today) const {
    return person.IsInfected() && person.natural_immunity_scalar <= GetNaturalImmunity(today - person.infected_day);
}

bool sim::VariantProbabilities::IsPersonVaxImmune(const sim::Person &person, int today) const {
    return person.is_vaccinated &&
           person.vaccine_immunity_scalar <= GetVaxImmunity(today - person.vaccination_day);
}
