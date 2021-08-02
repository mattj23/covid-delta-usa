#include "probabilities.hpp"

int sim::Probabilities::GetAlphaIncubation() {
    auto value = chance_(generator_);
    for (int i = 0; i < consts::kAlphaIncubation.size(); ++i) {
        if (consts::kAlphaIncubation[i] > value) return i;
    }
    return consts::kAlphaIncubation.size();
}

double sim::Probabilities::GetAlphaInfectivity(int days_from_symptoms) {
    if (days_from_symptoms < consts::kAlphaInfectivityStart) return 0;
    if (days_from_symptoms > consts::kAlphaInfectivityEnd) return 0;
    return consts::kAlphaInfectivity[days_from_symptoms + (-consts::kAlphaInfectivityStart)];
}

int sim::Probabilities::GetTestingLag(double factor) {
    return (int)((double)GetAlphaIncubation() * factor);
}

double sim::Probabilities::GetAlphaVaxEfficacy(int days_from_vax) {
    return consts::kAlphaVaxEfficacy[std::max(std::min(27, days_from_vax), 0)];
}

int sim::Probabilities::GetIncubation(sim::data::Variant variant) {
    if (variant == data::Variant::Alpha) return GetAlphaIncubation();
    if (variant == data::Variant::Delta) return (int)std::round(delta_incubation_ratio * GetAlphaIncubation());
    return 0;
}

double sim::Probabilities::GetInfectivity(sim::data::Variant variant, int days_from_symptoms) {
    if (variant == data::Variant::Alpha) return GetAlphaInfectivity(days_from_symptoms);
    if (variant == data::Variant::Delta) return delta_infectivity_ratio * GetAlphaInfectivity(days_from_symptoms);
    return 0;
}
