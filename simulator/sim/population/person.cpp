#include "person.hpp"

void sim::Person::Reset() {
    variant = Variant::None;
    infected_day = 0;
    symptom_onset = 0;
    test_day = 0;
    natural_immunity_scalar = 0;
    vaccine_immunity_scalar = 0;
    is_vaccinated = false;
    vaccination_day = 0;
}
