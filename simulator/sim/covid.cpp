#include "covid.hpp"

void sim::to_json(nlohmann::json &j, const DailySummary &r) {
    j = nlohmann::json{
            {"day", r.day},
            {"vaccine_saves", r.vaccine_saves},
            {"natural_saves", r.natural_saves},
            {"total_infections", r.total_infections},
            {"total_vaccinated", r.total_vaccinated},
            {"never_infected", r.never_infected},
            {"reinfections", r.reinfections},
            {"virus_carriers", r.virus_carriers},
            {"vaccinated_infections", r.vaccinated_infections},
            {"total_delta_infections", r.total_delta_infections},
            {"total_alpha_infections", r.total_alpha_infections},
            {"population_infectiousness", r.population_infectiousness},
    };
}
