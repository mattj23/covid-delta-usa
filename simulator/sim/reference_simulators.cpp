#include "reference_simulators.hpp"

sim::MethodBase::MethodBase(long population) {
    for (long i = 0; i < population; ++i) {
        pop_.emplace_back();
    }
}

long sim::MethodBase::TotalInfected() const {
    return std::count_if(pop_.begin(), pop_.end(), [](Person p) { return p.IsInfected(); });
}

void sim::NaiveMethod::Run(int initial_infected, int days, double contact_prob) {
    int today = 0;

    // Seed the initial infected
    for (int i = 0; i < initial_infected; ++i) {
        pop_[i].variant = data::Variant::Alpha;
        pop_[i].infected_day = today;
        pop_[i].symptom_onset = prob_.GetAlphaIncubation() + today;
    }

    auto normalized_contact_prob = contact_prob / (double)pop_.size();

    while (today < days) {
        for (const auto &p : pop_) {
            // If this person isn't infected, we don't need to continue
            if (!p.IsInfected()) continue;

            // Compute how infective this person is
            auto infection_p = sim::Probabilities::GetAlphaInfectivity(today - p.symptom_onset);

            // Compute against the rest of the population
            for (auto &k : pop_) {
                // No infected-to-infected transmission considered
                if (k.IsInfected()) continue;
                if (!prob_.UniformChance(normalized_contact_prob)) continue;

                // Roll the dice against the infection probability
                if (prob_.UniformChance(infection_p)) {
                    // Person caught covid
                    k.variant = data::Variant::Alpha;
                    k.infected_day = today;
                    k.symptom_onset = prob_.GetAlphaIncubation() + today;
                }
            }
        }

        printf("%i, %li\n", today, TotalInfected());
        today++;
    }


}

sim::NaiveMethod::NaiveMethod(long population) : MethodBase(population) { }

sim::OptimizedMethod::OptimizedMethod(long population) : MethodBase(population) { }

void sim::OptimizedMethod::Run(int initial_infected, int days, double contact_prob) {
    int today = 0;

    std::unordered_set<size_t> infective;

    // Seed the initial infected
    for (int i = 0; i < initial_infected; ++i) {
        pop_[i].variant = data::Variant::Alpha;
        pop_[i].infected_day = today;
        pop_[i].symptom_onset = prob_.GetAlphaIncubation() + today;
        infective.insert(i);
    }

    auto normalized_contact_prob = contact_prob / (double)pop_.size();
    std::binomial_distribution<int> contact_dist{(int)pop_.size(), normalized_contact_prob};
    std::uniform_int_distribution<int> selector_dist{0, (int)pop_.size()};

    std::vector<size_t> no_longer_infectious;

    while (today < days) {
        for (size_t index : infective) {
            // This is the potential infected carrier
            const auto& p = pop_[index];

            // How infectious are they today
            auto infection_p = sim::Probabilities::GetAlphaInfectivity(today - p.symptom_onset);

            // Check if this guy has passed the point of being infectious
            if (infection_p <= 0 && today > p.symptom_onset) {
                no_longer_infectious.push_back(index);
                continue;
            }

            // Randomly determine how many contacts this person had during the past day, we can
            // move onto the next person if we don't have any
            auto contact_count = contact_dist(prob_.GetGenerator());
            if (!contact_count) continue;

            // Now we'll iterate through that number of contacts, picking someone from the population at random
            // to act as the person who had contact with this carrier.
            for (int i = 0; i < contact_count; ++i) {
                auto selected = selector_dist(prob_.GetGenerator());

                // If the contacted person is already infected, we can move on
                if (pop_[selected].IsInfected()) continue;

                // At this point we know the contacted person is vulnerable to infection, so we roll the dice based
                // on how infectious the carrier is today
                if (prob_.UniformChance(infection_p)) {
                    pop_[selected].variant = data::Variant::Alpha;
                    pop_[selected].infected_day = today;
                    pop_[selected].symptom_onset = prob_.GetAlphaIncubation() + today;
                    infective.insert(selected);
                }
            }
        }

        // Now we can remove the people we marked to be taken off the infectious list
        for (auto index : no_longer_infectious) {
            infective.erase(index);
        }
        no_longer_infectious.clear();

        printf("%i, %li\n", today, TotalInfected());
        today++;
    }

}
