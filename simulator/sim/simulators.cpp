#include "simulators.hpp"

sim::Simulator::Simulator(const data::ProgramOptions &options, std::shared_ptr<const VariantDictionary> variants)
    : options_(options), variants_(variants) {}

sim::DailySummary sim::Simulator::GetDailySummary(const sim::Population &population, bool expensive) const {
    sim::DailySummary step{};
    step.day = population.today;

    step.total_infections = population.total_infections * population.Scale();
    step.total_vaccinated = population.total_vaccinated * population.Scale();
    step.vaccine_saves = population.vaccine_saves * population.Scale();
    step.natural_saves = population.natural_saves * population.Scale();
    step.never_infected = population.never_infected * population.Scale();
    step.total_delta_infections = population.total_delta_infections * population.Scale();
    step.total_alpha_infections = population.total_alpha_infections * population.Scale();
    step.reinfections = population.reinfections * population.Scale();
    step.vaccinated_infections = population.vaccinated_infections * population.Scale();
    step.virus_carriers = static_cast<int>(population.infectious_indices.size()) * population.Scale();

    // Expensive summary statistics
    if (expensive) {
        step.population_infectiousness = 0;
        for (auto index : population.infectious_indices) {
            const auto &person = population.people[index];
            step.population_infectiousness +=
                variants_->at(person.variant)->GetInfectivity(population.today - person.symptom_onset);
        }
    }

    return step;
}

void sim::Simulator::InfectPerson(sim::Population &population, size_t person_index,
                                  const VariantProbabilities &variant) {

    auto &person = population.people[person_index];
    if (person.variant == Variant::None) {
        population.never_infected--;
    } else {
        population.reinfections++;
    }

    if (person.is_vaccinated) {
        population.vaccinated_infections++;
    }

    person.variant = variant.GetVariant();
    person.infected_day = population.today;
    person.symptom_onset = population.today + variant.GetRandomIncubation(prob_.GetGenerator());
    person.natural_immunity_scalar = (float)prob_.UniformScalar();
    population.infectious_indices.insert(person_index);

    // Swap the person with the person at the end of the infectious indices and advance the pointer
    if (population.infectious_ptr_ != person_index) {
    }


    population.infectious_ptr_++;

    population.total_infections++;

    if (person.variant == Variant::Delta)
        population.total_delta_infections++;
    if (person.variant == Variant::Alpha)
        population.total_alpha_infections++;
}

void sim::Simulator::ApplyVaccines(sim::Population &population,
                                   const std::unordered_map<int, data::VaccineHistory> &vaccines) {
    // The issue that we have with the vaccine history is that it's taking into account "completed" vaccinations, which
    // means the patient has received the second shot.  Because we can't track vaccinations individually, the
    // approximation chosen here is to look at completed vaccinations 21 days in the future and apply those vaccinations
    // today. The immunity will begin to ramp up according to the efficacy curves.

    auto shifted = population.today + 21;
    auto vax = vaccines.find(shifted);
    if (vax == vaccines.end())
        return;

    const auto &today_data = vax->second;
    int to_be_vaxxed = today_data.total_completed_vax / population.Scale();

    while ((population.people.size() - population.unvaxxed_indices.size()) < to_be_vaxxed) {
        // Pick someone at random to get a vax_history
        std::uniform_int_distribution<int> dist(0, static_cast<int>(population.unvaxxed_indices.size()));
        int index = dist(prob_.GetGenerator());
        auto person_index = population.unvaxxed_indices[index];

        if (index < population.unvaxxed_indices.size() - 1)
            std::swap(population.unvaxxed_indices[index], population.unvaxxed_indices.back());
        population.unvaxxed_indices.pop_back();

        auto &person = population.people[person_index];
        if (person.is_vaccinated)
            continue;
        if (person.IsInfected() && (population.today - person.infected_day) < 30)
            continue;

        person.is_vaccinated = true;
        person.vaccination_day = population.today;
        person.vaccine_immunity_scalar = (float)prob_.UniformScalar();
        population.total_vaccinated++;
    }
}

std::vector<sim::DailySummary> sim::Simulator::InitializePopulation(
    sim::Population &population, const std::unordered_map<int, data::InfectedHistory> &history,
    const std::unordered_map<int, data::VaccineHistory> &vaccines,
    const std::vector<data::VariantRecord> &variant_history, std::optional<date::sys_days> up_to) {

    population.Reset();
    std::vector<DailySummary> summaries;

    size_t infected_pointer = 0;

    // Get the min and max days of interest
    population.today = std::numeric_limits<int>::max();
    int max_day = std::numeric_limits<int>::min();
    for (const auto &[day, _] : history) {
        population.today = std::min(population.today, day);
        max_day = std::max(max_day, day);
    }

    if (up_to.has_value()) {
        max_day = data::ToReferenceDate(up_to.value());
    }

    while (population.today < max_day) {
        auto h = history.find(population.today);
        if (h == history.end())
            continue;

        auto variant_fractions = data::GetVariantFractions(population.today, variant_history);

        // The number of infections we need
        long scaled_infections = (long)(h->second.total_infections / population.Scale());
        long total_to_add = scaled_infections - (long)infected_pointer;

        for (const auto &[variant, fraction] : variant_fractions) {
            int to_add = (int)std::round(fraction * (double)total_to_add);

            for (int k = 0; k < to_add; ++k) {
                InfectPerson(population, infected_pointer, *variants_->at(variant));
                infected_pointer++;
            }
        }

        ApplyVaccines(population, vaccines);

        // If the options are set to export the full history, we do it here
        if (options_.full_history) {
            summaries.push_back(GetDailySummary(population, options_.expensive_stats));
        }

        population.today++;
    }

    // Remove anyone who's no longer infectious
    for (int i = 0; i < infected_pointer; ++i) {
        const auto &person = population.people[i];
        int days_from_symptoms = population.today - person.symptom_onset;
        if (days_from_symptoms > 0 && variants_->at(person.variant)->GetInfectivity(days_from_symptoms) <= 0) {
            population.infectious_indices.erase(i);
        }
    }

    return summaries;
}

sim::DailySummary sim::Simulator::SimulateDay(sim::Population &population) {
    std::vector<size_t> no_longer_infectious;
    std::vector<std::tuple<size_t, Variant>> to_infect;

    auto normalized_contact = contact_probability_ / static_cast<int>(population.people.size());
    std::binomial_distribution<int> self_contact_dist(static_cast<int>(population.people.size()), normalized_contact);
    std::uniform_int_distribution<int> selector_dist(0, static_cast<int>(population.people.size()));

    // First, calculate the new infections, which will be applied in a later step
    for (size_t carrier_index : population.infectious_indices) {
        const auto &carrier = population.people[carrier_index];

        // How infectious are they today
        const auto &variant_info = variants_->at(carrier.variant);
        auto infection_p = variant_info->GetInfectivity(population.today - carrier.symptom_onset);

        // Check if this guy has passed the point of being infectious
        if (infection_p <= 0 && population.today > carrier.symptom_onset) {
            no_longer_infectious.push_back(carrier_index);
            continue;
        }

        // Randomly determine how many contacts this person had during the past day, we can
        // move onto the next person if we don't have any
        auto contact_count = self_contact_dist(prob_.GetGenerator());
        if (!contact_count)
            continue;

        // Now we'll iterate through that number of contacts, picking someone from the population at random
        // to act as the person who had contact with this carrier.
        for (int i = 0; i < contact_count; ++i) {
            // Randomly pick a member of the population
            auto contact_index = selector_dist(prob_.GetGenerator());
            const auto &contact = population.people[contact_index];

            // If the carrier's roll for infection doesn't succeed, continue
            if (!prob_.UniformChance(infection_p))
                continue;

            // At this point the carrier has successfully rolled to infect the contact. Now we will see if the contact
            // has an immunity which can prevent the infection.
            // Check if they have natural immunity
            if (variant_info->IsPersonNatImmune(contact, population.today)) {
                population.natural_saves++;
                continue;
            }

            // Check if they have vaccine immunity
            if (variant_info->IsPersonVaxImmune(contact, population.today)) {
                population.vaccine_saves++;
                continue;
            }

            // At this point we know the contacted person is vulnerable to infection, so we roll the dice based
            // on how infectious the carrier is today
            to_infect.emplace_back(contact_index, carrier.variant);
        }
    }

    // Remove people from the cache who are no long infectious
    for (auto index : no_longer_infectious) {
        population.infectious_indices.erase(index);
    }

    // Add the newly infected
    for (const auto &[selected, variant] : to_infect) {
        InfectPerson(population, selected, *variants_->at(variant));
    }

    auto result = GetDailySummary(population, options_.expensive_stats);

    population.today++;
    return result;
}
