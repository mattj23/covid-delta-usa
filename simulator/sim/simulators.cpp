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
    step.virus_carriers = population.CurrentlyInfectious();

    // Expensive summary statistics
    if (expensive) {
        step.population_infectiousness = 0;
        for (int i = 0; i < population.EndOfInfectious(); i++) {
            const auto &person = population.people[i];
            step.population_infectiousness +=
                variants_->at(person.variant)->GetInfectivity(population.today - person.symptom_onset);
        }
        step.population_infectiousness *= population.Scale();
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
    population.AddToInfected(person_index);

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
    auto search_position = population.EndOfInfectious();

        // Scan forward
    while (to_be_vaxxed > population.total_vaccinated) {
        auto& person = population.people[search_position];
        if (!person.is_vaccinated) {
            if (!person.IsInfected() || (population.today - person.infected_day > 30)) {
                person.is_vaccinated = true;
                person.vaccination_day = population.today;
                person.vaccine_immunity_scalar = (float)prob_.UniformScalar();
                population.total_vaccinated++;
            }
        }

        search_position++;
        if (search_position >= population.people.size()) {
            break;
        }
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
    std::vector<size_t> to_remove;
    for (int i = (int)infected_pointer - 1; i >= 0; --i) {
        const auto &person = population.people[i];
        int days_from_symptoms = population.today - person.symptom_onset;
        if (days_from_symptoms > 0 && variants_->at(person.variant)->GetInfectivity(days_from_symptoms) <= 0) {
            to_remove.push_back(i);
        }
    }

    for (auto index : to_remove) {
        population.RemoveFromInfected(index);
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
#ifdef PERF_MEASURE
    loop_timer.Start();
#endif

#pragma omp parallel
{
#ifdef PERF_MEASURE
    PerfTimer t_alloc;
    t_alloc.Start();
#endif
    Probabilities prob;
    std::vector<size_t> local_no_longer_infectious;
    std::vector<std::tuple<size_t, Variant>> local_to_infect;
#ifdef PERF_MEASURE
    t_alloc.Stop();
#endif

#pragma omp for
    for (int carrier_index = 0; carrier_index < population.EndOfInfectious(); carrier_index++) {
        const auto &carrier = population.people[carrier_index];

        // How infectious are they today
        const auto &variant_info = variants_->at(carrier.variant);
        auto infection_p = variant_info->GetInfectivity(population.today - carrier.symptom_onset);

        // Check if this guy has passed the point of being infectious
        if (infection_p <= 0 && population.today > carrier.symptom_onset) {
            local_no_longer_infectious.push_back(carrier_index);
            continue;
        }

        // Randomly determine how many contacts this person had during the past day, we can
        // move onto the next person if we don't have any
        auto contact_count = self_contact_dist(prob.GetGenerator());
        if (!contact_count)
            continue;

        // Now we'll iterate through that number of contacts, picking someone from the population at random
        // to act as the person who had contact with this carrier.
        for (int i = 0; i < contact_count; ++i) {
            // Randomly pick a member of the population
            auto contact_index = selector_dist(prob.GetGenerator());
            const auto &contact = population.people[contact_index];
            if (contact_index < population.EndOfInfectious()) continue;

            // If the carrier's roll for infection doesn't succeed, continue
            if (!prob.UniformChance(infection_p))
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

            local_to_infect.emplace_back(contact_index, carrier.variant);
        }
    }

    #pragma omp critical (sim_day_merge)
    {
        no_longer_infectious.insert(no_longer_infectious.end(), local_no_longer_infectious.begin(), local_no_longer_infectious.end());
        to_infect.insert(to_infect.end(), local_to_infect.begin(), local_to_infect.end());
#ifdef PERF_MEASURE
        alloc += t_alloc.Elapsed();
#endif
    }

}

#ifdef PERF_MEASURE
    loop_timer.Stop();
    remove_timer.Start();
#endif

    // Remove people from the cache who are no long infectious. This has to be done from largest to smallest, in order
    // to prevent the mechanism from moving a person at the end of the list to somewhere else
    std::sort(no_longer_infectious.begin(), no_longer_infectious.end(), std::greater<>());
    for (auto index : no_longer_infectious) {
        population.RemoveFromInfected(index);
    }

#ifdef PERF_MEASURE
    remove_timer.Stop();
    infect_timer.Start();
#endif

    // Add the newly infected. This has to be done from smallest to largest, to prevent the infectious_ptr_ from
    // advancing beyond the people to be infected at the front of the list, sending them off to elsewhere
    std::sort(to_infect.begin(), to_infect.end());
    size_t last_infected = population.people.size() + 1;
    for (const auto &[selected, variant] : to_infect) {
        // This mechanism prevents the same person from being infected multiple times, which won't work because someone
        // else is in that index after the swap
        if (selected == last_infected) continue;

        InfectPerson(population, selected, *variants_->at(variant));
        last_infected = selected;
    }
#ifdef PERF_MEASURE
    infect_timer.Stop();
#endif

    auto result = GetDailySummary(population, options_.expensive_stats);

    population.today++;
    return result;
}
