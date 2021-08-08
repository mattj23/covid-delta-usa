#include "simulators.hpp"

sim::StateSimulator::StateSimulator(long population, int scale, const VariantDictionary *variants) {
    scale_ = scale;
    long scaled_population = population / scale_;
    for (long i = 0; i < scaled_population; ++i) {
        pop_.emplace_back();
    }
    variants_ = variants;

    // At this point we can set the selector distribution
    selector_dist_ = std::make_unique<std::uniform_int_distribution<int>>(0, pop_.size());
}

void sim::StateSimulator::Reset() {
    vaccine_saves_ = 0;
    total_infections_ = 0;
    total_vaccinated_ = 0;
    never_infected_ = (int)pop_.size();
    total_delta_infections_ = 0;
    total_alpha_infections_ = 0;

    today_ = 0;

    for (auto &person : pop_) {
        person.Reset();
    }

    infectious_.clear();
    SynchronizeUnVaxxedCache();
}

std::vector<sim::data::StepResult>
sim::StateSimulator::InitializePopulation(const std::unordered_map<int, data::InfectedHistory> &history,
                                          const std::unordered_map<int, data::VaccineHistory> &vaccines,
                                          const std::vector<data::VariantRecord> &variant_history,
                                          std::optional<date::sys_days> up_to) {
    // Start by resetting the population completely
    Reset();

    // Prepare the results
    std::vector<data::StepResult> records;

    size_t infected_pointer = 0;

    // Get the min and max days of interest
    int working_day = std::numeric_limits<int>::max();
    int max_day = std::numeric_limits<int>::min();
    for (const auto &[day, _] : history) {
        working_day = std::min(working_day, day);
        max_day = std::max(max_day, day);
    }

    if (up_to.has_value()) {
        max_day = data::ToReferenceDate(up_to.value());
    }

    while (working_day < max_day) {
        auto h = history.find(working_day);
        if (h == history.end())
            continue;

        // Todo: Remove the working_day value and replace it with today_ everywhere it's used, and check that nothing
        // breaks when that happens
        today_ = working_day;
        auto variant_fractions = data::GetVariantFractions(working_day, variant_history);

        // The number of infections we need
        long scaled_infections = (long)(h->second.total_infections / scale_);
        long total_to_add = scaled_infections - (long)infected_pointer;

        for (const auto &[variant, fraction] : variant_fractions) {
            int to_add = (int)std::round(fraction * (double)total_to_add);

            for (int k = 0; k < to_add; ++k) {
                InfectPerson(infected_pointer, *variants_->at(variant));
                infected_pointer++;
            }
        }

        // Remove anyone who's no longer infectious
        for (int i = 0; i < infected_pointer; ++i) {
            const auto &person = pop_[i];
            int days_from_symptoms = today_ - person.symptom_onset;
            if (days_from_symptoms > 0 && variants_->at(person.variant)->GetInfectivity(days_from_symptoms) <= 0) {
                infectious_.erase(i);
            }
        }

        ApplyVaccines(vaccines);

        // If the options are set to export the full history, we do it here
        if (options_ != nullptr && options_->full_history) {
            records.push_back(GetStepResult());
//            printf("init: %i %u %u -> %u\n", records.back().year, records.back().month, records.back().day, records.back().total_infections);
        }

        working_day++;
    }

    return records;
}

void sim::StateSimulator::InfectPerson(size_t person_index, const VariantProbabilities &variant) {
    auto &person = pop_[person_index];
    if (person.variant == data::Variant::None) {
        never_infected_--;
    } else {
        reinfections_++;
    }

    if (person.vaccinated.has_value()) {
        vaccinated_infections_++;
    }

    person.variant = variant.GetVariant();
    person.infected_day = today_;
    person.symptom_onset = today_ + variant.GetRandomIncubation(prob_.GetGenerator());
    person.natural_immunity_scalar = (float)prob_.UniformScalar();
    infectious_.insert(person_index);
    total_infections_++;

    if (person.variant == data::Variant::Delta)
        total_delta_infections_++;
    if (person.variant == data::Variant::Alpha)
        total_alpha_infections_++;
}

void sim::StateSimulator::SynchronizeUnVaxxedCache() {
    un_vaxxed_.clear();
    for (size_t i = 0; i < pop_.size(); ++i) {
        if (!pop_[i].vaccinated.has_value())
            un_vaxxed_.push_back(i);
    }
}

void sim::StateSimulator::SetProbabilities(double p_self) {
    auto normalized_contact_prob = p_self / (double)pop_.size();
    self_contact_dist_ = std::make_unique<std::binomial_distribution<int>>((int)pop_.size(), normalized_contact_prob);
}

sim::data::StepResult sim::StateSimulator::SimulateDay() {
    std::vector<size_t> no_longer_infectious;
    std::vector<std::tuple<size_t, data::Variant>> to_infect;

    // First, calculate the new infections, which will be applied in a later step
    for (size_t carrier_index : infectious_) {
        const auto &carrier = pop_[carrier_index];

        // How infectious are they today
        const auto &variant_info = variants_->at(carrier.variant);
        auto infection_p = variant_info->GetInfectivity(today_ - carrier.symptom_onset);

        // Check if this guy has passed the point of being infectious
        if (infection_p <= 0 && today_ > carrier.symptom_onset) {
            no_longer_infectious.push_back(carrier_index);
            continue;
        }

        // Randomly determine how many contacts this person had during the past day, we can
        // move onto the next person if we don't have any
        auto contact_count = (*self_contact_dist_)(prob_.GetGenerator());
        if (!contact_count)
            continue;

        // Now we'll iterate through that number of contacts, picking someone from the population at random
        // to act as the person who had contact with this carrier.
        for (int i = 0; i < contact_count; ++i) {
            // Randomly pick a member of the population
            auto contact_index = (*selector_dist_)(prob_.GetGenerator());
            const auto &contact = pop_[contact_index];

            // If the carrier's roll for infection doesn't succeed, continue
            if (!prob_.UniformChance(infection_p))
                continue;

            // At this point the carrier has successfully rolled to infect the contact. Now we will see if the contact
            // has an immunity which can prevent the infection.
            // Check if they have natural immunity
            if (variant_info->IsPersonNatImmune(contact, today_)) {
                natural_saves_++;
                continue;
            }

            // Check if they have vaccine immunity
            if (variant_info->IsPersonVaxImmune(contact, today_)) {
                vaccine_saves_++;
                continue;
            }

            // At this point we know the contacted person is vulnerable to infection, so we roll the dice based
            // on how infectious the carrier is today
            to_infect.emplace_back(contact_index, carrier.variant);
        }
    }

    // Remove people from the cache who are no long infectious
    for (auto index : no_longer_infectious) {
        infectious_.erase(index);
    }

    // Add the newly infected
    for (const auto &[selected, variant] : to_infect) {
        InfectPerson(selected, *variants_->at(variant));
    }

    auto result = GetStepResult();

    today_++;
    return result;
}

void sim::StateSimulator::ApplyVaccines(const std::unordered_map<int, data::VaccineHistory> &vaccines) {
    // The issue that we have with the vaccine history is that it's taking into account "completed" vaccinations, which
    // means the patient has received the second shot.  Because we can't track vaccinations individually, the
    // approximation chosen here is to look at completed vaccinations 21 days in the future and apply those vaccinations
    // today. The immunity will begin to ramp up according to the efficacy curves.

    auto shifted = today_ + 21;
    auto vax = vaccines.find(shifted);
    if (vax == vaccines.end())
        return;

    const auto &today_data = vax->second;
    int to_be_vaxxed = today_data.total_completed_vax / scale_;

    while ((pop_.size() - un_vaxxed_.size()) < to_be_vaxxed) {
        // Pick someone at random to get a vax_history
        std::uniform_int_distribution<int> dist(0, (int)un_vaxxed_.size());
        int index = dist(prob_.GetGenerator());
        auto person_index = un_vaxxed_[index];

        if (index < un_vaxxed_.size() - 1)
            std::swap(un_vaxxed_[index], un_vaxxed_.back());
        un_vaxxed_.pop_back();

        auto &person = pop_[person_index];
        if (person.vaccinated.has_value())
            continue;
        if (person.IsInfected() && (today_ - person.infected_day) < 30)
            continue;

        person.vaccinated = today_;
        person.vaccine_immunity_scalar = (float)prob_.UniformScalar();
        total_vaccinated_++;
    }
}

sim::data::StepResult sim::StateSimulator::GetStepResult() const {
    data::StepResult step{};
    date::year_month_day d = data::ToSysDays(today_);
    step.year = (int)d.year();
    step.month = d.month().operator unsigned int();
    step.day = d.day().operator unsigned int();

    step.total_infections = total_infections_ * scale_;
    step.total_vaccinated = total_vaccinated_ * scale_;
    step.vaccine_saves = vaccine_saves_ * scale_;
    step.natural_saves = natural_saves_ * scale_;
    step.never_infected = never_infected_ * scale_;
    step.total_delta_infections = total_delta_infections_ * scale_;
    step.total_alpha_infections = total_alpha_infections_ * scale_;
    step.reinfections = reinfections_ * scale_;
    step.vaccinated_infections = vaccinated_infections_ * scale_;
    step.virus_carriers = static_cast<int>(infectious_.size());

    // Expensive summary statistics
    if (false) {
        step.population_infectiousness = 0;
        for (auto index : infectious_) {
            const auto &person = pop_[index];
            step.population_infectiousness +=
                variants_->at(person.variant)->GetInfectivity(today_ - person.symptom_onset);
        }
    }

    return step;
}
