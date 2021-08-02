#include "simulators.hpp"

sim::StateSimulator::StateSimulator(long population, int scale) {
    scale_ = scale;
    long scaled_population = population / scale_;
    for (long i = 0; i < scaled_population; ++i) {
        pop_.emplace_back();
    }

    // At this point we can set the selector distribution
    selector_dist_ = std::make_unique<std::uniform_int_distribution<int>>(0, pop_.size());
}

void sim::StateSimulator::ResetPopulationTo(date::sys_days reset_date) {
    today_ = data::ToReferenceDate(reset_date);
    vaccine_saves_ = 0;

    for (auto &p : pop_) {
        if (p.infected_day >= today_) {
            p.symptom_onset = 0;
            p.variant = data::Variant::None;
            p.infected_day = 0;
        }

        if (p.vaccinated.has_value() && p.vaccinated.value() >= today_) {
            p.vaccinated = {};
        }
    }
    SynchronizeInfectedCache();
}

void sim::StateSimulator::InitializePopulation(const std::unordered_map<int, data::InfectedHistory> &history,
                                               const std::vector<data::VariantRecord> &variants,
                                               std::optional<date::sys_days> up_to) {
    // Start by resetting the population completely
    ResetPopulationTo(data::kReferenceZeroDate);
    SynchronizeUnVaxxedCache();

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
        // The total number of infected people is the size of the infectious_ set, and for the initialization step can
        // also be assumed to be the value of the infected_pointer. At each step in the history the number of infected
        // should be the number of people who have TESTED positive.  That means that when we "infect" someone while
        // initializing the history, we set the date tested to today (the data should be corrected to count a case as
        // the date that the test was done) and then back-calculate a guess at when the infection actually occurred.
        auto h = history.find(working_day);
        if (h == history.end()) continue;

        today_ = working_day;
        auto variant_fractions = data::GetVariantFractions(working_day, variants);

        // The number of infections we need
        long scaled_infections = (long)(h->second.total_infections / scale_);
        long total_to_add = scaled_infections - (long)infected_pointer;

        for (const auto &[variant, fraction] : variant_fractions) {
            int to_add = (int)std::round(fraction * (double)total_to_add);

            for (int k = 0; k < to_add; ++k) {
                auto &p = pop_[infected_pointer];
                p.variant = variant;
                p.infected_day = today_;
                p.symptom_onset = p.infected_day + prob_.GetIncubation(variant);
                p.test_day = today_ + prob_.GetTestingLag();

                infected_pointer++;
            }

        }

        working_day++;
    }

    SynchronizeInfectedCache();
}

void sim::StateSimulator::SynchronizeInfectedCache() {
    infectious_.clear();
    for (size_t i = 0; i < pop_.size(); ++i) {
        if (pop_[i].IsInfected()) infectious_.insert(i);
    }
}

void sim::StateSimulator::SynchronizeUnVaxxedCache() {
    un_vaxxed_.clear();
    for (size_t i = 0; i < pop_.size(); ++i) {
        if (!pop_[i].vaccinated.has_value()) un_vaxxed_.push_back(i);
    }
}

void sim::StateSimulator::SetProbabilities(double p_self) {
    auto normalized_contact_prob = p_self / (double)pop_.size();
    self_contact_dist_ = std::make_unique<std::binomial_distribution<int>>((int)pop_.size(), normalized_contact_prob);
}

void sim::StateSimulator::SimulateDay() {
    std::vector<size_t> no_longer_infectious;
    std::vector<std::tuple<size_t, data::Variant>> to_infect;

    // First, calculate the new infections, which will be applied in a later step
    for (size_t index : infectious_) {
        const auto& p = pop_[index];

        // How infectious are they today
        auto infection_p = sim::Probabilities::GetInfectivity(p.variant, today_ - p.symptom_onset);

        // Check if this guy has passed the point of being infectious
        if (infection_p <= 0 && today_ > p.symptom_onset) {
            no_longer_infectious.push_back(index);
            continue;
        }

        // Randomly determine how many contacts this person had during the past day, we can
        // move onto the next person if we don't have any
        auto contact_count = (*self_contact_dist_)(prob_.GetGenerator());
        if (!contact_count) continue;

        // Now we'll iterate through that number of contacts, picking someone from the population at random
        // to act as the person who had contact with this carrier.
        for (int i = 0; i < contact_count; ++i) {
            // Randomly pick a member of the population
            auto selected = (*selector_dist_)(prob_.GetGenerator());

            // If the contacted person is already infected, we can move on
            if (pop_[selected].IsInfected()) continue;

            // At this point we know the contacted person is vulnerable to infection, so we roll the dice based
            // on how infectious the carrier is today
            if (prob_.UniformChance(infection_p)) {
                to_infect.emplace_back(selected, p.variant);
            }
        }
    }

    // Remove people from the cache who are no long infectious
    for (auto index : no_longer_infectious) {
        infectious_.erase(index);
    }

    // Add the newly infected
    for (const auto &[selected, variant] : to_infect) {
        auto &person = pop_[selected];

        // Vaccine saving throw?
        if (person.vaccinated.has_value()) {
            auto efficacy = Probabilities::GetAlphaVaxEfficacy(today_ - person.vaccinated.value());
            if (prob_.UniformChance(efficacy)) {
                vaccine_saves_++;
                continue;
            }
        }

        person.variant = variant;
        person.infected_day = today_;
        person.symptom_onset = prob_.GetAlphaIncubation() + today_;
        person.test_day = today_ + prob_.GetTestingLag();
        infectious_.insert(selected);
    }

    today_++;
}

size_t sim::StateSimulator::TotalInfectious() {
    return infectious_.size();
}

size_t sim::StateSimulator::TestedPositive(int on_day) {
    return std::count_if(pop_.begin(), pop_.end(),
                         [&on_day](const Person& p){ return p.IsInfected() && p.test_day <= on_day; }) * scale_;
}

size_t sim::StateSimulator::TotalInfected(int on_day) {
    return std::count_if(pop_.begin(), pop_.end(),
                         [&on_day](const Person& p){ return p.IsInfected() && p.infected_day <= on_day; }) * scale_;
}

size_t sim::StateSimulator::TotalVaccinated(int on_day) {
    return std::count_if(pop_.begin(), pop_.end(),
                         [&on_day](const Person& p){
        return p.vaccinated.has_value() && p.vaccinated.value() <= on_day;
    }) * scale_;
}

size_t sim::StateSimulator::TotalWithDelta(int on_day) {
    return std::count_if(pop_.begin(), pop_.end(),
                         [&on_day](const Person& p){
                             return p.variant == data::Variant::Delta && p.infected_day <= on_day;
                         }) * scale_;
}

void sim::StateSimulator::ApplyVaccines(const std::unordered_map<int, data::VaccineHistory> &vaccines, int for_date) {
    // Look ahead two weeks and add the vaccinations now
    auto shifted = for_date + 14;
    auto vax = vaccines.find(shifted);
    if (vax == vaccines.end()) return;

    const auto& today_data = vax->second;
    int to_be_vaxxed = today_data.total_completed_vax / scale_;
    int vaccinated = 0;

    while ((pop_.size() - un_vaxxed_.size()) < to_be_vaxxed) {
        auto s = un_vaxxed_.size();
        std::uniform_int_distribution<int> dist(0, (int)s);
        int index = dist(prob_.GetGenerator());
        auto person_index = un_vaxxed_[index];
        if (index < un_vaxxed_.size() - 1)
            std::swap(un_vaxxed_[index], un_vaxxed_.back());
        un_vaxxed_.pop_back();

        // Pick someone at random to get a vax_history
        auto &person = pop_[person_index];
        if (person.vaccinated.has_value()) continue;
        if (person.IsInfected() && (for_date - person.infected_day) < 30) continue;

        person.vaccinated = for_date;
    }
}

void sim::StateSimulator::InitializeVaccines(const std::unordered_map<int, data::VaccineHistory> &vaccines,
                                             std::optional<date::sys_days> up_to) {
    // Get the min and max days of interest
    int working_day = std::numeric_limits<int>::max();
    int max_day = std::numeric_limits<int>::min();
    for (const auto &[day, _] : vaccines) {
        working_day = std::min(working_day, day);
        max_day = std::max(max_day, day);
    }

    if (up_to.has_value()) {
        max_day = data::ToReferenceDate(up_to.value());
    }

    while (working_day < max_day) {
        ApplyVaccines(vaccines, working_day);
        working_day++;
    }

}

void sim::StateSimulator::ApplyTodaysVaccines(const std::unordered_map<int, data::VaccineHistory> &vaccines) {
    ApplyVaccines(vaccines, today_);
}






