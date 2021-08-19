#include "contact_prob.hpp"

sim::ContactProbabilitySearch::ContactProbabilitySearch(const sim::data::ProgramInput &input,
                                                        std::shared_ptr<const sim::VariantDictionary> variants)
    : input_(input), variants_(variants) {}

sim::ContactResult sim::ContactProbabilitySearch::FindContactProbability(int day) {
    auto start_date = data::ToSysDays(day);
    auto state_info = input_.state_info.at(input_.state);

    Simulator simulator(input_.options, variants_);
    Population ref_pop(state_info.population, input_.population_scale, state_info.ages);
    Population work_pop = ref_pop;

    // The starting guess for the contact probability is the value that was supplied

    const auto &infected_history = input_.infected_history.at(input_.state);
    std::vector<int> expected;
    for (int i = 0; i < kCheckDays; ++i) {
        int i0 = infected_history.at(sim::data::ToReferenceDate(start_date + date::days{i - 1})).total_infections;
        int i1 = infected_history.at(sim::data::ToReferenceDate(start_date + date::days{i})).total_infections;
        expected.push_back(i1 - i0);
    }

    // Initialize the population from the beginning
    auto init_start = std::chrono::system_clock::now();
    auto init_result = simulator.InitializePopulation(ref_pop, infected_history, input_.vax_history.at(input_.state),
                                                      input_.variant_history.at(input_.state), start_date);

    // Get the upper and lower bounds
    auto step0 = GetResultFromBounds(ref_pop, work_pop, expected, simulator, start_date, 2.0, 0.5);

    double upper = step0.prob + 3 * step0.stdev;
    double lower = step0.prob - 3 * step0.stdev;
    auto result = GetResultFromBounds(ref_pop, work_pop, expected, simulator, start_date, upper, lower);


    return result;
}

sim::ContactResult
sim::ContactProbabilitySearch::GetResultFromBounds(const sim::Population &reference_pop, sim::Population &working_pop,
                                                   const std::vector<int> &expected, sim::Simulator &simulator,
                                                   date::sys_days start_date, double upper, double lower) {
    double step = (upper - lower) / input_.run_count;
    std::vector<double> xs; // Contact probabilities
    std::vector<double> ys; // Errors

    total_timer.Start();

    for (int run = 0; run < input_.run_count; ++run) {
        std::vector<int> results;

        copy_timer.Start();
        working_pop.CopyFrom(reference_pop);
        copy_timer.Stop();

        // Setting the contact probability
        double contact_prob = lower + (step * run);
        simulator.SetProbabilities(contact_prob);
        int last_infections = working_pop.TotalInfections();

        auto today = start_date;
        while (today < start_date + date::days{kCheckDays}) {
            // Add the newly vaccinated
            if (!input_.vax_history.empty()) {
                vax_timer.Start();
                simulator.ApplyVaccines(working_pop, input_.vax_history.at(input_.state));
                vax_timer.Stop();
            }

            // Simulate the day's new infections and record them
            sim_timer.Start();
            simulator.SimulateDay(working_pop);
            sim_timer.Stop();
            auto new_infections = working_pop.TotalInfections() - last_infections;
            results.push_back(new_infections);
            last_infections = working_pop.TotalInfections();

            // Increment the clock
            today += date::days{1};
        }

        double error{};
        for (int i = 0; i < results.size(); ++i) {
            error += (results[i] - expected[i]);
        }
        error = error / (kCheckDays);
        xs.push_back(contact_prob);
        ys.push_back(error);
    }

    // Compute the line of best fit
    double n = static_cast<double>(xs.size());
    double sum_x = 0;
    double sum_xy = 0;
    double sum_x2 = 0;
    double sum_y = 0;

    for (size_t i = 0; i < xs.size(); ++i) {
        sum_x += xs[i];
        sum_xy += xs[i] * ys[i];
        sum_x2 += xs[i] * xs[i];
        sum_y += ys[i];
    }

    double mean_y = sum_y / n;
    double mean_x = sum_x / n;
    double ss_xx = sum_x2 - (1.0 / n) * sum_x * sum_x;
    double ss_xy = sum_xy - (1.0 / n) * sum_x * sum_y;
    double slope = ss_xy / ss_xx;
    double y0 = mean_y - slope * mean_x;

    // Compute the value at zero
    double x0 = -y0 / slope;

    // Compute the adjusted standard deviation
    double variance = 0;
    for (size_t i = 0; i < xs.size(); ++i) {
        auto value = ys[i] - (slope * xs[i] + y0);
        variance += std::pow(value, 2.);
    }
    variance = variance / n;
    auto stdev = std::sqrt(variance);

    total_timer.Stop();
//    printf("[time] total = %li\n", total_timer.Elapsed());
//    printf("[time] copy = %li\n", copy_timer.Elapsed());
//    printf("[time] vax = %li\n", vax_timer.Elapsed());
//    printf("[time] sim = %li\n", sim_timer.Elapsed());

    return {x0, stdev / slope};
}

void sim::to_json(nlohmann::json &j, const sim::ContactSearchResultSet &o) {
    j = nlohmann::json{{"days", o.days}, {"probabilities", o.probabilities}, {"stdevs", o.stdevs}};
}
