#include "contact_prob.hpp"

sim::ContactProbabilitySearch::ContactProbabilitySearch(const sim::data::ProgramInput& input,
                                                        std::shared_ptr<const sim::VariantDictionary> variants) :
                                                        input_(input), variants_(variants) {



}

sim::ContactResult sim::ContactProbabilitySearch::FindContactProbability(int day) {
    auto start_date = data::ToSysDays(day);
    auto state_info = input_.state_info.at(input_.state);

    Simulator simulator(input_.options, variants_);
    Population reference_population(state_info.population, input_.population_scale);
    Population working_population(state_info.population, input_.population_scale);

    // The starting guess for the contact probability is the value that was supplied
    constexpr int kCheckDays = 3;

    const auto& infected_history = input_.infected_history.at(input_.state);
    std::vector<int> expected;
    for (int i = 0; i < kCheckDays; ++i) {
        int i0 = infected_history.at(sim::data::ToReferenceDate(start_date + date::days{i-1})).total_infections;
        int i1 = infected_history.at(sim::data::ToReferenceDate(start_date + date::days{i})).total_infections;
        expected.push_back(i1 - i0);
    }

    // Initialize the population from the beginning
    auto init_start = std::chrono::system_clock::now();
    auto init_result = simulator.InitializePopulation(reference_population, infected_history,
                                                      input_.vax_history.at(input_.state),
                                                      input_.variant_history.at(input_.state),
                                                      start_date);

    std::vector<double> xs;     // Contact probabilities
    std::vector<double> ys;     // Errors

    // Get the upper and lower bounds
    double upper = 2.0;
    double lower = 0.5;
    double step = (upper - lower) / input_.run_count;

    for (int run = 0; run < input_.run_count; ++run) {
        std::vector<int> results;
        working_population.CopyFrom(reference_population);

        // Setting the contact probability
        double contact_prob = lower + (step * run);
        simulator.SetProbabilities(contact_prob);
        int last_infections = working_population.TotalInfections();

        auto today = start_date;
        while (today < start_date + date::days{kCheckDays}) {
            // Add the newly vaccinated
            if (!input_.vax_history.empty()) {
                simulator.ApplyVaccines(working_population, input_.vax_history.at(input_.state));
            }

            // Simulate the day's new infections and record them
            simulator.SimulateDay(working_population);
            auto new_infections = working_population.TotalInfections() - last_infections;
            results.push_back(new_infections);
            last_infections = working_population.TotalInfections();

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
        printf("%f, %f\n", xs[i], ys[i]);
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

    return {x0, stdev / slope};
}
