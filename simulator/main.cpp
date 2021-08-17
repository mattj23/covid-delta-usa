#include <cstdio>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <memory>
#include "date.h"

#include <omp.h>

#include "sim/data.hpp"
#include "sim/variant_probabilities.hpp"
#include "sim/simulators.hpp"
#include "sim/contact_prob.hpp"

using sim::VariantDictionary;

void Simulate(const sim::data::ProgramInput &input, std::shared_ptr<const sim::VariantDictionary> variants);
void FindContactProb(const sim::data::ProgramInput &input, std::shared_ptr<const sim::VariantDictionary> variants);

int main(int argc, char **argv) {
    using sim::Variant;

    std::string data_file = (argc > 1) ? argv[1] : "/tmp/input_data.json";
    printf("Covid Simulation\n");
    printf(" * input file: %s\n", data_file.c_str());

    auto input = sim::data::LoadData(data_file);
    auto variants = std::make_shared<sim::VariantDictionary>();
    (*variants)[Variant::Alpha] = std::make_unique<sim::VariantProbabilities>(input.world_properties.alpha, Variant::Alpha);
    (*variants)[Variant::Delta] = std::make_unique<sim::VariantProbabilities>(input.world_properties.delta, Variant::Delta);

    if (input.options.mode == sim::data::ProgramMode::Simulate) {
        Simulate(input, variants);
    } else if (input.options.mode == sim::data::ProgramMode::FindContactProb) {
        FindContactProb(input, variants);
    }

    return 0;
}

void FindContactProb(const sim::data::ProgramInput &input, std::shared_ptr<const sim::VariantDictionary> variants) {
    sim::ContactProbabilitySearch search(input, variants);
    sim::ContactSearchResultSet results;
    printf(" * finding contact probabilities\n");

    auto working_day = input.start_day;
    while (working_day <= input.end_day) {
        date::year_month_day ymd = working_day;
        printf(" * contact prob for %i-%u-%u\n",
               (int)ymd.year(),
               ymd.month().operator unsigned int(),
               ymd.day().operator unsigned int());
        auto ref_day = sim::data::ToReferenceDate(working_day);
        auto result = search.FindContactProbability(ref_day);
        results.days.push_back(ref_day);
        results.probabilities.push_back(result.prob);
        results.stdevs.push_back(result.stdev);

        working_day += date::days{std::max(1, input.contact_day_interval)};
    }

    nlohmann::json encoded = results;
    std::ofstream output{input.output_file.c_str()};
    output << encoded << std::endl;
    output.close();
}

void Simulate(const sim::data::ProgramInput &input, std::shared_ptr<const sim::VariantDictionary> variants) {
    auto state_info = input.state_info.at(input.state);
    sim::Simulator simulator(input.options, variants);
    sim::Population reference_population(state_info.population, input.population_scale);
    sim::Population population(state_info.population, input.population_scale);


    // Initialize the population from the beginning
    auto init_result = simulator.InitializePopulation(reference_population,
                                                      input.infected_history.at(input.state),
                                                      input.vax_history.at(input.state),
                                                      input.variant_history.at(input.state),
                                                      input.start_day);

    printf(" * starting simulation\n");

    auto start = std::chrono::system_clock::now();
    std::vector<sim::data::StateResult> results;
    for (int run = 0; run < input.run_count; ++run) {
        population.CopyFrom(reference_population);
        results.emplace_back();
        results.back().name = input.state;

        if (!init_result.empty()) {
            // If the option for exporting the full history is on, we copy the data from the initialization phase into
            // the results storage
            for (const auto& r : init_result) {
                results.back().results.push_back(r);
            }
        } else {
            // If the option for exporting the full history is off, we at least need to export the day before the first
            // simulation day so that differentiated statistics can be computed
            results.back().results.push_back(simulator.GetDailySummary(population, input.options.expensive_stats));
        }

        // Setting the contact probability
        simulator.SetProbabilities(input.contact_probability);

        auto today = input.start_day;
        while (today < input.end_day) {
            // Add the newly vaccinated
            if (!input.vax_history.empty()) simulator.ApplyVaccines(population,
                                                                     input.vax_history.at(input.state));

            // Simulate the day's events
            results.back().results.push_back(simulator.SimulateDay(population));

            // Increment the clock
            today += date::days{1};
        }
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    printf(" * simulation complete in %lu ms\n", elapsed);

    nlohmann::json encoded = results;
    std::ofstream output{input.output_file.c_str()};
    output << encoded << std::endl;
    output.close();
}