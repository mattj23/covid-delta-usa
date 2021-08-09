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
    auto state_info = input.state_info.at(input.state);
    sim::Simulator simulator(input.options, variants);
    sim::Population reference_population(state_info.population, input.population_scale);

    printf(" * starting contact probability simulation\n");

    // The starting guess for the contact probability is the value that was supplied
    constexpr int kCheckDays = 3;

    const auto& infected_history = input.infected_history.at(input.state);
    std::vector<int> expected;
    for (int i = 0; i < kCheckDays; ++i) {
        int i0 = infected_history.at(sim::data::ToReferenceDate(input.start_day + date::days{i-1})).total_infections;
        int i1 = infected_history.at(sim::data::ToReferenceDate(input.start_day + date::days{i})).total_infections;
        expected.push_back(i1 - i0);
    }

    // Initialize the population from the beginning
    auto init_start = std::chrono::system_clock::now();
    auto init_result = simulator.InitializePopulation(reference_population, infected_history,
                                                   input.vax_history.at(input.state),
                                                   input.variant_history.at(input.state),
                                                   input.start_day);
    auto init_end = std::chrono::system_clock::now();
    auto init_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(init_end - init_start).count();
    printf(" * init in %lu ms\n", init_elapsed);


    auto start = std::chrono::system_clock::now();
    std::vector<double> errors;
    for (int run = 0; run < input.run_count; ++run) {
        std::vector<int> results;
        sim::Population working_population(reference_population);

        // Setting the contact probability
        simulator.SetProbabilities(input.contact_probability);
        int last_infections = working_population.TotalInfections();

        auto today = input.start_day;
        while (today < input.start_day + date::days{kCheckDays}) {
            // Add the newly vaccinated
            if (!input.vax_history.empty()) {
                simulator.ApplyVaccines(working_population, input.vax_history.at(input.state));
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
        errors.push_back(error);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    printf(" * contact %.3f simulated in %lu ms\n", input.contact_probability, elapsed);

    nlohmann::json encoded = errors;
    std::ofstream output{input.output_file.c_str()};
    output << encoded << std::endl;
    output.close();
}

void Simulate(const sim::data::ProgramInput &input, std::shared_ptr<const sim::VariantDictionary> variants) {
    auto state_info = input.state_info.at(input.state);
    sim::Simulator simulator(input.options, variants);
    sim::Population reference_population(state_info.population, input.population_scale);

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
        auto population = std::make_unique<sim::Population>(reference_population);
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
            results.back().results.push_back(simulator.GetDailySummary(*population, input.options.expensive_stats));
        }

        // Setting the contact probability
        simulator.SetProbabilities(input.contact_probability);

        auto today = input.start_day;
        while (today < input.end_day) {
            // Add the newly vaccinated
            if (!input.vax_history.empty()) simulator.ApplyVaccines(*population,
                                                                     input.vax_history.at(input.state));

            // Simulate the day's events
            results.back().results.push_back(simulator.SimulateDay(*population));

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