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

int main(int argc, char **argv) {
    using sim::data::Variant;

    std::string data_file = (argc > 1) ? argv[1] : "/tmp/input_data.json";
    printf("Covid Simulation\n");
    printf(" * input file: %s\n", data_file.c_str());

    auto input = sim::data::LoadData(data_file);

    std::unordered_map<Variant, std::unique_ptr<sim::VariantProbabilities>> variants;
    variants[Variant::Alpha] = std::make_unique<sim::VariantProbabilities>(input.world_properties.alpha, Variant::Alpha);
    variants[Variant::Delta] = std::make_unique<sim::VariantProbabilities>(input.world_properties.delta, Variant::Delta);

    const auto &state_info = input.state_info[input.state];

    auto state = std::make_shared<sim::StateSimulator>(state_info.population, input.population_scale, &variants);

    auto start = std::chrono::system_clock::now();

    printf(" * starting simulation\n");
    std::vector<sim::data::StateResult> results;
    for (int run = 0; run < input.run_count; ++run) {
        // Initialize the population from the beginning
        state->InitializePopulation(input.infected_history[input.state],
                                    input.variant_history[input.state],
                                    input.start_day);

        if (!input.vax_history.empty()) {
            state->InitializeVaccines(input.vax_history[input.state], input.start_day);
        }

        state->SetProbabilities(input.contact_probability);

        results.emplace_back();
        results.back().name = input.state;

        auto today = input.start_day;
        while (today < input.end_day) {
            results.back().results.push_back(state->GetStepResult());

            // Add the newly vaccinated
            if (!input.vax_history.empty()) state->ApplyTodaysVaccines(input.vax_history[input.state]);

            // Simulate the day's events
            state->SimulateDay();

            // Save the results of the run

            // Increment the clock
            today += date::days{1};
        }
        results.back().results.push_back(state->GetStepResult());
    }

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    printf(" * simulation complete in %lu ms\n", elapsed);

    nlohmann::json encoded = results;
    std::ofstream output{input.output_file.c_str()};
    output << encoded << std::endl;
    output.close();

    return 0;
}
