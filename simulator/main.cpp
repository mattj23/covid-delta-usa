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

//    PerfTimer timer;
//    timer.Start();
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

        printf(" > result (%0.2f)\n", result.prob);

        working_day += date::days{std::max(1, input.contact_day_interval)};
    }

//    timer.Stop();
//    printf("[time] simulation = %0.4f s\n", static_cast<double>(search.sim_timer.Elapsed()) / 1.0e6);
//    printf("[time] copy       = %0.4f s\n", static_cast<double>(search.copy_timer.Elapsed()) / 1.0e6);
//    printf("[time] vaccine    = %0.4f s\n", static_cast<double>(search.vax_timer.Elapsed()) / 1.0e6);
//    printf("[time] total      = %0.4f s\n", static_cast<double>(search.total_timer.Elapsed()) / 1.0e6);
//    printf("[time] all total  = %0.4f s\n", static_cast<double>(timer.Elapsed()) / 1.0e6);

    nlohmann::json encoded = results;
    std::ofstream output{input.output_file.c_str()};
    output << encoded << std::endl;
    output.close();
}

void Simulate(const sim::data::ProgramInput &input, std::shared_ptr<const sim::VariantDictionary> variants) {
    auto state_info = input.state_info.at(input.state);
    sim::Simulator simulator(input.options, variants);
    sim::Population reference_population(state_info.population, input.population_scale, state_info.ages);
    sim::Population population = reference_population;
    printf(" * starting simulation (pop=%i at 1:%i scale)\n", reference_population.people.size(), input.population_scale);

    // Initialize the population from the beginning
    PerfTimer timer;
    timer.Start();
    auto init_result = simulator.InitializePopulation(reference_population,
                                                      input.infected_history.at(input.state),
                                                      input.vax_history.at(input.state),
                                                      input.variant_history.at(input.state),
                                                      input.start_day);
    timer.Stop();
    printf(" * initialization took %0.4f s\n", static_cast<double>(timer.Elapsed()) / 1.0e6);

    timer.Reset();
    timer.Start();
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

    timer.Stop();
    printf(" * %i runs in %0.4f s\n", input.run_count, static_cast<double>(timer.Elapsed()) / 1.0e6);

#ifdef PERF_MEASURE
    printf("[time] loop   = %0.4f s\n", static_cast<double>(simulator.loop_timer.Elapsed()) / 1.0e6);
    printf("[time] alloc  = %0.4f s\n", static_cast<double>(simulator.alloc) / 1.0e6);
    printf("[time] remove = %0.4f s\n", static_cast<double>(simulator.remove_timer.Elapsed()) / 1.0e6);
    printf("[time] infect = %0.4f s\n", static_cast<double>(simulator.infect_timer.Elapsed()) / 1.0e6);
#endif

//    nlohmann::json encoded = results;
//    std::ofstream output{input.output_file.c_str()};
//    output << encoded << std::endl;
//    output.close();
}