#include <cstdio>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <memory>

#include "sim/data.hpp"
#include "sim/simulators.hpp"

int main(int argc, char **argv) {
    using std::chrono::year;
    using std::chrono::days;

    std::string data_file = (argc > 1) ? argv[1] : "/tmp/input_data.json";
    printf("Covid Simulation\n");
    printf(" * input file: %s\n", data_file.c_str());

    auto input = sim::data::LoadData(data_file);

    sim::Probabilities::SetDeltaIncubation(input.delta_incubation_ratio);
    sim::Probabilities::SetDeltaInfectivity(input.delta_infectivity_ratio);

    const auto &state_info = input.state_info[input.state];

    auto state = std::make_shared<sim::StateSimulator>(state_info.population, input.population_scale);

    // If we have a direct infection history we'll use that to initialize the model, otherwise we will use the
    // positive test history and have to initialize past the start date and reset back to it
    printf(" * initializing from infections data\n");
    state->InitializePopulation(input.infected_history[input.state],
                                input.variant_history,
                                input.start_day);

    if (!input.vax_history.empty()) {
        printf(" * initializing population vaccines\n");
        state->InitializeVaccines(input.vax_history[input.state], input.start_day);
    }

    printf(" * starting simulation\n");
    std::vector<sim::data::StateResult> results;
    for (int run = 0; run < input.run_count; ++run) {
        results.emplace_back();
        results.back().name = input.state;

        state->ResetPopulationTo(input.start_day);
        state->SetProbabilities(input.contact_probability);

        auto today = input.start_day;
        while (today < input.end_day) {
            // Add the newly vaccinated
            if (!input.vax_history.empty()) state->ApplyTodaysVaccines(input.vax_history[input.state]);

            // Simulate the day's events
            state->SimulateDay();

            // Save the results of the run
            std::chrono::year_month_day converted_date = today;
            sim::data::StepResult step{};
            auto ref_date = sim::data::ToReferenceDate(today);
            step.year = (int)converted_date.year();
            step.month = converted_date.month().operator unsigned int();
            step.day = converted_date.day().operator unsigned int();
            step.infections = state->TotalInfected(ref_date);
            step.new_infections = step.infections - state->TotalInfected(ref_date - 1);
            step.positive_tests = state->TestedPositive(ref_date);
            step.vaccinated = state->TotalVaccinated(ref_date);
            step.vaccine_saves = state->VaccineSaves();
            step.delta = state->TotalWithDelta(ref_date);
            results.back().results.push_back(step);

            // Increment the clock
            today++;
        }

    }
    printf(" * simulation complete\n");

    nlohmann::json encoded = results;
    std::ofstream output{input.output_file.c_str()};
    output << encoded << std::endl;
    output.close();

    return 0;
}
