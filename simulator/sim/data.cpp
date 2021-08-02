#include <sstream>
#include "data.hpp"

void sim::data::from_json(const nlohmann::json &j, sim::data::KnownCaseHistory &d) {
    j.at("total_known_cases").get_to(d.total_known_cases);
}

void sim::data::from_json(const nlohmann::json &j, sim::data::StateInfo &d) {
    j.at("population").get_to(d.population);
    j.at("adjacent").get_to(d.adjacent);
}

void sim::data::from_json(const nlohmann::json &j, sim::data::InfectedHistory &d) {
    j.at("total_cases").get_to(d.total_cases);
    j.at("total_infections").get_to(d.total_infections);
}

void sim::data::from_json(const nlohmann::json &j, sim::data::VaccineHistory &d) {
    j.at("total_completed_vax").get_to(d.total_completed_vax);
}

sim::data::ProgramInput sim::data::LoadData(const std::string& file_name) {
    std::ifstream state_file(file_name);
    nlohmann::json j;
    state_file >> j;
    state_file.close();

    return j.get<sim::data::ProgramInput>();
}

void sim::data::to_json(nlohmann::json &j, const sim::data::StepResult &r) {
    j = nlohmann::json{
            {"year", r.year},
            {"month", r.month},
            {"day", r.day},
            {"infections", r.infections},
            {"new_infections", r.new_infections},
            {"positive_tests", r.positive_tests},
            {"has_delta", r.delta},
            {"vaccinated", r.vaccinated},
            {"vaccine_saves", r.vaccine_saves}
    };
}

void sim::data::to_json(nlohmann::json &j, const sim::data::StateResult &r) {
    j = nlohmann::json{
            {"name", r.name},
            {"results", r.results}
    };
}

date::sys_days sim::data::FromString(const std::string& s) {
    int day, month, year;
    sscanf(s.c_str(), "%4d-%2d-%2d", &year, &month, &day);
    return date::year{year}/month/day;
}

void sim::data::from_json(const nlohmann::json &j, sim::data::ProgramInput &i) {
    using std::unordered_map;
    using std::string;
    using date::sys_days;

    auto start_text = j.at("start_day").get<std::string>();
    auto end_text = j.at("end_day").get<std::string>();
    i.start_day = FromString(start_text);
    i.end_day = FromString(end_text);
    j.at("state").get_to(i.state);
    j.at("contact_probability").get_to(i.contact_probability);
    j.at("delta_infectivity_ratio").get_to(i.delta_infectivity_ratio);
    j.at("delta_incubation_ratio").get_to(i.delta_incubation_ratio);
    j.at("population_scale").get_to(i.population_scale);
    j.at("run_count").get_to(i.run_count);
    j.at("state_info").get_to(i.state_info);
    j.at("output_file").get_to(i.output_file);
    j.at("variant_history").get_to(i.variant_history);

    auto infected_history = j.at("infected_history").get<unordered_map<string, unordered_map<string, InfectedHistory>>>();
    auto test_data = j.at("test_history").get<unordered_map<string, unordered_map<string, KnownCaseHistory>>>();
    auto vax_data = j.at("vax_history").get<unordered_map<string, unordered_map<string, VaccineHistory>>>();

    // Convert infection data
    for (const auto &[state, state_data] : infected_history) {
        for (const auto &[date_string, data] : state_data) {
            int d = ToReferenceDate(FromString(date_string));
            i.infected_history[state][d] = data;
        }
    }

    // Convert test history data
    for (const auto &[state, state_data] : test_data) {
        for (const auto &[date_string, data] : state_data) {
            int d = ToReferenceDate(FromString(date_string));
            i.known_case_history[state][d] = data;
        }
    }

    // Convert vax_history data
    for (const auto &[state, state_data] : vax_data) {
        for (const auto &[date_string, data] : state_data) {
            int d = ToReferenceDate(FromString(date_string));
            i.vax_history[state][d] = data;
        }
    }
}

void sim::data::from_json(const nlohmann::json &j, sim::data::VariantRecord &v) {
    j.at("variants").get_to(v.variants);
    v.date = ToReferenceDate(FromString(j.at("date").get<std::string>()));
}

std::unordered_map<sim::data::Variant, double> sim::data::GetVariantFractions(int date,
                                                                              const std::vector<VariantRecord> &variants) {
    std::unordered_map<sim::data::Variant, double> results;
    for (const auto &row : variants) {
        if (row.date < date) {
            // This is the row that's valid for us
            for (const auto &[name, fraction] : row.variants) {
                Variant variant;
                if (name == "alpha") variant = Variant::Alpha;
                if (name == "delta") variant = Variant::Delta;
                results[variant] = fraction;
            }
        }
    }

    if (results.empty()) {
        results[Variant::Alpha] = 1.0;
    }

    return results;
}
