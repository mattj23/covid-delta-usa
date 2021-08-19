#include <sstream>
#include "data.hpp"

void sim::data::from_json(const nlohmann::json &j, sim::data::KnownCaseHistory &d) {
    j.at("total_known_cases").get_to(d.total_known_cases);
}

void sim::data::from_json(const nlohmann::json &j, sim::data::StateInfo &d) {
    j.at("population").get_to(d.population);
    j.at("adjacent").get_to(d.adjacent);
    j.at("ages").get_to(d.ages);
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
    j.at("contact_day_interval").get_to(i.contact_day_interval);
    j.at("population_scale").get_to(i.population_scale);
    j.at("run_count").get_to(i.run_count);
    j.at("state_info").get_to(i.state_info);
    j.at("output_file").get_to(i.output_file);
    j.at("variant_history").get_to(i.variant_history);
    j.at("world_properties").get_to(i.world_properties);
    j.at("options").get_to(i.options);

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

std::unordered_map<sim::Variant, double> sim::data::GetVariantFractions(int date,
                                                                              const std::vector<VariantRecord> &variants) {
    std::unordered_map<sim::Variant, double> results;

    for (const auto &row : variants) {
        if (date <= row.date) {
            // This is the row that's valid for us
            for (const auto &[name, fraction] : row.variants) {
                Variant variant;
                if (name == "alpha") variant = Variant::Alpha;
                if (name == "delta") variant = Variant::Delta;
                results[variant] = fraction;
            }
            return results;
        }
    }

    if (results.empty()) {
        results[Variant::Alpha] = 1.0;
    }

    return results;
}

void sim::data::from_json(const nlohmann::json &j, sim::data::DiscreteFunction &f) {
    j.at("values").get_to(f.values);
    j.at("offset").get_to(f.offset);
}

void sim::data::from_json(const nlohmann::json &j, sim::data::VariantProperties &v) {
    j.at("incubation").get_to(v.incubation);
    j.at("infectivity").get_to(v.infectivity);
    j.at("vax_immunity").get_to(v.vax_immunity);
    j.at("natural_immunity").get_to(v.natural_immunity);
}

void sim::data::from_json(const nlohmann::json &j, sim::data::WorldProperties &w) {
    j.at("alpha").get_to(w.alpha);
    j.at("delta").get_to(w.delta);
}

void sim::data::from_json(const nlohmann::json &j, sim::data::ProgramOptions &o) {
    j.at("full_history").get_to(o.full_history);
    j.at("expensive_stats").get_to(o.expensive_stats);
    j.at("mode").get_to(o.mode);
}


double sim::data::DiscreteFunction::operator()(int day) const {
    auto shifted = day + offset;
    shifted = std::min((int)values.size()-1, shifted);
    shifted = std::max(0, shifted);
    return values[shifted];
}
