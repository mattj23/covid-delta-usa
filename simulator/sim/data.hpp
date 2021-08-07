#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <nlohmann/json.hpp>

#include "../date.h"


namespace sim::data {
    /** @brief Reference date for standard integer date representations
     * @summary All dates in integer form are integer day offsets from a reference date, January 1, 2019. This constant
     * records that date for everywhere that a conversion is performed.
     */
    constexpr date::sys_days kReferenceZeroDate = date::year{2019}/01/01;

    /** @enum Enumeration for the different covid variants, or None for someone who is not infected
     *
     */
    enum class Variant {
        None,
        Alpha,
        Delta
    };

    struct StateInfo {
        int population;
        std::vector<std::string> adjacent;
    };

    struct KnownCaseHistory {
        int total_known_cases;
    };

    struct InfectedHistory {
        int total_infections;
        int total_cases;
    };

    struct VaccineHistory {
        int total_completed_vax;
    };

    void from_json(const nlohmann::json& j, StateInfo &d);
    void from_json(const nlohmann::json& j, KnownCaseHistory &d);
    void from_json(const nlohmann::json& j, InfectedHistory &d);
    void from_json(const nlohmann::json& j, VaccineHistory &d);

    struct StepResult {
        int year;
        uint month;
        uint day;
        int total_infections;
        int total_vaccinated;
        int never_infected;
        int total_delta_infections;
        int total_alpha_infections;
        int reinfections;
        int vaccine_saves;
        int natural_saves;
        int vaccinated_infections;
        int virus_carriers;

        double population_infectiousness;
    };

    struct StateResult {
        std::string name{};
        std::vector<StepResult> results;
    };

    void to_json(nlohmann::json &j, const StateResult& r);
    void to_json(nlohmann::json &j, const StepResult& r);

    struct VariantRecord {
        int date;
        std::unordered_map<std::string, double> variants;
    };

    std::unordered_map<Variant, double> GetVariantFractions(int date, const std::vector<VariantRecord> &variants);

    void from_json(const nlohmann::json &j, VariantRecord &v);

    struct DiscreteFunction {
        std::vector<double> values;
        int offset;

        double operator()(int day) const;
    };

    void from_json(const nlohmann::json &j, DiscreteFunction &f);

    struct VariantProperties {
        std::vector<double> incubation;
        DiscreteFunction infectivity;
        DiscreteFunction vax_immunity;
        DiscreteFunction natural_immunity;
    };

    void from_json(const nlohmann::json &j, VariantProperties &v);

    struct WorldProperties {
        VariantProperties alpha;
        VariantProperties delta;
    };

    void from_json(const nlohmann::json &j, WorldProperties &w);

    struct ProgramOptions {
        bool full_history = false;
        bool expensive_stats = false;
    };

    void from_json(const nlohmann::json &j, ProgramOptions &o);

    struct ProgramInput {
        date::sys_days start_day;
        date::sys_days end_day;
        std::string state;
        std::string output_file;
        double contact_probability;
        int population_scale;
        int run_count;
        ProgramOptions options;
        WorldProperties world_properties;
        std::unordered_map<std::string, std::unordered_map<int, InfectedHistory>> infected_history;
        std::unordered_map<std::string, std::unordered_map<int, KnownCaseHistory>> known_case_history;
        std::unordered_map<std::string, std::unordered_map<int, VaccineHistory>> vax_history;
        std::unordered_map<std::string, StateInfo> state_info;
        std::unordered_map<std::string, std::vector<VariantRecord>> variant_history;
    };

    void from_json(const nlohmann::json &j, ProgramInput &i);

    date::sys_days FromString(const std::string& s);

    ProgramInput LoadData(const std::string& file_name);

    inline int ToReferenceDate(const date::sys_days &day) {
        return (int)(day - kReferenceZeroDate).count();
    }

    inline date::sys_days ToSysDays(int day) {
        return kReferenceZeroDate + date::days{day};
    }

}