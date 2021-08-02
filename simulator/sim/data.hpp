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
        ulong infections;
        ulong new_infections;
        ulong positive_tests;
        ulong vaccine_saves;
        ulong vaccinated;
        ulong delta;
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

    struct ProgramInput {
        date::sys_days start_day;
        date::sys_days end_day;
        std::string state;
        std::string output_file;
        double contact_probability;
        double delta_infectivity_ratio;
        double delta_incubation_ratio;
        int population_scale;
        int run_count;
        std::unordered_map<std::string, std::unordered_map<int, InfectedHistory>> infected_history;
        std::unordered_map<std::string, std::unordered_map<int, KnownCaseHistory>> known_case_history;
        std::unordered_map<std::string, std::unordered_map<int, VaccineHistory>> vax_history;
        std::unordered_map<std::string, StateInfo> state_info;
        std::vector<VariantRecord> variant_history;
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