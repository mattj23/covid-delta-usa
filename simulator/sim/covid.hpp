#pragma once
#include <nlohmann/json.hpp>
#include "../date.h"

namespace sim {
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

    struct DailySummary {
    public:
        int day;

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

    void to_json(nlohmann::json &j, const DailySummary &r);

}