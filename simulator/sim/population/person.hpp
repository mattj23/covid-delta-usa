#pragma once

#include <cinttypes>
#include <optional>
#include "../covid.hpp"

namespace sim {

    /** @struct Person
     *
     * @summary This struct is a data-only representation of a single member of a population
     */
    struct Person {
    public:
        /** @summary What variant of SARS-CoV-2 the individual is carrying.
         */
        Variant variant{Variant::None};

        /** @summary On what day was the individual infected with the variant they're currently carrying
         */
        int infected_day{};

        /** @summary On what day does the individual's symptoms manifest
         */
        int symptom_onset{};

        int test_day{};

        /** @summary A scalar value unique to this individual, randomly generated at the time of infection, which
         * when combined with the natural immunity curves determines the binary immunity.
         */
        float natural_immunity_scalar{};

        /** @summary A scalar value unique to this individual, randomly generated at the time of vaccinations, which
         * when combined with the vaccine immunity curve determines the binary immunity.
         */
        float vaccine_immunity_scalar{};

        /** @summary An optional day on which this individual was vaccinated. If the individual is not vaccinated, this
         * optional does not have a value.
         */
        std::optional<int> vaccinated{};

        /** @summary Gets whether or not the individual is carrying a variant
         * @return
         */
        [[nodiscard]] inline bool IsInfected() const { return variant != Variant::None; }

        /** @summary Resets all of the individual's parameters to the defaults
         */
        inline void Reset() {
            variant = Variant::None;
            infected_day = 0;
            symptom_onset = 0;
            test_day = 0;
            natural_immunity_scalar = 0;
            vaccine_immunity_scalar = 0;
            vaccinated = {};
        }
    };

}