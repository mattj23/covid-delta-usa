#pragma once

#include <cinttypes>
#include <optional>

namespace sim {

    struct Person {
    public:
        data::Variant variant{data::Variant::None};
        int infected_day{};
        int symptom_onset{};
        int test_day{};
        float natural_immunity_scalar{};
        float vaccine_immunity_scalar{};
        std::optional<int> vaccinated{};

        [[nodiscard]] inline bool IsInfected() const { return variant != data::Variant::None; }
    };

}