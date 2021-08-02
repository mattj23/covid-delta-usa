#pragma once

#include <cinttypes>
#include <optional>

namespace sim {

    struct Person {
    public:
        data::Variant variant{data::Variant::None};
        bool is_immune{};
        int infected_day{};
        int symptom_onset{};
        int test_day{};
        std::optional<int> vaccinated{};

        inline bool IsInfected() const { return variant != data::Variant::None; }
    };

}