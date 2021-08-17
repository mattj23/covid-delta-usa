#pragma once

#include "covid.hpp"
#include "data.hpp"
#include "simulators.hpp"
#include "variant_probabilities.hpp"

namespace sim {
    struct ContactResult {
        double prob;
        double stdev;
    };

    class ContactProbabilitySearch {
    public:
        ContactProbabilitySearch(const data::ProgramInput& input, std::shared_ptr<const sim::VariantDictionary> variants);


        ContactResult FindContactProbability(int day);

    private:
        const data::ProgramInput &input_;
        std::shared_ptr<const VariantDictionary> variants_;

    };

}