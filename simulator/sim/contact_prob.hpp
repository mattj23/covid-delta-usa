#pragma once

#include "covid.hpp"
#include "data.hpp"
#include "simulators.hpp"
#include "variant_probabilities.hpp"
#include "timer.hpp"

namespace sim {
constexpr int kCheckDays = 3;

struct ContactResult {
    double prob;
    double stdev;
};

struct ContactSearchResultSet {
    std::vector<int> days;
    std::vector<double> probabilities;
    std::vector<double> stdevs;
};

void to_json(nlohmann::json &j, const ContactSearchResultSet &o);

class ContactProbabilitySearch {
  public:
    ContactProbabilitySearch(const data::ProgramInput &input, std::shared_ptr<const sim::VariantDictionary> variants);

    ContactResult FindContactProbability(int day);

  private:
    const data::ProgramInput &input_;
    std::shared_ptr<const VariantDictionary> variants_;

    ContactResult GetResultFromBounds(const Population &reference_pop, Population &working_pop,
                                      const std::vector<int> &expected, sim::Simulator &simulator,
                                      date::sys_days start_date, double upper, double lower);
};

} // namespace sim