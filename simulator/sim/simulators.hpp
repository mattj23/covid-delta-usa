#pragma once
#include "data.hpp"
#include "timer.hpp"
#include "population/person.hpp"
#include "population/population.hpp"
#include "probabilities.hpp"
#include "variant_probabilities.hpp"
#include <limits>
#include <memory>
#include <optional>
#include <unordered_set>
#include <vector>

//#define PERF_MEASURE

namespace sim {

class Simulator {
  public:
    Simulator(const data::ProgramOptions &options, std::shared_ptr<const VariantDictionary> variants);

    [[nodiscard]] DailySummary GetDailySummary(const sim::Population &population, bool expensive) const;

    std::vector<DailySummary> InitializePopulation(sim::Population &population,
                                                   const std::unordered_map<int, data::InfectedHistory> &history,
                                                   const std::unordered_map<int, data::VaccineHistory> &vaccines,
                                                   const std::vector<data::VariantRecord> &variant_history,
                                                   std::optional<date::sys_days> up_to = {});

    void InfectPerson(sim::Population &population, size_t person_index, const VariantProbabilities &variant);

    void ApplyVaccines(sim::Population &population, const std::unordered_map<int, data::VaccineHistory> &vaccines);

    inline void SetProbabilities(double p_self) { contact_probability_ = p_self; }

    DailySummary SimulateDay(sim::Population &population);

#ifdef PERF_MEASURE
    PerfTimer loop_timer;
    PerfTimer remove_timer;
    PerfTimer infect_timer;
    long alloc{};
#endif

  private:
    double contact_probability_{};
    std::shared_ptr<const VariantDictionary> variants_;
    data::ProgramOptions options_;

    Probabilities prob_{};
};


} // namespace sim