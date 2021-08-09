#pragma once
#include "data.hpp"
#include "population/person.hpp"
#include "population/population.hpp"
#include "probabilities.hpp"
#include "variant_probabilities.hpp"
#include <limits>
#include <memory>
#include <optional>
#include <unordered_set>
#include <vector>

namespace sim {
using VariantDictionary = std::unordered_map<Variant, std::unique_ptr<VariantProbabilities>>;

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

  private:
    double contact_probability_{};
    std::shared_ptr<const VariantDictionary> variants_;
    data::ProgramOptions options_;

    Probabilities prob_{};
};

//class StateSimulator {
//  public:
//    StateSimulator(long population, int scale, const VariantDictionary *variants);
//
//    /** @brief Reset all of the population data, bringing them back to their state pre-pandemic
//     */
//    void Reset();
//
//    /** @brief Initialize the population from a dataset of assumed daily infections
//     *
//     * @param history
//     * @param up_to
//     */
//    std::vector<data::StepResult> InitializePopulation(const std::unordered_map<int, data::InfectedHistory> &history,
//                                                       const std::unordered_map<int, data::VaccineHistory> &vaccines,
//                                                       const std::vector<data::VariantRecord> &variant_history,
//                                                       std::optional<date::sys_days> up_to = {});
//
//    void SetProbabilities(double p_self);
//    inline void SetOptions(const data::ProgramOptions *options) { options_ = options; }
//
//    data::StepResult SimulateDay();
//
//    void ApplyVaccines(const std::unordered_map<int, data::VaccineHistory> &vaccines);
//
//    inline int TotalInfections() const { return total_infections_ * scale_; }
//    inline int TotalVaccinated() const { return total_vaccinated_ * scale_; }
//    inline int NeverInfected() const { return never_infected_ * scale_; }
//    inline int TotalDeltaInfections() const { return total_delta_infections_ * scale_; }
//    inline int TotalAlphaInfections() const { return total_alpha_infections_ * scale_; }
//    inline int Reinfections() const { return reinfections_ * scale_; }
//
//    inline int VaccineSaves() const { return vaccine_saves_ * scale_; }
//    inline int NaturalSaves() const { return natural_saves_ * scale_; }
//
//    data::StepResult GetStepResult() const;
//
//    void Reseed();
//
//  private:
//    int scale_{};
//    int today_{};
//    int vaccine_saves_{};
//    int natural_saves_{};
//
//    int total_infections_{};
//    int total_vaccinated_{};
//    int never_infected_{};
//    int total_delta_infections_{};
//    int total_alpha_infections_{};
//    int reinfections_{};
//    int vaccinated_infections_{};
//
//    Probabilities prob_{};
//    std::vector<Person> pop_{};
//    std::unordered_set<size_t> infectious_{};
//    std::vector<size_t> un_vaxxed_{};
//
//    std::binomial_distribution<int> self_contact_dist_{};
//    std::uniform_int_distribution<int> selector_dist_{};
//
//    const VariantDictionary *variants_{};
//    const data::ProgramOptions *options_{};
//
//    void SynchronizeUnVaxxedCache();
//    void InfectPerson(size_t person_index, const VariantProbabilities &variant);
//};

} // namespace sim