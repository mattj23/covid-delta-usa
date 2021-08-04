#pragma once
#include "data.hpp"
#include "person.hpp"
#include "probabilities.hpp"
#include "variant_probabilities.hpp"
#include <limits>
#include <optional>
#include <unordered_set>
#include <vector>

namespace sim {
class StateSimulator {
  public:
    StateSimulator(long population, int scale);

    /** @brief Reset all of the population data, bringing them back to their state at a particular date
     *
     * The reset date is inclusive, meaning infections on the reset date itself will be reset. The reset date thus
     * becomes the first date ready for a simulation.
     */
    void ResetPopulationTo(date::sys_days date);

    /** @brief Initialize the population from a dataset of assumed daily infections
     *
     * @param history
     * @param up_to
     */
    void InitializePopulation(const std::unordered_map<int, data::InfectedHistory> &history,
                              const std::vector<data::VariantRecord> &variant_history,
                              const std::unordered_map<data::Variant, VariantProbabilities> &variants,
                              std::optional<date::sys_days> up_to = {});

    void InitializeVaccines(const std::unordered_map<int, data::VaccineHistory> &vaccines,
                            std::optional<date::sys_days> up_to = {});

    void SetProbabilities(double p_self);

    void SimulateDay(const std::unordered_map<data::Variant, VariantProbabilities> &variants);

    void ApplyVaccines(const std::unordered_map<int, data::VaccineHistory> &vaccines, int for_date);
    void ApplyTodaysVaccines(const std::unordered_map<int, data::VaccineHistory> &vaccines);

    size_t TotalInfectious();
    size_t TestedPositive(int on_day);
    size_t TotalInfected(int on_day);
    size_t TotalVaccinated(int on_day);
    size_t TotalWithDelta(int on_day);

    inline int VaccineSaves() const { return vaccine_saves_; }

  private:
    int scale_{};
    int today_{};
    int vaccine_saves_{};
    Probabilities prob_{};
    std::vector<Person> pop_{};
    std::unordered_set<size_t> infectious_{};
    std::vector<size_t> un_vaxxed_{};

    std::unique_ptr<std::binomial_distribution<int>> self_contact_dist_{};
    std::unique_ptr<std::uniform_int_distribution<int>> selector_dist_{};

    /** @brief clear the infectious_ cache set and re-fill it by going through the entire population
     *
     */
    void SynchronizeInfectedCache();
    void SynchronizeUnVaxxedCache();
};

} // namespace sim