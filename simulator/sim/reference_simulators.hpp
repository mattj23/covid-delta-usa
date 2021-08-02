#pragma once

#include <vector>
#include <unordered_set>
#include <algorithm>
#include "probabilities.hpp"
#include "person.hpp"

namespace sim {
    class MethodBase {
    public:
        explicit MethodBase(long population);

        virtual void Run(int initial_infected, int days, double contact_prob) = 0;

    protected:
        Probabilities prob_{};
        std::vector<Person> pop_{};

        [[nodiscard]] long TotalInfected() const;
    };

    /**
     * This is an extremely naive, blatantly simple and (hopefully) obviously correct implementation of a covid
     * alpha variant spread in an isolated population.  It's meant to be the reference implementation to check
     * optimized versions against.
     */
    class NaiveMethod : public MethodBase {
    public:
        explicit NaiveMethod(long population);

        void Run(int initial_infected, int days, double contact_prob) override;
    };

    class OptimizedMethod : public MethodBase {
    public:
        explicit OptimizedMethod(long population);

        void Run(int initial_infected, int days, double contact_prob) override;
    };
}