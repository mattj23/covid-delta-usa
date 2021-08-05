#pragma once

#include <array>
#include <random>
#include <functional>
#include "data.hpp"

namespace sim {
    class Probabilities {
    public:
        /**
         * Simulates a true or false chance of something happening according to a uniform distribution. If the
         * randomly generated value is less than the probability supplied the function will return true. Thus small
         * probability values have a small chance of occurring, and large ones a large chance of occurring.
         * @param probability a double between 0. and 1. representing the probability of the event occurring
         * @return true if the event "occured", false if not
         */
        inline bool UniformChance(double probability) { return chance_(generator_) <= probability; }

        /**
         * Generates a random, uniformly distributed scalar value that will range between 0 and 1
         * @return double between 0 and 1
         */
        inline double UniformScalar() { return chance_(generator_); }

        inline std::mt19937_64& GetGenerator() { return generator_; }

    private:
        std::mt19937_64 generator_{std::random_device{}()};
        std::uniform_real_distribution<double> chance_{0, 1};
    };

}


