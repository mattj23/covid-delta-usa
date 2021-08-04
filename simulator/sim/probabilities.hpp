#pragma once

#include <array>
#include <random>
#include <functional>
#include "data.hpp"

namespace sim {
    namespace consts {
        constexpr std::array<double, 20> kAlphaIncubation = { 0.017, 0.135, 0.310, 0.475, 0.608, 0.708, 0.782, 0.836, 0.876, 0.905, 0.927, 0.944, 0.956, 0.965, 0.973, 0.978, 0.983, 0.986, 0.989, 0.991};

        constexpr double kAlphaInfectivity[] = { 0.001, 0.005, 0.015, 0.035, 0.066, 0.103, 0.135, 0.151, 0.145, 0.122, 0.091, 0.060, 0.035, 0.019, 0.009, 0.004, 0.002};
        constexpr int kAlphaInfectivityStart = -7;
        constexpr int kAlphaInfectivityEnd = 9;

        constexpr double kAlphaVaxEfficacy[] = {0.03, 0.07, 0.10, 0.13, 0.17, 0.20, 0.23, 0.27, 0.30, 0.33, 0.37, 0.40, 0.43, 0.47, 0.50, 0.53, 0.56, 0.60, 0.63, 0.66, 0.70, 0.73, 0.76, 0.80, 0.83, 0.86, 0.90, 0.93 };
    }


    class Probabilities {
    public:
        /**
         * Simulates a true or false chance of something happening according to a uniform distribution. If the
         * randomly generated value is less than the probability supplied the function will return true. Thus small
         * probability values have a small chance of occurring, and large ones a large chance of occurring.
         * @param probability a double between 0. and 1. representing the probability of the event occurring
         * @return true if the event "occured", false if not
         */
        bool UniformChance(double probability) {
            return chance_(generator_) <= probability;
        }

        std::mt19937_64& GetGenerator() { return generator_; }

    private:
        std::mt19937_64 generator_{std::random_device{}()};
        std::uniform_real_distribution<double> chance_{0, 1};
    };

}


