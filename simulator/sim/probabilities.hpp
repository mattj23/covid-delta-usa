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
        static double GetAlphaInfectivity(int days_from_symptoms);
        static double GetAlphaVaxEfficacy(int days_from_vax);

        static double GetInfectivity(data::Variant variant, int days_from_symptoms);

        int GetIncubation(data::Variant variant);

        int GetAlphaIncubation();

        /** @brief Generate a random test lag related to the alpha variant incubation time
         *
         * @param factor a value to scale the randomly sampled alpha incubation time for test lag
         * @return a number of days which represent the lag between the date of infection and the date of the test
         */
        int GetTestingLag(double factor=1.5);

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

        inline static void SetDeltaIncubation(double ratio) { delta_incubation_ratio = ratio; }
        inline static void SetDeltaInfectivity(double ratio) { delta_infectivity_ratio = ratio; }

    private:
        std::random_device rd_{};
        std::mt19937_64 generator_{rd_()};
        std::uniform_real_distribution<double> chance_{0, 1};
        std::weibull_distribution<double> alpha_incubation_{2.04, 1.89};

        inline static double delta_incubation_ratio;
        inline static double delta_infectivity_ratio;
    };

}


