#include <gtest/gtest.h>
#include <random>
#include "../sim/population/population.hpp"


TEST(PopulationTests, InfectiousStressTests) {
    std::mt19937_64 generator{std::random_device{}()};

    sim::Population pop(1000, 1);
    int iterations = 0;
    int infectious = 0;
    while (++iterations < 10000) {

        // Add a random number of infectious people
        auto not_infectious = pop.people.size() - pop.EndOfInfectious();
        std::uniform_int_distribution<size_t> dist_infect(0, (size_t)std::min((int)not_infectious, 100));
        auto to_infect = dist_infect(generator);
        for (size_t i = 0; i < to_infect; ++i) {
            std::uniform_int_distribution<size_t> select(pop.EndOfInfectious(), pop.people.size());
            auto index = select(generator);
            infectious++;
            pop.people[index].variant = sim::Variant::Alpha;
            pop.AddToInfected(index);
        }

        // Remove a random number of infectious people
        std::uniform_int_distribution<size_t> dist_rmv(0, (size_t)std::min((int)pop.EndOfInfectious(), 100));
        auto to_disinfect = dist_rmv(generator);
        for (size_t i = 0; i < to_disinfect; ++i) {
            std::uniform_int_distribution<size_t> select(0, pop.EndOfInfectious() - 1);
            auto index = select(generator);
            infectious--;
            pop.people[index].variant = sim::Variant::None;
            pop.RemoveFromInfected(index);
        }

        // Verify that only infectious people are at a position less than the pointer end
        for (size_t i = 0; i < pop.EndOfInfectious(); ++i) {
            EXPECT_EQ(sim::Variant::Alpha, pop.people[i].variant);
        }

        // Verify that only non-infectious people are at a position greater or equal to the pointer end
        for (size_t i = pop.EndOfInfectious(); i < pop.people.size(); ++i) {
            EXPECT_EQ(sim::Variant::None, pop.people[i].variant);
        }

        // Verify that the count of infectious people matches expectations
        EXPECT_EQ(infectious, pop.CurrentlyInfectious());
    }
}
