cmake_minimum_required(VERSION 3.16)
project(delta_sim)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(OpenMP REQUIRED)

set(INSTALL_GTEST OFF)
add_subdirectory(googletest)

set(TARGET_SOURCE
        date.h
        sim/timer.hpp
        sim/covid.hpp
        sim/covid.cpp
        sim/contact_prob.hpp
        sim/contact_prob.cpp
        sim/probabilities.hpp
        sim/probabilities.cpp
        sim/variant_probabilities.hpp
        sim/variant_probabilities.cpp
        sim/population/person.hpp
        sim/population/person.cpp
        sim/population/population.hpp
        sim/population/population.cpp
        sim/data.hpp
        sim/data.cpp
        sim/simulators.hpp
        sim/simulators.cpp)

add_executable(delta_sim main.cpp ${TARGET_SOURCE})
target_link_libraries(delta_sim PRIVATE nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

add_executable(gtest_run tests/population_tests.cpp ${TARGET_SOURCE})
target_link_libraries(gtest_run PRIVATE gtest gtest_main nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)
