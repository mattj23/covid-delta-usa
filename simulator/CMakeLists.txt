cmake_minimum_required(VERSION 3.16)
project(delta_sim)

set(CMAKE_CXX_STANDARD 20)
find_package(nlohmann_json 3.2.0 REQUIRED)

set(TARGET_SOURCE
        date.h
        sim/probabilities.hpp
        sim/probabilities.cpp
        sim/person.hpp
        sim/person.cpp
        sim/reference_simulators.hpp
        sim/reference_simulators.cpp
        sim/data.hpp
        sim/data.cpp
        sim/simulators.hpp
        sim/simulators.cpp)

add_executable(delta_sim main.cpp ${TARGET_SOURCE})
add_executable(ref_sim reference.cpp ${TARGET_SOURCE})
target_link_libraries(delta_sim PRIVATE nlohmann_json::nlohmann_json)
