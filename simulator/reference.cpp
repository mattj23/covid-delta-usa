#include <iostream>
#include <cstdio>
#include "sim/reference_simulators.hpp"
#include <chrono>
#include <thread>

int main(int argc, char **argv) {
    using namespace std::chrono_literals;

    auto start = std::chrono::system_clock::now();

    sim::OptimizedMethod optimized{100000};
    optimized.Run(100, 100, 1.75);

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//    std::cout << elapsed.count() << " ms" << std::endl;

//    std::this_thread::sleep_for(100ms);
}

