#pragma once
#include <chrono>

class PerfTimer {
public:
    inline void Start() {
        if (is_running_) return;
        start_ = std::chrono::system_clock::now();
        is_running_ = true;
    }

    inline void Stop() {
        if (!is_running_) return;
        auto end = std::chrono::system_clock::now();
        elapsed_ += std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
        is_running_ = false;
    }

    inline long Elapsed() {
        if (is_running_) Stop();
        return elapsed_;
    }

    inline void Reset() {
        is_running_ = false;
        elapsed_ = 0;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> start_{};
    bool is_running_{};
    long elapsed_{};
};