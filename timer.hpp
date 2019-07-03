#ifndef _TIMER_H_
#define _TIMER_H_

#include <chrono>

template <typename T>
class AccumulateTimer {
    private:
        bool started;
        std::chrono::time_point<std::chrono::system_clock> start_time;
        double total_time;

    public:
        AccumulateTimer() : started(false), total_time(0.0) {}

        void start() {
            if (!started) {
                start_time = std::chrono::system_clock::now();
                started = true;
            }
        }

        void stop() {
            if (started) {
                std::chrono::time_point<std::chrono::system_clock> stop_time;
                stop_time = std::chrono::system_clock::now();
                total_time += std::chrono::duration_cast<T>(stop_time - start_time).count();
                started = false;
            }
        }

        double get_total_time() {
            return total_time;
        }
};

#endif
