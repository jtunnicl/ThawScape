#ifndef _TIMER_H_
#define _TIMER_H_

#include <chrono>

/// \brief A timer object that accumulates time across any number of start-stops
template <typename T>
class AccumulateTimer {
    private:
        bool started;  ///< Stores whether the timer has been started or not
        std::chrono::time_point<std::chrono::system_clock> start_time;  ///< The time the timer was started at
        double total_time;  ///< The total time the timer has been running for

    public:
        /// \brief Create an AccumulateTimer that is not started
        AccumulateTimer() : started(false), total_time(0.0) {}

        /// \brief Start the timer
        void start() {
            if (!started) {
                start_time = std::chrono::system_clock::now();
                started = true;
            }
        }

        /// \brief Stop the timer
        void stop() {
            if (started) {
                std::chrono::time_point<std::chrono::system_clock> stop_time;
                stop_time = std::chrono::system_clock::now();
                total_time += std::chrono::duration_cast<T>(stop_time - start_time).count();
                started = false;
            }
        }

        /// \brief Get the total time the timer has been running for
        ///
        /// This will only be accurate if the timer is stopped. If the timer is running then
        /// it will not include the time spent running in the current open interval.
        double get_total_time() {
            return total_time;
        }
};

#endif
