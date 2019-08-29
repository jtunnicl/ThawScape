#include <iostream>
#include "utility.h"
#include "parameters.h"
#include "model_time.h"


ModelTime::ModelTime(int year_, int day_, int hour_, int minute_, int end_year_, int end_day_) {
    year = year_;
    day = day_;
    hour = hour_;
    minute = minute_;
    LocalTime = 0;
    UT = 0;
    end_year = end_year_;
    end_day = end_day_;
}

ModelTime::ModelTime(const Parameters& params) : ModelTime() {
    year = params.get_year();
    day = params.get_day();
    hour = params.get_hour();
    minute = params.get_minute();
    end_year = params.get_end_year();
    end_day = params.get_end_day();
}

int ModelTime::get_year() const {
    return year;
}

int ModelTime::get_day() const {
    return day;
}

int ModelTime::get_hour() const {
    return hour;
}

int ModelTime::get_minute() const {
    return minute;
}

int ModelTime::get_end_year() const {
    return end_year;
}

int ModelTime::get_end_day() const {
    return end_day;
}

bool ModelTime::keep_going() const {
    return (year < end_year || (year == end_year && day < end_day));
}

void ModelTime::increment(int timestep) {
    if (timestep >= 24) {
        Util::Error("Timestep must be less than 24", 1);
    }

    if ( (hour + timestep) <= 24) {
        hour = hour + timestep;
    }
    else {                      // Clicking over 24 hours - day change, possible year change
        if ( day + 1 > LAST_DAY ) {  // Goes to Oct 10th, freezeup and minimum change after that
             year += 1;
             day = FIRST_DAY;         // Reset to May 25th
        }
        else {
            day += 1;
        }
        hour = hour + timestep - 24;
    }
}

void ModelTime::print() const {
    std::cout << "Year: " << year << "; Day: " << day << "; Hour: " << hour << std::endl;
}
