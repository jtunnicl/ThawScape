#include <iostream>
#include "time_fcn.h"

#define FIRST_DAY 144  // May 25th
#define LAST_DAY 282   // Oct 10th, freezeup and minimum change after that


time_fcn::time_fcn(int year_, int day_, int hour_, int minute_, int end_year_, int end_day_) {
    year = year_;
    day = day_;
    hour = hour_;
    minute = minute_;
    LocalTime = 0;
    UT = 0;
    end_year = end_year_;
    end_day = end_day_;
}

int time_fcn::get_year() {
    return year;
}

int time_fcn::get_day() {
    return day;
}

int time_fcn::get_hour() {
    return hour;
}

int time_fcn::get_minute() {
    return minute;
}

int time_fcn::get_end_year() {
    return end_year;
}

int time_fcn::get_end_day() {
    return end_day;
}

bool time_fcn::keep_going() {
    return (year <= end_year && day < end_day);
}

void time_fcn::increment(int timestep) {
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

void time_fcn::print() {
    std::cout << "Year: " << year << "; Day: " << day << "; Hour: " << hour << std::endl;
}
