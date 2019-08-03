#include <iostream>
#include "catch2/catch.hpp"
#include "model_time.h"


TEST_CASE("ModelTime class", "[model_time]") {
    int year = 2015;
    int day = LAST_DAY - 1;
    int hour = 12;
    int minute = 1;
    int end_year = 2016;
    int end_day = FIRST_DAY + 1;
    ModelTime mytime(year, day, hour, minute, end_year, end_day);
    REQUIRE(year == mytime.get_year()); 
    REQUIRE(day == mytime.get_day()); 
    REQUIRE(hour == mytime.get_hour()); 
    REQUIRE(minute == mytime.get_minute()); 
    REQUIRE(end_year == mytime.get_end_year()); 
    REQUIRE(end_day == mytime.get_end_day()); 
    mytime.print();
    std::cout << "end year is " << end_year << ", end day is " << end_day << std::endl;

    SECTION("keep going during year prior to end_year") {
        REQUIRE(mytime.keep_going());
    }

    SECTION("increment within same day") {
        mytime.increment(1);
        REQUIRE(year == mytime.get_year()); 
        REQUIRE(day == mytime.get_day()); 
        REQUIRE(hour + 1 == mytime.get_hour()); 
        REQUIRE(minute == mytime.get_minute()); 
        REQUIRE(end_year == mytime.get_end_year()); 
        REQUIRE(end_day == mytime.get_end_day()); 
        mytime.print();
    }

    SECTION("roll over to next day") {
        mytime.increment(13);
        REQUIRE(year == mytime.get_year()); 
        REQUIRE(LAST_DAY == mytime.get_day()); 
        REQUIRE(1 == mytime.get_hour()); 
        REQUIRE(minute == mytime.get_minute()); 
        REQUIRE(end_year == mytime.get_end_year()); 
        REQUIRE(end_day == mytime.get_end_day()); 
        mytime.print();
    }

    SECTION("roll over to next year") {
        mytime.increment(15);
        mytime.increment(15);
        mytime.increment(15);
        REQUIRE(end_year == mytime.get_year());
        REQUIRE(FIRST_DAY == mytime.get_day());
        REQUIRE(9 == mytime.get_hour()); 
        REQUIRE(minute == mytime.get_minute()); 
        REQUIRE(end_year == mytime.get_end_year()); 
        REQUIRE(end_day == mytime.get_end_day()); 
        mytime.print();

        SECTION("keep going in end year prior to end day") {
            REQUIRE(mytime.keep_going());
        }
    }

    SECTION("increment to end day") {
        mytime.increment(20);
        mytime.increment(20);
        mytime.increment(20);
        mytime.increment(12);
        REQUIRE(end_year == mytime.get_year());
        REQUIRE(end_day == mytime.get_day());
        REQUIRE(12 == mytime.get_hour()); 
        REQUIRE(minute == mytime.get_minute()); 
        REQUIRE(end_year == mytime.get_end_year()); 
        REQUIRE(end_day == mytime.get_end_day()); 
        mytime.print();

        SECTION("stop on end day") {
            REQUIRE(!mytime.keep_going());
        }
    }

    SECTION("increment beyond end day") {
        mytime.increment(20);
        mytime.increment(20);
        mytime.increment(20);
        mytime.increment(20);
        mytime.increment(16);
        REQUIRE(end_year == mytime.get_year());
        REQUIRE(end_day + 1 == mytime.get_day());
        REQUIRE(12 == mytime.get_hour()); 
        REQUIRE(minute == mytime.get_minute()); 
        REQUIRE(end_year == mytime.get_end_year()); 
        REQUIRE(end_day == mytime.get_end_day()); 
        mytime.print();

        SECTION("stop beyond end day") {
            REQUIRE(!mytime.keep_going());
        }
    }
}
