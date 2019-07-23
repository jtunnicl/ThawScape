#include <iostream>
#include "catch2/catch.hpp"
#include "global_defs.h"
#include "raster.h"


TEST_CASE("Raster class", "[raster]") {
    Raster test;
    REQUIRE(test.get_size_x() == 0);
    REQUIRE(test.get_size_y() == 0);

    SECTION("Initialise with value") {
        test = Raster(3, 2, 3.4);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) {
                REQUIRE(test(i, j) == Approx(3.4));
            }
        }
    }

    SECTION("Load Raster from file") {
        test.load("test_raster.asc");

        REQUIRE(test.get_size_x() == 4);
        REQUIRE(test.get_size_y() == 6);
        REQUIRE(test.get_xllcorner() == Approx(-3.1));
        REQUIRE(test.get_yllcorner() == Approx(4.1));
        REQUIRE(test.get_deltax() == Approx(5.0));
        REQUIRE(test.get_nodata() == Approx(-99999.0));

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 6; j++) {
                real_type id = static_cast<real_type>(i);
                real_type jd = static_cast<real_type>(j);
                real_type expected = id + 1 + (jd + 1) / 10;
                REQUIRE(test(i, j) == Approx(expected));
            }
        }
    }

    SECTION("Resize Raster") {
        int nx = 5;
        int ny = 3;
        test.resize(nx, ny);
        REQUIRE(test.get_size_x() == nx);
        REQUIRE(test.get_size_y() == ny);

        SECTION("Setting data") {
            test(0, 0) = 10;
            test(0, 1) = 1;
            test(0, 2) = 15;
            test(1, 0) = 12;
            test(1, 1) = 22;
            test(1, 2) = 11;
            test(2, 0) = 0.5;
            test(2, 1) = 7;
            test(2, 2) = 9;
            test(3, 0) = 29;
            test(3, 1) = 25;
            test(3, 2) = 4;
            test(4, 0) = 6;
            test(4, 1) = 2;
            test(4, 2) = 39;
            REQUIRE(test(0, 0) == Approx(10));
            REQUIRE(test(0, 1) == Approx(1));
            REQUIRE(test(0, 2) == Approx(15));
            REQUIRE(test(1, 0) == Approx(12));
            REQUIRE(test(1, 1) == Approx(22));
            REQUIRE(test(1, 2) == Approx(11));
            REQUIRE(test(2, 0) == Approx(0.5));
            REQUIRE(test(2, 1) == Approx(7));
            REQUIRE(test(2, 2) == Approx(9));
            REQUIRE(test(3, 0) == Approx(29));
            REQUIRE(test(3, 1) == Approx(25));
            REQUIRE(test(3, 2) == Approx(4));
            REQUIRE(test(4, 0) == Approx(6));
            REQUIRE(test(4, 1) == Approx(2));
            REQUIRE(test(4, 2) == Approx(39));

            SECTION("Sorting data") {
                // expected values in order
                std::vector<real_type> expected_vals = {0.5, 1, 2, 4, 6, 7, 9, 10, 11, 12, 15, 22, 25, 29, 39};
                std::vector<int> expected_i      = {2,   0, 4, 3, 4, 2, 2, 0,  1,  1,  0,  1,  3,  3,  4};
                std::vector<int> expected_j      = {0,   1, 1, 2, 0, 1, 2, 0,  2,  0,  2,  1,  1,  0,  2};

                // do the ordering
                test.sort_data();

                // check the result
                int nxy = nx * ny;
                for (int t = 0; t < nxy; t++) {
                    int i, j;

                    // check i, j are correct
                    test.get_sorted_ij(t, i, j);
                    REQUIRE(i == expected_i[t]);
                    REQUIRE(j == expected_j[t]);

                    // check the value is correct
                    real_type value = test(i, j);
                    REQUIRE(value == Approx(expected_vals[t]));
                }
            }
        }
    }
}
