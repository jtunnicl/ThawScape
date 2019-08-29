#include <iostream>
#include "catch2/catch.hpp"
#include "grid_neighbours.h"

TEST_CASE("GridNeighbours class", "[grid_neighbours]") {
    GridNeighbours nebs;
    int size_x = 20;
    int size_y = 30;

    SECTION("Intialise with constructor") {
        nebs = GridNeighbours(size_x, size_y);
    }

    SECTION("Initialise with setup function") {
        nebs.setup(size_x, size_y);
    }

    // check all returned neighbours are within bounds
    for (int i = 0; i < size_x; i++) {
        int nup = nebs.iup(i);
        REQUIRE(nup >= 0);
        REQUIRE(nup < size_x);
        int ndown = nebs.idown(i);
        REQUIRE(ndown >= 0);
        REQUIRE(ndown < size_x);
    }
    for (int j = 0; j < size_y; j++) {
        int nup = nebs.jup(j);
        REQUIRE(nup >= 0);
        REQUIRE(nup < size_y);
        int ndown = nebs.jdown(j);
        REQUIRE(ndown >= 0);
        REQUIRE(ndown < size_y);
    }
}
