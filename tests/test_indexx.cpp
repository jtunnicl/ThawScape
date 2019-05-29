#include <iostream>
#include <random>
#include <vector>
#include "catch2/catch.hpp"
#include "indexx.hpp"


TEST_CASE("Indexx class", "[indexx]") {
    // input array
    int nx = 5;
    int ny = 3;
    std::vector< std::vector<float> > array = std::vector<std::vector<float>>(nx, std::vector<float>(ny));
    array[0][0] = 10;
    array[0][1] = 1;
    array[0][2] = 15;
    array[1][0] = 12;
    array[1][1] = 22;
    array[1][2] = 11;
    array[2][0] = 0.5;
    array[2][1] = 7;
    array[2][2] = 9;
    array[3][0] = 29;
    array[3][1] = 25;
    array[3][2] = 4;
    array[4][0] = 6;
    array[4][1] = 2;
    array[4][2] = 39;

    // expected values in order
    std::vector<float> expected_vals = {0.5, 1, 2, 4, 6, 7, 9, 10, 11, 12, 15, 22, 25, 29, 39};
    std::vector<int> expected_i      = {2,   0, 4, 3, 4, 2, 2, 0,  1,  1,  0,  1,  3,  3,  4};
    std::vector<int> expected_j      = {0,   1, 1, 2, 0, 1, 2, 0,  2,  0,  2,  1,  1,  0,  2};

    // do the ordering
    Indexx<float> array_indexx(nx, ny);
    array_indexx.update_array(array);

    // check the result
    int nxy = nx * ny;
    for (int t = 0; t < nxy; t++) {
        int i, j;

        // check i, j are correct
        array_indexx.get_ij(t, i, j);
        REQUIRE(i == expected_i[t]);
        REQUIRE(j == expected_j[t]);

        // check the value is correct
        float value = array[i][j];
        REQUIRE(value == expected_vals[t]);
    }
}
