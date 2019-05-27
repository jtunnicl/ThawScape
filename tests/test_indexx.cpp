
#include <iostream>
#include <random>
#include <vector>
#include "catch2/catch.hpp"
#include "indexx.hpp"


TEST_CASE("Indexx class", "[Indexx]") {
    int nx = 5;
    int ny = 3;
    std::vector< std::vector<float> > array = std::vector<std::vector<float>>(nx, std::vector<float>(ny));
    std::cout << "SIZE = " << array.size() << " x " << array[0].size() << std::endl;
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(1, 1000);
    std::cout << "Input:" << std::endl;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            array[i][j] = distribution(generator);
            std::cout << i << ", " << j << " = " << array[i][j] << std::endl;
        }
    }

    Indexx<float> array_indexx(nx, ny);
    array_indexx.update_array(array);

    std::cout << "Output:" << std::endl;
    int ii;
    int jj;
    // indices of first value
    array_indexx.get_ij(0, ii, jj);
    // checking indices are good
    REQUIRE((ii >= 0 && ii < nx));
    REQUIRE((jj >= 0 && jj < ny));
    // first value
    float value1 = array[ii][jj];
    std::cout << ii << ", " << jj << " = " << value1 << std::endl;

    for (int t = 1; t < nx * ny; t++) {
        // indices of next value
        array_indexx.get_ij(t, ii, jj);
        // checking indices are good
        REQUIRE((ii >= 0 && ii < nx));
        REQUIRE((jj >= 0 && jj < ny));
        // next value
        float value2 = array[ii][jj];
        std::cout << ii << ", " << jj << " = " << value2 << std::endl;
        // checking ordering
        REQUIRE(value1 <= value2);
        // moving on
        value1 = value2;
    }
}
