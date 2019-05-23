
#include <iostream>
#include <random>
#include <vector>
#include "../indexx.hpp"


int main() {
    int nx = 4;
    int ny = 2;
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

    Indexx<float> array_indexx(nx, ny, array);

    std::cout << "Output:" << std::endl;
    int ii;
    int jj;
    // indices of first value
    array_indexx.get_ij(0, ii, jj);
    // checking indices are good
    if (ii < 0 || ii >= nx) {
        std::cerr << "Index error for i" << std::endl;
        return 1;
    }
    if (jj < 0 || jj > ny) {
        std::cerr << "Index error for j" << std::endl;
        return 1;
    }
    // first value
    float value1 = array[ii][jj];
    std::cout << ii << ", " << jj << " = " << value1 << std::endl;

    for (int t = 1; t < nx * ny; t++) {
        // indices of next value
        array_indexx.get_ij(t, ii, jj);
        // checking indices are good
        if (ii < 0 || ii >= nx) {
            std::cerr << "Index error for i" << std::endl;
            return 1;
        }
        if (jj < 0 || jj > ny) {
            std::cerr << "Index error for j" << std::endl;
            return 1;
        }
        // next value
        float value2 = array[ii][jj];
        std::cout << ii << ", " << jj << " = " << value2 << std::endl;
        // checking ordering
        if (value1 > value2) {
            std::cerr << "Values not ordered correctly (" << value1 << ", " << value2 << ")" << std::endl;
            return 1;
        }
        // moving on
        value1 = value2;
    }

    std::cout << "All good" << std::endl;

    return 0;
}
