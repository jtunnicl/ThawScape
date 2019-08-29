#include <cmath>
#include <string>
#include <iostream>
#include "raster.h"
#include "global_defs.h"


int main(int argc, char** argv) {
    // arguments
    if (argc < 3) {
        std::cerr << "Usage: compare2rasters <raster_file_1> <raster_file_2> [tolerance]" << std::endl;
        return 1;
    }
    std::string file1 = std::string(argv[1]);
    std::string file2 = std::string(argv[2]);
    real_type tolerance = 0.1;
    if (argc > 3) {
        tolerance = std::stoi(argv[3]);
    }
    std::cout << "Comparing rasters: " << file1 << " and " << file2 << " with tolerance " << tolerance << std::endl;

    // load rasters
    Raster raster1(file1);
    Raster raster2(file2);

    // check they are the same size
    if (raster1.get_size_x() != raster2.get_size_x()) {
        std::cerr << "x dimensions differ" << std::endl;
        return 1;
    }
    if (raster1.get_size_y() != raster2.get_size_y()) {
        std::cerr << "y dimensions differ" << std::endl;
        return 1;
    }

    // calculate Euclidean distance
    real_type sum = 0.0;
    for (int i = 0; i < raster1.get_size_x(); i++) {
        for (int j = 0; j < raster1.get_size_y(); j++) {
            real_type diff = raster1(i, j) - raster2(i, j);
            sum += diff * diff;
        }
    }
    real_type separation = sqrt(sum);
    std::cout << "Separation between rasters: " << separation << std::endl;

    // check separation against tolerance
    if (separation > tolerance) {
        std::cerr << "Rasters differ" << std::endl;
        return 1;
    }
    else {
        std::cout << "Rasters are the same" << std::endl;
    }

    return 0;
}
