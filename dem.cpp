#include <cmath>
#include "global_defs.h"
#include "raster.h"
#include "grid_neighbours.h"
#include "dem.h"


DEM::DEM() : Raster(), slope_(), aspect_() {}

DEM::DEM(int size_x_, int size_y_, real_type value) : Raster(size_x_, size_y_, value) {}

DEM::DEM(const std::string &filename) : Raster(filename) {}

/// See http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-slope-works.htm
/// and http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-aspect-works.htm
void DEM::compute_slope_and_aspect(const GridNeighbours& nebs) {
    real_type deltax = get_deltax();
    int size_x = get_size_x();
    int size_y = get_size_y();
    slope_.resize(size_x, size_y);
    aspect_.resize(size_x, size_y);
    for (int i = 0; i < size_x; i++) {
        for (int j = 0; j < size_y; j++) {
            real_type dzdx = ( ( this->operator()(nebs.iup(i), nebs.jdown(j)) +
                        2 * this->operator()(nebs.iup(i), j) + this->operator()(nebs.iup(i), nebs.jup(j)) ) -
                    ( this->operator()(nebs.idown(i), nebs.jdown(j)) + 2 * this->operator()(nebs.idown(i), j) +
                      this->operator()(nebs.idown(i), nebs.jup(j)) ) ) /
                8 / deltax;
            real_type dzdy = ( ( this->operator()(nebs.idown(i), nebs.jup(j)) +
                        2 * this->operator()(i, nebs.jup(j)) + this->operator()(nebs.iup(i), nebs.jup(j)) ) -
                    ( this->operator()(nebs.idown(i), nebs.jdown(j)) + 2 * this->operator()(i, nebs.jdown(j)) +
                      this->operator()(nebs.iup(i), nebs.jdown(j)) ) ) /
                8 / deltax;
            aspect_(i, j) = atan2(dzdy, dzdx);                             // n.b. Aspect in Radians
            slope_(i, j) = sqrt(pow(dzdx, 2) + pow(dzdy, 2));              // n.b. Slope in Radians
        }
    }
}

/// You must call compute_slope_and_aspect before trying to access slope elements.
/// No checking is performed.
real_type DEM::slope(const int i, const int j) const {
    return slope_(i, j);
}

/// You must call compute_slope_and_aspect before trying to access aspect elements
/// No checking is performed.
real_type DEM::aspect(const int i, const int j) const {
    return aspect_(i, j);
}
