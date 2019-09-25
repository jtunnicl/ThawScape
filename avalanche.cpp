#include <algorithm>
#include "dem.h"
#include "raster.h"
#include "grid_neighbours.h"
#include "avalanche.h"


Avalanche::Avalanche(DEM& topo_, Raster& sed_track_, GridNeighbours& nebs_) :
        topo(topo_), sed_track(sed_track_), nebs(nebs_), initialised(false) {
    size_x = topo.get_size_x();
    size_y = topo.get_size_y();
}

void Avalanche::initialise() {
    size_x = topo.get_size_x();
    size_y = topo.get_size_y();

	// Code in Init subroutine:
	//	thresh = 0.577 * deltax;   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
	//  thresh_diag = thresh * sqrt2;
	thresh = 0.577 * topo.get_deltax();   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
	thresh_diag = thresh * sqrt2;

    initialised = true;
}

void Avalanche::run() {
    if (!initialised) {
        initialise();
    }

    // sort by elevation
    topo.sort_data();

	// NEED TO ASSESS WHETHER PIXEL HAS SEDIMENT, BEFORE FAILURE CALCS?

    int t = 0;
    // Landsliding, proceeding from high elev to low
    // TODO: actually we are proceeding from low to high here
    while (t < size_x * size_y)
    {
        int i, j;
        topo.get_sorted_ij(t, i, j);
        real_type clifftop = 0;

        if (topo(nebs.iup(i), j) - topo(i, j) > thresh) {
            clifftop = topo(nebs.iup(i), j);    // Height of overhanging pixel
            topo(nebs.iup(i), j) = std::max((topo(i, j) + thresh), (topo(nebs.iup(i), j) - sed_track(nebs.iup(i), j)));
        }

        if (topo(nebs.idown(i), j) - topo(i, j) > thresh)
            topo(nebs.idown(i), j) = topo(i, j) + thresh;
        if (topo(i, nebs.jup(j)) - topo(i, j) > thresh)
            topo(i, nebs.jup(j)) = topo(i, j) + thresh;
        if (topo(i, nebs.jdown(j)) - topo(i, j) > thresh)
            topo(i, nebs.jdown(j)) = topo(i, j) + thresh;
        if (topo(nebs.iup(i), nebs.jup(j)) - topo(i, j) > (thresh_diag))
            topo(nebs.iup(i), nebs.jup(j)) = topo(i, j) + thresh_diag;
        if (topo(nebs.iup(i), nebs.jdown(j)) - topo(i, j) > (thresh_diag))
            topo(nebs.iup(i), nebs.jdown(j)) = topo(i, j) + thresh_diag;

        if (topo(nebs.idown(i), nebs.jup(j)) - topo(i, j) > (thresh_diag))
            topo(nebs.idown(i), nebs.jup(j)) = topo(i, j) + thresh_diag;
        if (topo(nebs.idown(i), nebs.jdown(j)) - topo(i, j) > (thresh_diag))
            topo(nebs.idown(i), nebs.jdown(j)) = topo(i, j) + thresh_diag;

        t++;
    }
}
