#include <iostream>
#include "priority_flood.hpp"
#include "parameters.h"
#include "flood.h"

#define fillincrement 0.01

/// References to the input Rasters and vectors are stored since they are expected to be
/// managed elsewhere. One must call the inititialse function before running the flow
/// routing.
Flood::Flood(DEM& topo_, GridNeighbours& nebs_) :
        topo(topo_), nebs(nebs_), initialised(false), algorithm(0) {
    size_x = topo.get_size_x();
    size_y = topo.get_size_y();
}

void Flood::initialise(Parameters& params) {
    size_x = topo.get_size_x();
    size_y = topo.get_size_y();
    algorithm = params.get_flood_algorithm();

    if (algorithm == 1) {
        std::cout << "<Flood>: using Barnes' original_priority_flood" << std::endl;
        elevation = Array2D<real_type>(size_x, size_y, -9999.0);
    }
    else if (algorithm == 2) {
        std::cout << "<Flood>: using Barnes' priority_flood_epsilon" << std::endl;
        elevation = Array2D<real_type>(size_x, size_y, -9999.0);
    }
    else {
        std::cout << "<Flood>: using Pelletier's fillinpitsandflats" << std::endl;
        elevation = Array2D<real_type>();
    }

    initialised = true;
}

void Flood::run() {
    if ((algorithm == 1) || (algorithm == 2)) {
        run_barnes_flood();
    }
    else {
        run_fillinpitsandflats();
    }
}

void Flood::run_barnes_flood() {
	// update elev
	for (int i = 0; i < size_x; i++)
	{
		for (int j = 0; j < size_y; j++)
		{
			elevation(i, j) = topo(i, j);     //  Change indexing to suit flood subroutine.
		}
	}

	// perform flooding
    if (algorithm == 1) {
        original_priority_flood(elevation);
    }
    else if (algorithm == 2) {
        priority_flood_epsilon(elevation);
    }

	// update topo
	for (int i = 0; i < size_x; i++)
	{
		for (int j = 0; j < size_y; j++)
		{
			topo(i, j) = elevation(i, j);     //  Back to original indexing
		}
	}
}

void Flood::fillinpitsandflats(int i, int j) {
    real_type minv = topo(i, j);

    if (topo(nebs.iup(i), j) < minv) minv = topo(nebs.iup(i), j);
    if (topo(nebs.idown(i), j) < minv) minv = topo(nebs.idown(i), j);
    if (topo(i, nebs.jup(j)) < minv) minv = topo(i, nebs.jup(j));
    if (topo(i, nebs.jdown(j)) < minv) minv = topo(i, nebs.jdown(j));
    if (topo(nebs.iup(i), nebs.jup(j)) < minv) minv = topo(nebs.iup(i), nebs.jup(j));
    if (topo(nebs.idown(i), nebs.jup(j)) < minv) minv = topo(nebs.idown(i), nebs.jup(j));
    if (topo(nebs.idown(i), nebs.jdown(j)) < minv) minv = topo(nebs.idown(i), nebs.jdown(j));
    if (topo(nebs.iup(i), nebs.jdown(j)) < minv) minv = topo(nebs.iup(i), nebs.jdown(j));

    if ((topo(i, j) <= minv) && (i > 0) && (j > 0) && (i < size_x - 1) && (j < size_y - 1)) {
        topo(i, j) = minv + fillincrement;
        fillinpitsandflats(i, j);
        fillinpitsandflats(nebs.iup(i), j);
        fillinpitsandflats(nebs.idown(i), j);
        fillinpitsandflats(i, nebs.jup(j));
        fillinpitsandflats(i, nebs.jdown(j));
        fillinpitsandflats(nebs.iup(i), nebs.jup(j));
        fillinpitsandflats(nebs.idown(i), nebs.jup(j));
        fillinpitsandflats(nebs.idown(i), nebs.jdown(j));
        fillinpitsandflats(nebs.iup(i), nebs.jdown(j));
    }

}

void Flood::run_fillinpitsandflats() {
    for(int i = 0; i < size_x; i++) {
        for(int j = 0; j < size_y; j++) {
            fillinpitsandflats(i, j);
        }
    }
}
