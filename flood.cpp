#include <iostream>
#include "priority_flood.hpp"
#include "parameters.h"
#include "flood.h"
#include "utility.h"

#define fillincrement 0.01

Flood::Flood() : size_x(0), size_y(0), algorithm(2) {}

void Flood::initialise(DEM& topo, Parameters& params) {
    size_x = topo.get_size_x();
    size_y = topo.get_size_y();
    algorithm = params.get_flood_algorithm();

    if (algorithm == 0) {
        std::cout << "<Flood>: using Pelletier's fillinpitsandflats" << std::endl;
        elevation = Array2D<real_type>();
    }
    if (algorithm == 1) {
        std::cout << "<Flood>: using Barnes' original_priority_flood" << std::endl;
        elevation = Array2D<real_type>(size_x, size_y, -9999.0);
    }
    else if (algorithm == 2) {
        std::cout << "<Flood>: using Barnes' priority_flood_epsilon" << std::endl;
        elevation = Array2D<real_type>(size_x, size_y, -9999.0);
    }
    else {
        Util::Error("Unrecognised flood algorithm", 1);
    }
}

void Flood::run(DEM& topo, GridNeighbours& nebs) {
    if (topo.get_size_x() != size_x || topo.get_size_y() != size_y) {
        Util::Error("Must initialse flood object", 1);
    }

    if ((algorithm == 1) || (algorithm == 2)) {
        run_barnes_flood(topo);
    }
    else {
        run_fillinpitsandflats(topo, nebs);
    }
}

void Flood::run_barnes_flood(DEM& topo) {
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

void Flood::fillinpitsandflats(int i, int j, DEM& topo, GridNeighbours& nebs) {
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
        fillinpitsandflats(i, j, topo, nebs);
        fillinpitsandflats(nebs.iup(i), j, topo, nebs);
        fillinpitsandflats(nebs.idown(i), j, topo, nebs);
        fillinpitsandflats(i, nebs.jup(j), topo, nebs);
        fillinpitsandflats(i, nebs.jdown(j), topo, nebs);
        fillinpitsandflats(nebs.iup(i), nebs.jup(j), topo, nebs);
        fillinpitsandflats(nebs.idown(i), nebs.jup(j), topo, nebs);
        fillinpitsandflats(nebs.idown(i), nebs.jdown(j), topo, nebs);
        fillinpitsandflats(nebs.iup(i), nebs.jdown(j), topo, nebs);
    }

}

void Flood::run_fillinpitsandflats(DEM& topo, GridNeighbours& nebs) {
    for(int i = 0; i < size_x; i++) {
        for(int j = 0; j < size_y; j++) {
            fillinpitsandflats(i, j, topo, nebs);
        }
    }
}
