#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include "raster.h"
#include "dem.h"
#include "utility.h"
#include "mfd_flow_router.h"


MFDFlowRouter::MFDFlowRouter() :
        size_x(0), size_y(0) {
}

void MFDFlowRouter::initialise(Raster& flow) {
    size_x = flow.get_size_x();
    size_y = flow.get_size_y();

    // FA boundary values; zero otherwise
    int boundary_layer_pixels = 1;
    fa_bounds = Raster(size_x, size_y, 0.0);
    #pragma omp parallel for
	for (int i = 0; i < size_x; i++)
	{
        for (int m = 0; m < boundary_layer_pixels; m++) {
            fa_bounds(i, m) = flow(i, m);
            fa_bounds(i, size_y - m - 1) = flow(i, size_y - m - 1);
        }
    }
    #pragma omp parallel for
    for (int j = 0; j < size_y; j++)
    {
        for (int m = 0; m < boundary_layer_pixels; m++) {
            fa_bounds(m, j) = flow(m, j);
            fa_bounds(size_x - m - 1, j) = flow(size_x - m - 1, j);
        }
    }
}


void MFDFlowRouter::run(DEM& topo, Raster& flow, GridNeighbours& nebs) {
    // make sure initialise was called before proceeding
    if (topo.get_size_x() != size_x || topo.get_size_y() != size_y) {
        Util::Error("Must initialise flow router", 1);
    }

    // sort data after pit filling
    topo.sort_data();

    // Note that deltax is not used in this computation, so the flow raster represents simply the number of contributing unit cells upstream.
    // Initialise flow to ones everywhere
//    flow.set_data(1.0);
    // initialise flow with pixel area
    flow.set_data(topo.get_deltax() * topo.get_deltax());

    // loop over points starting from highest elevation to lowest
    int t = size_x * size_y;
    while (t > 0)
    {
        t--;
        int i, j;
        topo.get_sorted_ij(t, i, j);
        real_type flow1 = 0;
        real_type flow2 = 0;
        real_type flow3 = 0;
        real_type flow4 = 0;
        real_type flow5 = 0;
        real_type flow6 = 0;
        real_type flow7 = 0;
        real_type flow8 = 0;

        real_type tot = 0;
        if (topo(i, j) > topo(nebs.iup(i), j)) {
            flow1 = pow(topo(i, j) - topo(nebs.iup(i), j), 1.1);
            tot += flow1;
        }
        if (topo(i, j) > topo(nebs.idown(i), j)) {
            flow2 = pow(topo(i, j) - topo(nebs.idown(i), j), 1.1);
            tot += flow2;
        }
        if (topo(i, j) > topo(i, nebs.jup(j))) {
            flow3 = pow(topo(i, j) - topo(i, nebs.jup(j)), 1.1);
            tot += flow3;
        }
        if (topo(i, j) > topo(i, nebs.jdown(j))) {
            flow4 = pow(topo(i, j) - topo(i, nebs.jdown(j)), 1.1);
            tot += flow4;
        }
        if (topo(i, j) > topo(nebs.iup(i), nebs.jup(j))) {
            flow5 = pow((topo(i, j) - topo(nebs.iup(i), nebs.jup(j)))*oneoversqrt2, 1.1);
            tot += flow5;
        }
        if (topo(i, j) > topo(nebs.iup(i), nebs.jdown(j))) {
            flow6 = pow((topo(i, j) - topo(nebs.iup(i), nebs.jdown(j)))*oneoversqrt2, 1.1);
            tot += flow6;
        }
        if (topo(i, j) > topo(nebs.idown(i), nebs.jup(j))) {
            flow7 = pow((topo(i, j) - topo(nebs.idown(i), nebs.jup(j)))*oneoversqrt2, 1.1);
            tot += flow7;
        }
        if (topo(i, j) > topo(nebs.idown(i), nebs.jdown(j))) {
            flow8 = pow((topo(i, j) - topo(nebs.idown(i), nebs.jdown(j)))*oneoversqrt2, 1.1);
            tot += flow8;
        }

        if (tot != 0) {
            real_type reciptot = 1.0 / tot;
            flow1 *= reciptot;
            flow2 *= reciptot;
            flow3 *= reciptot;
            flow4 *= reciptot;
            flow5 *= reciptot;
            flow6 *= reciptot;
            flow7 *= reciptot;
            flow8 *= reciptot;
        }

        flow(nebs.iup(i), j) += flow(i, j) * flow1 + fa_bounds(i, j);     // final fa_bounds(i, j) applies only to edges; zero otherwise
        flow(nebs.idown(i), j) += flow(i, j) * flow2 + fa_bounds(i, j);
        flow(i, nebs.jup(j)) += flow(i, j) * flow3 + fa_bounds(i, j);
        flow(i, nebs.jdown(j)) += flow(i, j) * flow4 + fa_bounds(i, j);
        flow(nebs.iup(i), nebs.jup(j)) += flow(i, j) * flow5 + fa_bounds(i, j);
        flow(nebs.iup(i), nebs.jdown(j)) += flow(i, j) * flow6 + fa_bounds(i, j);
        flow(nebs.idown(i), nebs.jup(j)) += flow(i, j) * flow7 + fa_bounds(i, j);
        flow(nebs.idown(i), nebs.jdown(j)) += flow(i, j) * flow8 + fa_bounds(i, j);
    }
}
