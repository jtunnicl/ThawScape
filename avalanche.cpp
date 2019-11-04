#include <algorithm>
#include "raster.h"
#include "grid_neighbours.h"
#include "utility.h"
#include "avalanche.h"


Avalanche::Avalanche() : size_x(0), size_y(0) {}

void Avalanche::initialise(Raster& topo) {
    size_x = topo.get_size_x();
    size_y = topo.get_size_y();

	deltax = topo.get_deltax();
	deltax2 = deltax * deltax;

	// Code in Init subroutine:
	//	thresh = 0.577 * deltax;   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
	//  thresh_diag = thresh * sqrt2;
	thresh = 0.577 * topo.get_deltax();   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
	thresh_diag = thresh * sqrt2;
}

void Avalanche::run(Raster& topo, Raster& sed_track, Raster& incoming_watts, real_type melt, GridNeighbours& nebs) {
    if (size_x != topo.get_size_x() || size_y != topo.get_size_y()) {
        Util::Error("Must initialise Avalanche object", 1);
    }

	// NEED TO ASSESS WHETHER PIXEL HAS SEDIMENT, BEFORE FAILURE CALCS?

    int t = 0;
    // Landsliding, proceeding from low elev to high
    while (t < size_x * size_y)
    {
        int i, j;
        topo.get_sorted_ij(t, i, j);
        real_type clifftop = 0;

		//---------- Melt Happens Here ----------

		if ((i > 0) && (i < (size_x - 1)) && (j > 0) && (j < (size_y - 1)))  // Do not alter boundary elements
		{
			if (incoming_watts(i, j) != 0) {
				real_type elev_drop = 0;       // Decrease in elevation at central pixel, following ice melt
				real_type accommodation = 0;   // Volume available to fill below central pixel, in the immediate neighbourhood

				// Elevations within 9-element neighbourhood NW-N-NE-W-ctr-E-SW-S-SE
				std::vector<real_type> neighb{ topo(nebs.idown(i), nebs.jup(j)), topo(i, nebs.jup(j)), topo(nebs.iup(i), nebs.jup(j)),
						topo(nebs.idown(i), j), topo(i, j), topo(nebs.iup(i), j),
						topo(nebs.idown(i), nebs.jdown(j)), topo(i, nebs.jdown(j)), topo(nebs.iup(i), nebs.jdown(j)) };
				real_type lowestpixel = *min_element(neighb.begin(), neighb.end());

				// Ice mass lost, based on ablation at each face
				// incoming watts / meltrate / pixel area

				for (int m = 0; m < 8; m++) {
					if (topo(i, j) - neighb[m] > 0) accommodation += deltax2 * (topo(i, j) - neighb[m]);   // sum up all the volume available on pixels below the central pixel
				}

				elev_drop = incoming_watts(i, j) / melt / deltax2;
				if (elev_drop * deltax < accommodation)           // i.e. There is room to accommodate the failed mass in neighbouring cells
					topo(i, j) -= elev_drop;
				else
					topo(i, j) = lowestpixel;


				// Water lost in melt flows downstream (add to 'flow' raster)
			}
		}

		//---------- Avalanching Happens Here ----------

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

		// Update SedTrack, somewhere around here

        t++;
    }
}
