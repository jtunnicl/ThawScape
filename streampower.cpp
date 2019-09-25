#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <string>
#include <map>
#include <chrono>
#include <iomanip>

#include "global_defs.h"
#include "streampower.h"
#include "utility.h"
#include "priority_flood.hpp"
#include "model_time.h"
#include "timer.hpp"
#include "raster.h"
#include "mfd_flow_router.h"
#include "grid_neighbours.h"
#include "parameters.h"
#include "dem.h"
#include "solar_geometry.h"


/* allocate a real_type vector with subscript range v[nl..nh] */
std::vector<real_type> StreamPower::Vector(int nl, int nh)
{
	int size = nh - nl + 1 + NR_END;
	return std::vector<real_type>(size);
}

/* allocate an integer vector with subscript range v[nl..nh] */
std::vector<int> StreamPower::IVector(int nl, int nh)
{
	int size = nh - nl + 1 + NR_END;
	return std::vector<int>(size);
}


std::vector<std::vector<int>> StreamPower::IMatrix(int nrl, int nrh, int ncl, int nch)
{
	int rsize = nrh - nrl + 1 + NR_END;
	int csize = nch - ncl + 1 + NR_END;
	return std::vector<std::vector<int>>(rsize, std::vector<int>(csize));
}

std::vector<std::vector<real_type>> StreamPower::Matrix(int nrl, int nrh, int ncl, int nch)
{
	int rsize = nrh - nrl + 1 + NR_END;
	int csize = nch - ncl + 1 + NR_END;
	return std::vector<std::vector<real_type>>(rsize, std::vector<real_type>(csize));
}

real_type StreamPower::Ran3(std::default_random_engine& generator, std::uniform_real_distribution<real_type>& distribution)
{
	return distribution(generator);
}

real_type StreamPower::Gasdev(std::default_random_engine& generator, std::normal_distribution<real_type>& distribution)
{
	/*
		Assuming this is the same code from here: http://www.stat.berkeley.edu/~paciorek/diss/code/regression.binomial/gasdev.C
		We need to return a standard, normally distributed gaussian random number
	*/

	return distribution(generator);

}

void StreamPower::SetupGridNeighbors()
{
	int i, j;

	idown = std::vector<int>(lattice_size_x);
	iup = std::vector<int>(lattice_size_x);
	jup = std::vector<int>(lattice_size_y);
	jdown = std::vector<int>(lattice_size_y);

	for (i = 0; i <= lattice_size_x - 1; i++)
	{
		idown[i] = i - 1;
		iup[i] = i + 1;
	}
	idown[0] = 0;
	iup[lattice_size_x - 1] = lattice_size_x - 1;
	for (j = 0; j <= lattice_size_y - 1; j++)
	{
		jdown[j] = j - 1;
		jup[j] = j + 1;
	}
	jdown[0] = 0;
	jup[lattice_size_y - 1] = lattice_size_y - 1;
}

void StreamPower::SetTopo()
{
	topo = DEM(params.get_topo_file());
    lattice_size_x = topo.get_size_x();
    lattice_size_y = topo.get_size_y();
    xllcorner = topo.get_xllcorner();
    yllcorner = topo.get_yllcorner();
    deltax = topo.get_deltax();
    deltax2 = deltax * deltax;
    nodata = topo.get_nodata();

	// Landscape Elements
	veg = Raster(lattice_size_x, lattice_size_y, params.get_init_veg());
	veg_old = Raster(lattice_size_x, lattice_size_y);
	Sed_Track = Raster(lattice_size_x, lattice_size_y, params.get_init_sed_track()); // 2m of overburden to begin
	ExposureAge = Raster(lattice_size_x, lattice_size_y, params.get_init_exposure_age());  // Once over 20, ice is primed for melt
	ExposureAge_old = Raster(lattice_size_x, lattice_size_y);

	elevation = Array2D<real_type>(lattice_size_x, lattice_size_y, -9999.0);

	SetupGridNeighbors();
    nebs.setup(lattice_size_x, lattice_size_y);

	InitDiffusion();
}

void StreamPower::SetFA()
{
    flow = Raster(params.get_fa_file());
    lattice_size_x = flow.get_size_x();
    lattice_size_y = flow.get_size_y();
    xllcorner = flow.get_xllcorner();
    yllcorner = flow.get_yllcorner();
    deltax = flow.get_deltax();
    deltax2 = deltax * deltax;
    nodata = flow.get_nodata();
}

void StreamPower::Flood()
{
	// update elev
	for (int i = 0; i < lattice_size_x; i++)
	{
		for (int j = 0; j < lattice_size_y; j++)
		{
			elevation(i, j) = topo(i, j);     //  Change indexing to suit flood subroutine.
		}
	}

	// perform flooding
	original_priority_flood(elevation);

	// update topo
	for (int i = 0; i < lattice_size_x; i++)
	{
		for (int j = 0; j < lattice_size_y; j++)
		{
			topo(i, j) = elevation(i, j);     //  Back to original indexing
		}
	}

}

void StreamPower::InitDiffusion()
{
	//construct diffusional landscape for initial flow routing
	for (int step = 1; step <= 10; step++)
	{
        hillslope_diffusion.run();
		for (int i = 1; i <= lattice_size_x - 2; i++)
		{
			for (int j = 1; j <= lattice_size_y - 2; j++)
			{
				topo(i, j) += 0.1;
			}
		}
	}
}

void StreamPower::Avalanche(int i, int j)
{
	real_type thresh = 0.577 * deltax;   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
	real_type thresh_diag = thresh * sqrt2;

	// Code in Init subroutine:
	//	thresh = 0.577 * deltax;   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
	//  thresh_diag = thresh * sqrt2;
	// NEED TO ASSESS WHETHER PIXEL HAS SEDIMENT, BEFORE FAILURE CALCS?

	real_type clifftop = 0;

	if (topo(iup[i], j) - topo(i, j) > thresh) {
		clifftop = topo(iup[i], j);    // Height of overhanging pixel
		topo(iup[i], j) = std::max((topo(i, j) + thresh), (topo(iup[i], j) - Sed_Track(iup[i], j)));
	}
		//Sed_Track[iup[i]][j]

		//Sed_Track[iup[i]][j] =





	if (topo(idown[i], j) - topo(i, j) > thresh)
		topo(idown[i], j) = topo(i, j) + thresh;
	if (topo(i, jup[j]) - topo(i, j) > thresh)
		topo(i, jup[j]) = topo(i, j) + thresh;
	if (topo(i, jdown[j]) - topo(i, j) > thresh)
		topo(i, jdown[j]) = topo(i, j) + thresh;
	if (topo(iup[i], jup[j]) - topo(i, j) > (thresh_diag))
		topo(iup[i], jup[j]) = topo(i, j) + thresh_diag;
	if (topo(iup[i], jdown[j]) - topo(i, j) > (thresh_diag))
		topo(iup[i], jdown[j]) = topo(i, j) + thresh_diag;
	if (topo(idown[i], jup[j]) - topo(i, j) > (thresh_diag))
		topo(idown[i], jup[j]) = topo(i, j) + thresh_diag;
	if (topo(idown[i], jdown[j]) - topo(i, j) > (thresh_diag))
		topo(idown[i], jdown[j]) = topo(i, j) + thresh_diag;
}



void StreamPower::Init(std::string parameter_file)
{
    // load parameters
    params = Parameters(parameter_file);

    // create model time object
    ct = ModelTime(params);

    // load input data
	SetFA();
	SetTopo();

    // initialise some objects all data loaded
    mfd_flow_router.initialise();
    radiation_model.initialise();
}

void StreamPower::Start()
{
	real_type deltah, max;
	int i, j, t;
    if (params.get_save_topo()) {
        char fname[100];
        sprintf(fname, "erosion_%d.asc", 0);
        topo.save(fname);
    }
    if (params.get_save_flow()) {
        char fname[100];
        sprintf(fname, "flow_%d.asc", 0);
        flow.save(fname);
    }
	int tstep = 0;    // Counter for printing results to file
	std::cout << "U: " << params.get_U() << "; K: " << params.get_K() << "; D: " << params.get_D() << std::endl;

    // set up some timers
    AccumulateTimer<std::chrono::milliseconds> total_time;
    std::vector<std::string> timer_names {"Avalanche", "Flood", "Indexx", "MFDFlowRoute",
        "HillSlopeDiffusion", "Uplift", "SlopeAspect", "SolarCharacteristics", "Melt",
        "ChannelErosion"};
    std::map<std::string, AccumulateTimer<std::chrono::milliseconds> > timers;
    for (auto timer_name : timer_names) {
        timers[timer_name] = AccumulateTimer<std::chrono::milliseconds>();
    }

    total_time.start();
	while ( ct.keep_going() )
	{
		// Setup grid index with ranked topo values
        timers["Indexx"].start();
        topo.sort_data();
//        Sed_Track.sort_data();
        timers["Indexx"].stop();

		t = 0;
		// Landsliding, proceeding from high elev to low
        timers["Avalanche"].start();
		while ( t < lattice_size_x * lattice_size_y )
		{
            topo.get_sorted_ij(t, i, j);
			Avalanche( i , j );
			t++;
		}
        timers["Avalanche"].stop();

		// Pit filling
//        timers["Flood"].start();
//		Flood();
//        timers["Flood"].stop();

        // flow routing
        timers["MFDFlowRoute"].start();
        mfd_flow_router.run();
        timers["MFDFlowRoute"].stop();

		// Diffusive hillslope erosion
        timers["HillSlopeDiffusion"].start();
        hillslope_diffusion.run();
        timers["HillSlopeDiffusion"].stop();

		// Uplift
        timers["Uplift"].start();
        real_type U = params.get_U();
		for (i = 1; i <= lattice_size_x - 2; i++)
		{
			for (j = 1; j <= lattice_size_y - 2; j++)
			{
				topo(i, j) += U * params.get_ann_timestep();
			}
		}
        timers["Uplift"].stop();

        // Slope/Aspect
        timers["SlopeAspect"].start();
        topo.compute_slope_and_aspect(nebs);
        timers["SlopeAspect"].stop();

		// Update solar characteristics
        timers["SolarCharacteristics"].start();
        radiation_model.update_solar_characteristics(ct);
        timers["SolarCharacteristics"].stop();

		// Carry out melt on exposed pixels
        timers["Melt"].start();
        radiation_model.melt_exposed_ice();
        timers["Melt"].stop();

		//Channel erosion
        timers["ChannelErosion"].start();
        real_type K = params.get_K();
		max = 0;
		for (i = 1; i <= lattice_size_x - 2; i++)
		{
			for (j = 1; j <= lattice_size_y - 2; j++)
			{
                deltah = params.get_ann_timestep() * K * sqrt( flow(i, j)/1e6 ) * deltax * topo.slope(i, j);     // Fluvial erosion law;
				topo(i, j) -= deltah;
				//std::cout << "ann_ts: " << ann_timestep << ", K: " << K << ", flow: " << flow(i, j) / 1e6 << ", slope: " << slope(i, j) << std::endl;

				if ( topo(i, j) < 0 ) { topo(i, j) = 0; }
				if ( K * sqrt( flow(i, j)/1e6 ) * deltax > max ) { max = K * sqrt( flow(i, j)/1e6 ) * deltax; }
			}
		}
        timers["ChannelErosion"].stop();

		// Update current time
        ct.increment(params.get_timestep());

/*         // Rate adjustment, based on deltax
		if (max > 0.3*deltax / timestep)
		{
			time -= timestep;
			timestep /= 2.0;
			for (i = 2; i <= lattice_size_x - 1; i++)
			{
				for (j = 2; j <= lattice_size_y - 1; j++)
				{
					topo(i, j) = topoold(i, j) - U*timestep;
				}
			}
		}
		else
		{
			if (max < 0.03*deltax / timestep)
			{
				timestep *= 1.2;
			}
			for (j = 1; j <= lattice_size_y; j++)
			{
				for (i = 1; i <= lattice_size_x; i++)
				{
					topoold(i, j) = topo(i, j);
				}

			}

		}
		*/

        ct.print();

		// Write to file at intervals
		tstep += params.get_timestep();
		if (tstep >= params.get_printinterval()) {
            if (params.get_save_topo()) {
                char fname[100];
                sprintf(fname, "erosion_%i_%i_%i_%.3f.asc", ct.get_year(), ct.get_day(), ct.get_hour(), radiation_model.get_solar_altitude() );
                topo.save(fname);
            }
            if (params.get_save_flow()) {
                char fname[100];
                sprintf(fname, "flow_%i_%i_%i_%.3f.asc", ct.get_year(), ct.get_day(), ct.get_hour(), radiation_model.get_solar_altitude() );
                flow.save(fname);
            }
			tstep = 0;
		}
	}
    total_time.stop();

    // summarise timings
    std::cout << std::endl << "Post run diagnostics..." << std::endl;
    double total_time_secs = total_time.get_total_time() / 1000.0;
    std::cout << "  Total time (excluding input): " << total_time_secs << " s" << std::endl;
    std::cout << "  Component timings:" << std::endl;
    for (auto item : timers) {
        double item_time_secs = item.second.get_total_time() / 1000.0;
        std::cout << "    " << item.first << ": " << item_time_secs << " s";
        std::cout << " (" << item_time_secs / total_time_secs * 100 << " %)" << std::endl;
    }
    std::cout << std::endl;
}

std::vector<std::vector<real_type> > StreamPower::CreateRandomField()
{
	std::vector<std::vector<real_type>> mat = std::vector<std::vector<real_type>>(lattice_size_x, std::vector<real_type>(lattice_size_y));
	std::default_random_engine generator;
    if (params.get_fix_random_seed()) {
        Util::Warning("Fixing random seed - this should only be used for testing/debugging!");
        generator.seed(12345);
    }
	std::normal_distribution<real_type> distribution(0.0f, 1.0);
	for (int i = 0; i <= lattice_size_x-1; i++)
	{
		for (int j = 0; j <= lattice_size_y-1; j++)
		{
			mat[i][j] = 0.5 * Gasdev(generator, distribution);
		}
	}
	return mat;
}

StreamPower::StreamPower(int nx, int ny) : lattice_size_x(nx), lattice_size_y(ny), mfd_flow_router(topo, flow, nebs),
        hillslope_diffusion(topo, flow, nebs, params), radiation_model(topo, Sed_Track, flow, nebs, params)
{

}

StreamPower::~StreamPower() {}
