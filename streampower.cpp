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
#include "model_time.h"
#include "timer.hpp"
#include "raster.h"
#include "mfd_flow_router.h"
#include "grid_neighbours.h"
#include "parameters.h"
#include "dem.h"
#include "solar_geometry.h"
#include "radiation_model.h"
#include "avalanche.h"


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
    flood.initialise(params);
}

void StreamPower::Start()
{
	int i, j;
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
		// Landsliding, proceeding from high elev to low
        if (params.get_avalanche()) {
            timers["Avalanche"].start();
            avalanche.run();
            timers["Avalanche"].stop();
        }

		// Pit filling
        if (params.get_flood()) {
            timers["Flood"].start();
            flood.run();
            timers["Flood"].stop();
        }

        // flow routing
        if (params.get_flow_routing()) {
            timers["MFDFlowRoute"].start();
            mfd_flow_router.run();
            timers["MFDFlowRoute"].stop();
        }

		// Diffusive hillslope erosion
        if (params.get_diffusive_erosion()) {
            timers["HillSlopeDiffusion"].start();
            hillslope_diffusion.run();
            timers["HillSlopeDiffusion"].stop();
        }

		// Uplift
        if (params.get_uplift()) {
            timers["Uplift"].start();
            real_type U = params.get_U();
            #pragma omp parallel for
            for (i = 1; i < lattice_size_x - 1; i++)
            {
                for (j = 1; j < lattice_size_y - 1; j++)
                {
                    topo(i, j) += U * params.get_ann_timestep();
                }
            }
            timers["Uplift"].stop();
        }

        // Slope/Aspect
        if (params.get_melt_component() || params.get_channel_erosion()) {
            timers["SlopeAspect"].start();
            topo.compute_slope_and_aspect(nebs);
            timers["SlopeAspect"].stop();
        }

		// Update solar characteristics
        if (params.get_melt_component()) {
            timers["SolarCharacteristics"].start();
            radiation_model.update_solar_characteristics(ct);
            timers["SolarCharacteristics"].stop();

            // Carry out melt on exposed pixels
            timers["Melt"].start();
            radiation_model.melt_exposed_ice();
            timers["Melt"].stop();
        }

		// Channel erosion
        real_type maxe = 0;
        if (params.get_channel_erosion()) {
            timers["ChannelErosion"].start();
            real_type K = params.get_K();
            #pragma omp parallel for reduction(max: maxe)
            for (i = 1; i <= lattice_size_x - 2; i++)
            {
                for (j = 1; j <= lattice_size_y - 2; j++)
                {
                    real_type flow_sqrt = sqrt(flow(i, j) / 1e6);
                    real_type deltah = params.get_ann_timestep() * K * flow_sqrt * deltax * topo.slope(i, j);     // Fluvial erosion law;
                    topo(i, j) -= deltah;
                    //std::cout << "ann_ts: " << ann_timestep << ", K: " << K << ", flow: " << flow(i, j) / 1e6 << ", slope: " << slope(i, j) << std::endl;

                    if ( topo(i, j) < 0 ) {
                        topo(i, j) = 0;
                    }
                    if ( K * flow_sqrt * deltax > maxe ) {
                        maxe = K * flow_sqrt * deltax;
                    }
                }
            }
            timers["ChannelErosion"].stop();
        }

		// Update current time
        ct.increment(params.get_timestep());

/*         // Rate adjustment, based on deltax
		if (maxe > 0.3*deltax / timestep)
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
			if (maxe < 0.03*deltax / timestep)
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
	std::normal_distribution<real_type> distribution(real_type(0.0), 1.0);
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
        hillslope_diffusion(topo, flow, nebs, params), radiation_model(topo, Sed_Track, flow, nebs, params),
        avalanche(topo, Sed_Track, nebs), flood(topo, nebs)
{

}

StreamPower::~StreamPower() {}
