#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "streampower.h"
#include "utility.h"
#include "priority_flood.hpp"
#include "indexx.hpp"
#include "inih/INIReader.h"


/* allocate a float vector with subscript range v[nl..nh] */
std::vector<float> StreamPower::Vector(int nl, int nh)
{
	int size = nh - nl + 1 + NR_END;
	return std::vector<float>(size);
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

std::vector<std::vector<float>> StreamPower::Matrix(int nrl, int nrh, int ncl, int nch)
{
	int rsize = nrh - nrl + 1 + NR_END;
	int csize = nch - ncl + 1 + NR_END;
	return std::vector<std::vector<float>>(rsize, std::vector<float>(csize));
}

float StreamPower::Ran3(std::default_random_engine& generator, std::uniform_real_distribution<float>& distribution)
{
	return distribution(generator);
}

float StreamPower::Gasdev(std::default_random_engine& generator, std::normal_distribution<float>& distribution)
{
	/*
		Assuming this is the same code from here: http://www.stat.berkeley.edu/~paciorek/diss/code/regression.binomial/gasdev.C
		We need to return a standard, normally distributed gaussian random number
	*/

	return distribution(generator);

}

void StreamPower::Tridag(std::vector<float>& a, std::vector<float>& b, std::vector<float>& c, std::vector<float>& r, std::vector<float>& u, int n)
{
	unsigned long j;
	float bet;
	std::vector<float> gam;

	gam = std::vector<float>(n);
	u[0] = r[0] / (bet = b[0]);
	for (j = 1; j <= n-1; j++)
	{
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];
		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}
	for (j = (n - 1); j >= 1; j--)
	{
		u[j - 1] -= gam[j] * u[j];
	}
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

void StreamPower::SetTopo(std::vector<std::vector<float>> t)
{
	topo = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	topo2 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	topoold = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	slope = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	aspect = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	solar_raster = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	shade_raster = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	I_D = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));

	veg = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	veg_old = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	Sed_Track = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	ExposureAge = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	ExposureAge_old = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));

	flow = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	flow1 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	flow2 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	flow3 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	flow4 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	flow5 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	flow6 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	flow7 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	flow8 = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));

    topo_indexx = Indexx<float>(lattice_size_x, lattice_size_y);
    sed_indexx = Indexx<float>(lattice_size_x, lattice_size_y);

	elevation = Array2D<float>(lattice_size_x, lattice_size_y, -9999.0f);

	SetupGridNeighbors();

    std::cout << "DEBUG: size_z, size_y = " << lattice_size_x << ", " << lattice_size_y << std::endl;

	for (int i = 0; i < lattice_size_x; i++)     // Populate model grids
	{
		for (int j = 0; j < lattice_size_y; j++)
		{
			topo[i][j] = t[i][j];
			topoold[i][j] = topo[i][j];
			slope[i][j] = 0;
			aspect[i][j] = 0;
			solar_raster[i][j] = 0;
			veg[i][j] = init_veg;
			flow[i][j] = 1;
			Sed_Track[i][j] = init_sed_track;       // 2 m of overburden to begin			
			ExposureAge[i][j] = init_exposure_age;	// Once over 20, ice is primed for melt.
		}
	}
	InitDiffusion();
}

void StreamPower::SetFA(std::vector<std::vector<float>> f)
{
	FA = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	for (int i = 0; i < lattice_size_x; i++)     // Populate flow accumulation grid; only boundary components are retained, zeros otherwise.
	{
		for (int j = 0; j < lattice_size_y; j++)
		{
			FA[i][j] = f[i][j];
		}
	}
}

void StreamPower::Flood()
{
	// update elev
	for (int i = 0; i < lattice_size_x; i++)
	{
		for (int j = 0; j < lattice_size_y; j++)
		{
			elevation(i, j) = topo[i][j];     //  Change indexing to suit flood subroutine.
		}
	}

	// perform flooding
	original_priority_flood(elevation);

	// update topo
	for (int i = 0; i < lattice_size_x; i++)
	{
		for (int j = 0; j < lattice_size_y; j++)
		{
			topo[i][j] = elevation(i, j);     //  Back to original indexing
		}
	}

}

void StreamPower::MFDFlowRoute(int i, int j)
{
	float tot;

	// Note that deltax is not used in this computation, so the flow raster represents simply the number of contributing unit cells upstream.
	tot = 0;
	if (topo[i][j] > topo[iup[i]][j])
		tot += pow(topo[i][j] - topo[iup[i]][j], 1.1f);
	if (topo[i][j] > topo[idown[i]][j])
		tot += pow(topo[i][j] - topo[idown[i]][j], 1.1f);
	if (topo[i][j] > topo[i][jup[j]])
		tot += pow(topo[i][j] - topo[i][jup[j]], 1.1f);
	if (topo[i][j] > topo[i][jdown[j]])
		tot += pow(topo[i][j] - topo[i][jdown[j]], 1.1f);
	if (topo[i][j] > topo[iup[i]][jup[j]])
		tot += pow((topo[i][j] - topo[iup[i]][jup[j]])*oneoversqrt2, 1.1f);
	if (topo[i][j] > topo[iup[i]][jdown[j]])
		tot += pow((topo[i][j] - topo[iup[i]][jdown[j]])*oneoversqrt2, 1.1f);
	if (topo[i][j] > topo[idown[i]][jup[j]])
		tot += pow((topo[i][j] - topo[idown[i]][jup[j]])*oneoversqrt2, 1.1f);
	if (topo[i][j] > topo[idown[i]][jdown[j]])
		tot += pow((topo[i][j] - topo[idown[i]][jdown[j]])*oneoversqrt2, 1.1f);

	if (topo[i][j] >= topo[iup[i]][j])
		flow1[i][j] = pow(topo[i][j] - topo[iup[i]][j], 1.1f) / tot;
	else flow1[i][j] = 0;
	if (topo[i][j] >= topo[idown[i]][j])
		flow2[i][j] = pow(topo[i][j] - topo[idown[i]][j], 1.1f) / tot;
	else flow2[i][j] = 0;
	if (topo[i][j] >= topo[i][jup[j]])
		flow3[i][j] = pow(topo[i][j] - topo[i][jup[j]], 1.1f) / tot;
	else flow3[i][j] = 0;
	if (topo[i][j] >= topo[i][jdown[j]])
		flow4[i][j] = pow(topo[i][j] - topo[i][jdown[j]], 1.1f) / tot;
	else flow4[i][j] = 0;
	if (topo[i][j] >= topo[iup[i]][jup[j]])
		flow5[i][j] = pow((topo[i][j] - topo[iup[i]][jup[j]])*oneoversqrt2, 1.1f) / tot;
	else flow5[i][j] = 0;
	if (topo[i][j] >= topo[iup[i]][jdown[j]])
		flow6[i][j] = pow((topo[i][j] - topo[iup[i]][jdown[j]])*oneoversqrt2, 1.1f) / tot;
	else flow6[i][j] = 0;
	if (topo[i][j] >= topo[idown[i]][jup[j]])
		flow7[i][j] = pow((topo[i][j] - topo[idown[i]][jup[j]])*oneoversqrt2, 1.1f) / tot;
	else flow7[i][j] = 0;
	if (topo[i][j] >= topo[idown[i]][jdown[j]])
		flow8[i][j] = pow((topo[i][j] - topo[idown[i]][jdown[j]])*oneoversqrt2, 1.1f) / tot;
	else flow8[i][j] = 0;

	flow[iup[i]][j] += flow[i][j] * flow1[i][j] + FA[i][j];    // final FA(j, i) applies only to edges; zero otherwise
	flow[idown[i]][j] += flow[i][j] * flow2[i][j] + FA[i][j];
	flow[i][jup[j]] += flow[i][j] * flow3[i][j] + FA[i][j];
	flow[i][jdown[j]] += flow[i][j] * flow4[i][j] + FA[i][j];
	flow[iup[i]][jup[j]] += flow[i][j] * flow5[i][j] + FA[i][j];
	flow[iup[i]][jdown[j]] += flow[i][j] * flow6[i][j] + FA[i][j];
	flow[idown[i]][jup[j]] += flow[i][j] * flow7[i][j] + FA[i][j];
	flow[idown[i]][jdown[j]] += flow[i][j] * flow8[i][j] + FA[i][j];
}

void StreamPower::InitDiffusion()
{
	//construct diffusional landscape for initial flow routing
	for (int step = 1; step <= 10; step++)
	{
		HillSlopeDiffusion();
		for (int i = 1; i <= lattice_size_x - 2; i++)
		{
			for (int j = 1; j <= lattice_size_y - 2; j++)
			{
				topo[i][j] += 0.1;
				topoold[i][j] += 0.1;
			}
		}
	}
}

void StreamPower::HillSlopeDiffusion()
{
	int i, j, count;
	float term1;

	ax = std::vector<float>(lattice_size_x);
	ay = std::vector<float>(lattice_size_y);
	bx = std::vector<float>(lattice_size_x);
	by = std::vector<float>(lattice_size_y);
	cx = std::vector<float>(lattice_size_x);
	cy = std::vector<float>(lattice_size_y);
	ux = std::vector<float>(lattice_size_x);
	uy = std::vector<float>(lattice_size_y);
	rx = std::vector<float>(lattice_size_x);
	ry = std::vector<float>(lattice_size_y);
	//  D = 10.0;    Previous version used initiliazed value of 10M.
	count = 0;

	while (count < 5)
	{
		count++;
		for (i = 0; i < lattice_size_x; i++)
			for (j = 0; j < lattice_size_y; j++)
				topoold[i][j] = topo[i][j];
		for (i = 0; i < lattice_size_x; i++)
		{
			for (j = 0; j < lattice_size_y; j++)
			{
				term1 = D * ann_timestep / (deltax2);
				if (flow[i][j] < thresholdarea)
				{
					ay[j] = -term1;
					cy[j] = -term1;
					by[j] = 4 * term1 + 1;
					ry[j] = term1 * ( topo[iup[i]][j] + topo[idown[i]][j] ) + topoold[i][j];
				}
				else
				{
					by[j] = 1;
					ay[j] = 0;
					cy[j] = 0;
					ry[j] = topoold[i][j];
				}
				if (j == 0)
				{
					by[j] = 1;
					cy[j] = 0;
					ry[j] = topoold[i][j];
				}
				if (j == lattice_size_y-1)
				{
					by[j] = 1;
					ay[j] = 0;
					ry[j] = topoold[i][j];
				}
			}
			Tridag(ay, by, cy, ry, uy, lattice_size_y);
			for (j = 0; j < lattice_size_y; j++)
				topo[i][j] = uy[j];
		}
		for (i = 0; i < lattice_size_x; i++)
			for (j = 0; j < lattice_size_y; j++)
				topoold[i][j] = topo[i][j];
		for (j = 0; j < lattice_size_y; j++)
		{
			for (i = 0; i < lattice_size_x; i++)
			{
				term1 = D * timestep / ( deltax2 );
				if (flow[i][j] < thresholdarea)
				{
					ax[i] = -term1;
					cx[i] = -term1;
					bx[i] = 4 * term1 + 1;
					rx[i] = term1 * ( topo[i][jup[j]] + topo[i][jdown[j]] ) + topoold[i][j];
				}
				else
				{
					bx[i] = 1;
					ax[i] = 0;
					cx[i] = 0;
					rx[i] = topoold[i][j];
				}
				if (i == 0)
				{
					bx[i] = 1;
					cx[i] = 0;
					rx[i] = topoold[i][j];
				}
				if (i == lattice_size_x-1)
				{
					bx[i] = 1;
					ax[i] = 0;
					rx[i] = topoold[i][j];
				}
			}
			Tridag(ax, bx, cx, rx, ux, lattice_size_x);
			for (i = 0; i < lattice_size_x; i++)
				topo[i][j] = ux[i];
		}
	}
}

void StreamPower::Avalanche(int i, int j)
{

	// NEED TO ASSESS WHETHER PIXEL HAS SEDIMENT, BEFORE FAILURE CALCS?

	if (topo[iup[i]][j] - topo[i][j] > thresh)
		topo[iup[i]][j] = topo[i][j] + thresh;
	if (topo[idown[i]][j] - topo[i][j] > thresh)
		topo[idown[i]][j] = topo[i][j] + thresh;
	if (topo[i][jup[j]] - topo[i][j] > thresh)
		topo[i][jup[j]] = topo[i][j] + thresh;
	if (topo[i][jdown[j]] - topo[i][j] > thresh)
		topo[i][jdown[j]] = topo[i][j] + thresh;
	if (topo[iup[i]][jup[j]] - topo[i][j] > (thresh_diag))
		topo[iup[i]][jup[j]] = topo[i][j] + thresh_diag;
	if (topo[iup[i]][jdown[j]] - topo[i][j] > (thresh_diag))
		topo[iup[i]][jdown[j]] = topo[i][j] + thresh_diag;
	if (topo[idown[i]][jup[j]] - topo[i][j] > (thresh_diag))
		topo[idown[i]][jup[j]] = topo[i][j] + thresh_diag;
	if (topo[idown[i]][jdown[j]] - topo[i][j] > (thresh_diag))
		topo[idown[i]][jdown[j]] = topo[i][j] + thresh_diag;
}

void StreamPower::SlopeAspect(int i, int j) {

	float dzdx = ( ( topo[iup[i]][jdown[j]] + 2 * topo[iup[i]][j] + topo[iup[i]][jup[j]] ) -
		( topo[idown[i]][jdown[j]] + 2 * topo[idown[i]][j] + topo[idown[i]][jup[j]] ) ) /
		8 / deltax;
	float dzdy = ( ( topo[idown[i]][jup[j]] + 2 * topo[i][jup[j]] + topo[iup[i]][jup[j]] ) -
		( topo[idown[i]][jdown[j]] + 2 * topo[i][jdown[j]] + topo[iup[i]][jdown[j]] ) ) /
		8 / deltax;
	aspect[i][j] = atan2(dzdy, dzdx);                            // n.b. Aspect in Radians
	                                                  // N = 0; E = -1/2 pi; S = -pi/+pi; W = +1/2 pi

	slope[i][j] = sqrt(pow(dzdx, 2) + pow(dzdy, 2));              // n.b. Slope in Radians
}

void StreamPower::SunPosition()
{
	int Year1, Day1;
	float m, n, t1, eps, G, C, L, B, alpha, EOT, GHA, AST, LST;

	// Solar declination: angle of suns rays relative to equator. Max pos at summer equinox = +23.5
	r.declination = 23.45 * sin( 360. / 365. * (284. + ct.day) * degrad);

	B = ( 360. / 364. ) * ( ct.day - 81. ) * degrad;              // result in radians
	EOT = 0.165 * sin(2. * B ) - 0.126 * cos( B ) - 0.025 * sin( B ); 
	                           // Equation of Time [hr]
	LST = ct.hour + ( ct.minute/60 ) + ( r.stdmed - r.longitude ) / 15 + EOT - 0;     //  Last term is Daylight Saving (e.g. +1)
	                           // Local Solar Time, correcting for distance from nearest time zone meridian
	r.SHA = 15 * (LST - 12);   // Local Solar Hour

	m = sin(r.lattitude * degrad) * sin(r.declination * degrad);
	n = cos(r.lattitude * degrad) * cos(r.declination * degrad) * cos(r.SHA * degrad);
	r.altitude = asin(m + n) / degrad;

	m = sin(r.lattitude * degrad) * cos(r.declination * degrad) * cos(r.SHA * degrad);
	n = cos(r.lattitude * degrad) * sin(r.declination * degrad);
	r.azimuth = acos ( (m - n) / cos(r.altitude * degrad) ) / degrad - 180;   // in degrees
	if (LST > 12) r.azimuth = -r.azimuth;              // N = 0; E = -90; S = -180/+180; W = +90
}

void StreamPower::SolarInflux(){

// Calculate shading from surrounding terrain
	int i, j;
	float m1, m2;
	float cos_theta, I_DN, I_dH;
	// Solar position, in radians
	float azm = r.azimuth * degrad;            //  Anything in the 'r' object uses degrees; converted here to radians
	float alt = r.altitude * degrad;           //  Invert cos<>sin to obtain zenith angle. 0 is sun at zenith (flat terrain faces up); 90 sun is at the horizon (vertical terrain)
	float lat = r.lattitude * degrad;
	float gamma;

// Solar: Direct and Diffuse
	for (i = 2; i <= lattice_size_x - 1; i++)
	{
		for (j = 2; j <= lattice_size_y - 1; j++)
		{			
			shade_raster[i][j] = ((sin(alt) * cos(slope[i][j])) + (cos(alt) * sin(slope[i][j]) * cos(azm - aspect[i][j])));
			if (shade_raster[i][j] < 0) shade_raster[i][j] = 0;
			else shade_raster[i][j] = 1;                                 //  Temporary, binary on/off for shading effect. This can be improved!

			gamma = azm - aspect[i][j];
			m1 = cos(alt) * cos(gamma) * sin(slope[i][j]);
			m2 = sin(alt) * cos(slope[i][j]);
			cos_theta = m1 + m2;
			I_DN = 1085 * exp(-0.207 / sin(slope[i][j]));        // Direct solar flux striking the horizontal
			I_dH = 0.134 * I_DN;                                 // Diffuse solar flux, based on the ASHRAE Clear Day model
			I_D[i][j] = I_DN * cos_theta * shade_raster[i][j] + I_dH;   // Solar flux moderated by terrain slope and shading effect. Still missing reflected component
		}
	}
}

void StreamPower::MeltExposedIce() {
	// Do not melt or change boundary pixels; maybe 3 deep from edge.
	// Maps of elevation drop in each of 8 directions


}

void StreamPower::Init(std::string parameter_file)
{
    // load parameter file if it exists
    std::cout << "Loading parameter file: " << parameter_file << std::endl;
    INIReader reader(parameter_file);
    if (reader.ParseError() != 0) {
        Util::Error(std::string("Cannot load parameter file: ") + parameter_file, 1);
    }

	// Setup Key Model Variables

	U = reader.GetReal("model", "U", 0.010);        // 'Uplift', m yr^-1
	K = reader.GetReal("model", "K", 0.050);        // Stream Power
	D = reader.GetReal("model", "D", 1.500);        // Diffusion, yr^-1
	melt = reader.GetReal("model", "melt", 250);    // Reciprocal melt rate, for a given radiation input

	deltax = reader.GetReal("model", "deltax", 10.0);   // m; This gets reset after reading ASCII file
	deltax2 = deltax * deltax;                      // m2; Area of a pixel
	nodata = reader.GetReal("model", "nodata", -9999.0);
	xllcorner = reader.GetReal("model", "xllcorner", 0);
	yllcorner = reader.GetReal("model", "yllcorner", 0);

	timestep = reader.GetReal("time", "timestep", 1);   // Time step in hours 
	printinterval = reader.GetInteger("time", "printinterval", 96); // Output timestep, in hours
	ann_timestep = timestep / 8760;    //  Used in formula based on annual rates (e.g. 2 hrs, over 8760 hrs in 365 days)

	thresh = 0.577 * deltax;   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
	thresh_diag = thresh * sqrt2;
	thresholdarea = reader.GetReal("model", "thresholdarea", 0.1);  // Threshold for diffusion domain - to prevent diffusion in channels, etc.

	init_exposure_age = reader.GetReal("model", "init_exposure_age", 0);    // Variables used to initiate exposure time, sed depth and veg age rasters
	init_sed_track = reader.GetReal("model", "init_sed_track", 2);
	init_veg = reader.GetReal("model", "init_veg", 8);

	ct.year = reader.GetInteger("time", "year", 2010);
	ct.day = reader.GetInteger("time", "day", 78);  // 144;              // May 25th is the start of melt/rain season
	ct.hour = reader.GetInteger("time", "hour", 12);    // 24-hr clock
	ct.end_year = reader.GetInteger("time", "end_year", 2015);  // Model execution ends on the first day of this year
	duration = ct.end_year - ct.year;   // Model execution time, in years, keeping in mind melt season is 138 days

	r.lattitude = reader.GetReal("solar_geom", "lattitude", 0); // 67.3;
	r.longitude = reader.GetReal("solar_geom", "longitude", 0); // 134.9;         // Dempster Coordinates
	r.stdmed = reader.GetReal("solar_geom", "stdmed", 0); //9 * 15;         // Standard meridian of nearest time zone. LSTM = (UTC - 9H * 15 deg) n.b. Alaska Time Meridian = 135 deg W
	r.declination = reader.GetReal("solar_geom", "declination", 0.0);
	r.altitude = reader.GetReal("solar_geom", "altitude", 0.0);
	r.azimuth = reader.GetReal("solar_geom", "azimuth", 0.0);

    // input file names
    topo_file = reader.Get("input", "topo", "topo.asc");
    fa_file = reader.Get("input", "FA", "FA.asc");
    sed_file = reader.Get("input", "sed", "SedThickness.asc");
}


void StreamPower::LoadInputs()
{
	SetTopo(ReadArcInfoASCIIGrid(topo_file.c_str()));
	SetFA(ReadArcInfoASCIIGrid(fa_file.c_str()));
	//sp.SetFA(sp.ReadArcInfoASCIIGrid(sed_file.c_str()));   // Option to set sediment thickness
}


void StreamPower::Start()
{
	float deltah, max;
	int idum, i, j, t, step;
	char fname[100];
	sprintf(fname, "erosion_%d.txt", 0);
	PrintState(fname);
	int tstep = 0;    // Counter for printing results to file

	while ( ct.year < ct.end_year )
	{
		// Setup grid index with ranked topo values
        topo_indexx.update_array(topo);
        sed_indexx.update_array(Sed_Track);

		t = 0;	
		// Landsliding, proceeding from high elev to low
		while ( t < lattice_size_x * lattice_size_y )
		{
            topo_indexx.get_ij(t, i, j);
			Avalanche( i , j );
			t++;
		}
		for (j = 0; j < lattice_size_y; j++)
		{
			for (i = 0; i < lattice_size_x; i++)
			{
				topoold[i][j] = topo[i][j];
			}
		}

		// Pit filling
		Flood();

		// Setup grid index again with topo values
		for (j = 0; j < lattice_size_y; j++)
		{
			for (i = 0; i < lattice_size_x; i++)
			{
				// Boundary elements from FA raster
				if (i == 0 || i == (lattice_size_x-1) || j == 0 || j == (lattice_size_y-1) )
					FA[i][j] = FA[i][j];
				else FA[i][j] = 0;   // default value
			}
		}
        topo_indexx.update_array(topo);

		t = lattice_size_x * lattice_size_y;
		while (t > 0)
		{
			t--;
            topo_indexx.get_ij(t, i, j);
			if (i > 0 && i < (lattice_size_x - 1) && j > 0 && j < (lattice_size_y - 1))  // Don't access boundary elements
			{
				MFDFlowRoute(i, j);
			}
		}

		// Diffusive hillslope erosion
		HillSlopeDiffusion();


		// Uplift and Slope/Aspect Calcs
		for (i = 1; i <= lattice_size_x - 2; i++)
		{
			for (j = 1; j <= lattice_size_y - 2; j++)
			{
				topo[i][j] += U * ann_timestep;
				topoold[i][j] += U * ann_timestep;
                SlopeAspect(i, j);
			}
		}

		// Update solar characteristics
		SunPosition();
		SolarInflux();

		// Carry out melt on exposed pixels
		// -- Yet to be coded --

		//Channel erosion
		max = 0;
		for (i = 1; i <= lattice_size_x - 2; i++)
		{
			for (j = 1; j <= lattice_size_y - 2; j++)
			{
				deltah = ann_timestep * K * sqrt( flow[i][j] ) * deltax * slope[i][j];                  // Fluvial erosion law
				topo[i][j] -= deltah;
				if ( topo[i][j] < 0 ) { topo[i][j] = 0; }
				if ( K * sqrt( flow[i][j] ) * deltax > max ) { max = K * sqrt( flow[i][j] ) * deltax; }
			}
		}

		// Update current time
		if ( (ct.hour + timestep) <= 24) { ct.hour = ct.hour + timestep; }
		else {                 // Clicking over 24 hours - day change, possible year change
			if ( ct.day + 1 > 282 ) {                          // Goes to Oct 10th, freezeup and minimum change after that
				 ct.year += 1;
				 ct.day = 144;                                 // Reset to May 25th
			}
			else { ct.day += 1; }
			ct.hour = ( ct.hour + timestep ) - 24;
		}

/*         // Rate adjustment, based on deltax
		if (max > 0.3*deltax / timestep)
		{
			time -= timestep;
			timestep /= 2.0;
			for (i = 2; i <= lattice_size_x - 1; i++)
			{
				for (j = 2; j <= lattice_size_y - 1; j++)
				{
					topo[i][j] = topoold[i][j] - U*timestep;
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
					topoold[i][j] = topo[i][j];
				}

			}

		}
		*/

		std::cout << "Year: " << ct.year << "; Day: " << ct.day << "; Hour: " << ct.hour << std::endl;

		// Write to file at intervals
		tstep += timestep;
		if (tstep >= printinterval) {
			char fname[100];
			sprintf(fname, "erosion_%i_%i_%i.asc", ct.year, ct.day, ct.hour);
			PrintState(fname);
			tstep = 0;
		}
	}

}

void StreamPower::PrintState(char* fname)
{
	std::ofstream file, file2;

	file.open(fname);
	// write arcgrid format
	file << "ncols " << lattice_size_y << std::endl;
	file << "nrows " << lattice_size_x << std::endl;
	file << "xllcorner " << xllcorner << std::endl;
	file << "yllcorner " << yllcorner << std::endl;
	file << "cellsize " << deltax << std::endl;
	file << "NODATA_value " << nodata << std::endl;

	deltax2 = deltax * deltax;
	for (int i = 0; i <= lattice_size_x-1; i++)
	{
		for (int j = 0; j <= lattice_size_y-1; j++)
		{
			file << topo[i][j] << " ";
		}
		file << std::endl;
	}
	file.close();

	/*
	file2.open("params.txt", std::ofstream::app);
	file2 << ct.year << ", " << ct.day << ", " << ct.LocalTime << ", " << ct.UT << ", " << ct.hour 
		<< ", " << r.lattitude << ", " << r.longitude << ", " << r.SHA << ", " << r.stdmed << ", " 
		<< r.declination << ", " << r.altitude << ", " << r.azimuth << std::endl;
	*/

}

std::vector<std::vector<float>> StreamPower::CreateRandomField()
{
	std::vector<std::vector<float>> mat = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));
	std::default_random_engine generator;
	std::normal_distribution<float> distribution(0.0f, 1.0f);
	for (int i = 0; i <= lattice_size_x-1; i++)
	{
		for (int j = 0; j <= lattice_size_y-1; j++)
		{
			mat[i][j] = 0.5 * Gasdev(generator, distribution);
		}
	}
	return mat;
}

std::vector<std::vector<float>> StreamPower::ReadArcInfoASCIIGrid(const char* fname)
{
	std::ifstream in(fname);
	std::vector<std::vector<float>> raster;
	std::string line;

	if (in.fail()) 
	    Util::Warning("Well that didn't work ..!  Missing or invalid file");
	else 
		Util::Warning("Reading raster without any checks or guarantees ...");

	// read 6 lines of metadata
	std::string key;
	in >> key; in >> lattice_size_y; // ncols //NOTE: Pelltier's code was originally written for [x][y] indexing; Saga uses [y][x].
	in >> key; in >> lattice_size_x; // nrows
	in >> key; in >> xllcorner;
	in >> key; in >> yllcorner;
	in >> key; in >> deltax;
	in >> key; in >> nodata;

	raster = std::vector<std::vector<float>>(lattice_size_y, std::vector<float>(lattice_size_x));

	// read data
	for (int x = 0; x < lattice_size_x; x++)
	{
		for (int y = 0; y < lattice_size_y; y++)
		{
			in >> raster[x][y];
		}
	}

	Util::Info("Done reading raster");

	return raster;
}

StreamPower::StreamPower(int nx, int ny) : lattice_size_x(nx), lattice_size_y(ny)
{

}

StreamPower::~StreamPower() {}

std::vector<std::vector<float>> StreamPower::GetTopo()
{
	return topo;
}
