#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <map>
#include <chrono>

#include "streampower.h"
#include "utility.h"
#include "priority_flood.hpp"
#include "indexx.hpp"
#include "inih/INIReader.h"
#include "time_fcn.h"
#include "timer.hpp"


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
	for (j = 1; j < n; j++)
	{
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];
		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}
	for (j = (n-2); j > 1; j--)
	{
		u[j] -= gam[j + 1] * u[j + 1];
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
	topo = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	topo2 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	topoold = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	slope = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	aspect = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));

	// Radiation Model
	solar_raster = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	shade_raster = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	I_D = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	I_P = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	N_Ip = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	E_Ip = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	S_Ip = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	W_Ip = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	NE_Ip = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	SE_Ip = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	SW_Ip = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	NW_Ip = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));

	// Landscape Elements
	veg = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	veg_old = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	Sed_Track = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	ExposureAge = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	ExposureAge_old = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));

	// Indexing
    topo_indexx = Indexx<float>(lattice_size_x, lattice_size_y);
    sed_indexx = Indexx<float>(lattice_size_x, lattice_size_y);

	elevation = Array2D<float>(lattice_size_x, lattice_size_y, -9999.0f);

	SetupGridNeighbors();

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
			Sed_Track[i][j] = init_sed_track;       // 2 m of overburden to begin			
			ExposureAge[i][j] = init_exposure_age;	// Once over 20, ice is primed for melt.
		}
	}
	InitDiffusion();
}

void StreamPower::SetFA(std::vector<std::vector<float>> f)
{

	flow = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	flow1 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	flow2 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	flow3 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	flow4 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	flow5 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	flow6 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	flow7 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	flow8 = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	FA_Bounds = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));

	for (int i = 0; i < lattice_size_x; i++)
	{
		for (int j = 0; j < lattice_size_y; j++)
		{
			flow[i][j] = f[i][j];
			flow1[i][j] = 1;
			flow2[i][j] = 1;
			flow3[i][j] = 1;
			flow4[i][j] = 1;
			flow5[i][j] = 1;
			flow6[i][j] = 1;
			flow7[i][j] = 1;
			flow8[i][j] = 1;
			if (i == 0 || j == 0 || i == (lattice_size_x - 1) || j == (lattice_size_y - 1))
				FA_Bounds[i][j] = f[i][j];   // FA boundary values; zero otherwise
			else
				FA_Bounds[i][j] = 0;
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

	if (topo[i][j] > topo[iup[i]][j])
		flow1[i][j] = pow(topo[i][j] - topo[iup[i]][j], 1.1f) / tot;
	else flow1[i][j] = 0;
	if (topo[i][j] > topo[idown[i]][j])
		flow2[i][j] = pow(topo[i][j] - topo[idown[i]][j], 1.1f) / tot;
	else flow2[i][j] = 0;
	if (topo[i][j] > topo[i][jup[j]])
		flow3[i][j] = pow(topo[i][j] - topo[i][jup[j]], 1.1f) / tot;
	else flow3[i][j] = 0;
	if (topo[i][j] > topo[i][jdown[j]])
		flow4[i][j] = pow(topo[i][j] - topo[i][jdown[j]], 1.1f) / tot;
	else flow4[i][j] = 0;
	if (topo[i][j] > topo[iup[i]][jup[j]])
		flow5[i][j] = pow((topo[i][j] - topo[iup[i]][jup[j]])*oneoversqrt2, 1.1f) / tot;
	else flow5[i][j] = 0;
	if (topo[i][j] > topo[iup[i]][jdown[j]])
		flow6[i][j] = pow((topo[i][j] - topo[iup[i]][jdown[j]])*oneoversqrt2, 1.1f) / tot;
	else flow6[i][j] = 0;
	if (topo[i][j] > topo[idown[i]][jup[j]])
		flow7[i][j] = pow((topo[i][j] - topo[idown[i]][jup[j]])*oneoversqrt2, 1.1f) / tot;
	else flow7[i][j] = 0;
	if (topo[i][j] > topo[idown[i]][jdown[j]])
		flow8[i][j] = pow((topo[i][j] - topo[idown[i]][jdown[j]])*oneoversqrt2, 1.1f) / tot;
	else flow8[i][j] = 0;

	flow[iup[i]][j] += flow[i][j] * flow1[i][j] + FA_Bounds[i][j];     // final FA_Bounds[i][j] applies only to edges; zero otherwise
	flow[idown[i]][j] += flow[i][j] * flow2[i][j] + FA_Bounds[i][j];
	flow[i][jup[j]] += flow[i][j] * flow3[i][j] + FA_Bounds[i][j];
	flow[i][jdown[j]] += flow[i][j] * flow4[i][j] + FA_Bounds[i][j];
	flow[iup[i]][jup[j]] += flow[i][j] * flow5[i][j] + FA_Bounds[i][j];
	flow[iup[i]][jdown[j]] += flow[i][j] * flow6[i][j] + FA_Bounds[i][j];
	flow[idown[i]][jup[j]] += flow[i][j] * flow7[i][j] + FA_Bounds[i][j];
	flow[idown[i]][jdown[j]] += flow[i][j] * flow8[i][j] + FA_Bounds[i][j];
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
	aspect[i][j] = atan2(dzdy, dzdx);                             // n.b. Aspect in Radians
	slope[i][j] = sqrt(pow(dzdx, 2) + pow(dzdy, 2));              // n.b. Slope in Radians
}

void StreamPower::SunPosition()
{

	/*  Documentation for Sun Position and Solar influx come from the following sources:
	      ME 4131 THERMAL ENVIRONMENTAL ENGINEERING LABORATORY MANUAL, Appendix D
	      http://www.me.umn.edu/courses/me4131/LabManual/AppDSolarRadiation.pdf

		  Kumar, L., Skidmore, A. and Knowles, E. 1997. Modelling topographic variation 
		  in solar radiation in a GIS environment. Int J. Information Science, 11(5), p.475-497

		  Angus, R., Muneer, T. 1993. Sun position for daylight models: Precise algorithms for
		  determination. Lighting Research and Technology 25(2) 81-83.
	*/

	int Year1, Day1;
	float m, n, t1, eps, G, C, L, B, alpha, EOT, GHA, AST, LST;
    int day = ct.get_day();
    int hour = ct.get_hour();
    int minute = ct.get_minute();

	// Solar declination: angle of suns rays relative to equator. Max pos at summer equinox = +23.5
	r.declination = 23.45 * sin( 360. / 365. * (284. + day) * degrad);

	B = ( 360. / 364. ) * ( day - 81. ) * degrad;              // result in radians
	EOT = 0.165 * sin(2. * B ) - 0.126 * cos( B ) - 0.025 * sin( B ); 
	                           // Equation of Time [hr]
	LST = hour + ( minute/60 ) + ( r.stdmed - r.longitude ) / 15 + EOT - 0;     //  Last term is Daylight Saving (e.g. +1)
	                           // Local Solar Time, correcting for distance from nearest time zone meridian
	r.SHA = 15 * (LST - 12);   // Local Solar Hour  (negative before solar noon, positive after)

	m = sin(r.lattitude * degrad) * sin(r.declination * degrad);
	n = cos(r.lattitude * degrad) * cos(r.declination * degrad) * cos(r.SHA * degrad);
	r.altitude = asin(m + n) / degrad;

	m = sin(r.lattitude * degrad) * cos(r.declination * degrad) * cos(r.SHA * degrad);
	n = cos(r.lattitude * degrad) * sin(r.declination * degrad);
	r.azimuth = acos ( (m - n) / cos(r.altitude * degrad) ) / degrad - 180;   // in degrees

	if (r.SHA > 0) r.azimuth = -r.azimuth;              // N = 0; E = -90; S = -180/+180; W = +90

}

void StreamPower::SolarInflux(){

// Calculate shading from surrounding terrain
	int i, j, m;
	float m1, m2, m3, m4, m5, M, asp360, d80;
	float cos_theta, cos_i, I_o, I_r, tau_b, I_DN, I_dH;
	std::vector<float> cos_i80, asp_4;

	// Solar position, in radians
	float azm = r.azimuth * degrad;            //  Anything in the 'r' object uses degrees; converted here to radians
	float alt = r.altitude * degrad;           //  Invert cos<>sin to obtain zenith angle. 0 is sun at zenith (flat terrain faces up); 90 sun is at the horizon (vertical terrain)
	float lat = r.lattitude * degrad;
	float dec = r.declination * degrad;
	float sha = r.SHA * degrad;
	float gamma;

	cos_i80 = std::vector<float>(8);
	d80 = (80 * degrad);
	asp_4 = { 0, 90, 180, 270, 45, 135, 225, 315 };

	I_o = 1367 * (1 + 0.0344 * cos(360 * ct.get_day() / 365 * degrad));
	M = sqrt(1229. + pow((614. * sin(alt)), 2.)) - 614 * sin(alt);               // Air mass ratio  (Keith and Kreider 1978)
	tau_b = 0.56 * (exp(-0.65 * M) + exp(-0.095 * M));                               // Atmospheric transmittance for beam radiation

// Solar: Direct and Diffuse
	for (i = 2; i <= lattice_size_x - 1; i++)
	{
		for (j = 2; j <= lattice_size_y - 1; j++)
		{			
			shade_raster[i][j] = ((sin(alt) * cos(slope[i][j])) + (cos(alt) * sin(slope[i][j]) * cos(azm - aspect[i][j])));
			if (shade_raster[i][j] < 0) shade_raster[i][j] = 0;

			asp360 = aspect[i][j];                 // Change aspect coordinates for flux estimates:
			if ( asp360 < PI) asp360 += PI / 2;    // N = 0; E = 1/2 pi; S = pi; W = 1.5 pi
			if (asp360 < 0) asp360 += 2 * PI;

			// Solar radiation striking a tilted surface
			m1 = sin(lat) * cos(slope[i][j]);
			m2 = cos(lat) * sin(slope[i][j]) * cos( asp360 );
			m3 = cos(lat) * cos(slope[i][j]);
			m4 = sin(lat) * cos(slope[i][j]) * cos( asp360 );
			m5 = cos(dec) * sin(slope[i][j]) * sin( sha );
			// Incident angle of incoming beam radiation

			cos_i = sin(dec) * (m1 - m2) + cos(dec) * cos(sha) * (m3 + m4) + m5;

			I_P[i][j] = (I_o * tau_b) * cos_i;
			if (I_P[i][j] < 0) I_P[i][j] = 0;
			if (shade_raster[i][j] < 0) I_P[i][j] = 0;
			I_D[i][j] = I_o * ( 0.271 - 0.294 * tau_b ) * pow ( cos(slope[i][j] / 2 ),  2 ) * sin(alt);  // Diffuse insolation
			I_r = 0.2 * I_o * ( 0.271 + 0.706 * tau_b ) * pow ( sin(slope[i][j] / 2 ),  2 ) * sin(alt);  // Reflected insolation
			I_P[i][j] += I_D[i][j] + I_r;

			//  Maps of incident radiation on near-vertical (>80 deg) slopes in each cardinal direction
			//  8 Ip values and 'cos_i's; one for each direction

			m1 = sin(lat) * cos(d80);
			m3 = cos(lat) * cos(d80);
			m5 = cos(dec) * sin(d80) * sin(sha);

			for (m = 0; m < 8; m++)
			{
				m2 = cos(lat) * sin(d80) * cos(asp_4[m] * degrad);
				m4 = sin(lat) * cos(d80) * cos(asp_4[m] * degrad);
				cos_i80[m] = sin(dec) * (m1 - m2) + cos(dec) * cos(sha) * (m3 + m4) + m5;
			}

			N_Ip[i][j] = (I_o * tau_b) * cos_i80[0] * shade_raster[i][j];
			E_Ip[i][j] = (I_o * tau_b) * cos_i80[1] * shade_raster[i][j];
			S_Ip[i][j] = (I_o * tau_b) * cos_i80[2] * shade_raster[i][j];
			W_Ip[i][j] = (I_o * tau_b) * cos_i80[3] * shade_raster[i][j];
			NE_Ip[i][j] = (I_o * tau_b) * cos_i80[4] * shade_raster[i][j];
			SE_Ip[i][j] = (I_o * tau_b) * cos_i80[5] * shade_raster[i][j];
			SW_Ip[i][j] = (I_o * tau_b) * cos_i80[6] * shade_raster[i][j];
			NW_Ip[i][j] = (I_o * tau_b) * cos_i80[7] * shade_raster[i][j];
		}
	}
}

void StreamPower::MeltExposedIce(int i, int j) {

	float IceTop;

	// Do not melt or change boundary pixels; maybe 3 deep from edge.
	// Maps of elevation drop in each of 8 directions

	// Deal with melt-out first

	IceTop = topo[i][j] - Sed_Track[i][j];      // Sediment-Ice interface at i,j

	// Planar faces
	if (topo[iup[i]][j] < IceTop) {
		
	
	}

	// Radiation Input and resultant melt

	// I_D[i][j]


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
	K = reader.GetReal("model", "K", 0.001);        // Stream Power, yr^-1
	D = reader.GetReal("model", "D", 1.500);        // Diffusion, yr^-1
	melt = reader.GetReal("model", "melt", 250);    // Reciprocal melt rate, for a given radiation input

	deltax = reader.GetReal("model", "deltax", 10.0);   // m; This gets reset after reading ASCII file
	deltax2 = deltax * deltax;                      // m2; Area of a pixel
	nodata = reader.GetReal("model", "nodata", -9999.0);
	xllcorner = reader.GetReal("model", "xllcorner", 0);
	yllcorner = reader.GetReal("model", "yllcorner", 0);

	timestep = reader.GetReal("time", "timestep", 1);   // Time step in hours 
	printinterval = reader.GetInteger("time", "printinterval", 1); // Output timestep, in hours
	ann_timestep = timestep / 8760;    //  Used in formula based on annual rates (e.g. 2 hrs, over 8760 hrs in 365 days)

	thresh = 0.577 * deltax;   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
	thresh_diag = thresh * sqrt2;
	thresholdarea = reader.GetReal("model", "thresholdarea", 1e35);  // Threshold for diffusion domain - to prevent diffusion in channels, etc. (m2)

	init_exposure_age = reader.GetReal("model", "init_exposure_age", 0);    // Variables used to initiate exposure time, sed depth and veg age rasters
	init_sed_track = reader.GetReal("model", "init_sed_track", 2);
	init_veg = reader.GetReal("model", "init_veg", 8);

	int year = reader.GetInteger("time", "year", 2010);
	int day = reader.GetInteger("time", "day", 145);  // 144;              // May 25th is the start of melt/rain season
	int hour = reader.GetInteger("time", "hour", 12);    // 24-hr clock
	int minute = reader.GetInteger("time", "hour", 0);   // 0 in most cases
	int end_year = reader.GetInteger("time", "end_year", 2015);  // Model execution ends on the first day of this year
    int end_day = reader.GetInteger("time", "end_day", 1);
    ct = time_fcn(year, day, hour, minute, end_year, end_day);
	duration = end_year - year;   // Model execution time, in years, keeping in mind melt season is 138 days

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

    // should we fix the random number seed (for testing)
    fix_random_seed = reader.GetBoolean("random", "fix_seed", false);
}

void StreamPower::LoadInputs()
{
	SetFA(ReadArcInfoASCIIGrid(fa_file.c_str()));
	SetTopo(ReadArcInfoASCIIGrid(topo_file.c_str()));
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
	std::cout << "U: " << U << "; K: " << K << "; D: " << D << std::endl;

    // set up some timers
    AccumulateTimer<std::chrono::milliseconds> total_time;
    std::vector<std::string> timer_names {"Avalanche", "Flood", "Indexx", "MFDFlowRoute",
        "HillSlopeDiffusion", "UpliftAndSlopeAspect", "SolarCharacteristics", "Melt",
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
        topo_indexx.update_array(topo);
        sed_indexx.update_array(Sed_Track);
        timers["Indexx"].stop();

		t = 0;	
		// Landsliding, proceeding from high elev to low
        timers["Avalanche"].start();
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
        timers["Avalanche"].stop();

		// Pit filling
        timers["Flood"].start();
		Flood();
        timers["Flood"].stop();

		// Setup grid index again with topo values
        timers["Indexx"].start();
        topo_indexx.update_array(topo);
        timers["Indexx"].stop();


        timers["MFDFlowRoute"].start();
		t = lattice_size_x * lattice_size_y;
		while (t > 0)
		{
			t--;
            topo_indexx.get_ij(t, i, j);
            if ( ( i > 3 ) && ( i < (lattice_size_x - 3) ) && ( j > 3 ) && ( j < (lattice_size_y - 3) ) )  MFDFlowRoute(i, j);// Do not alter boundary elements
		}
        timers["MFDFlowRoute"].stop();


		// Diffusive hillslope erosion
        timers["HillSlopeDiffusion"].start();
		HillSlopeDiffusion();
        timers["HillSlopeDiffusion"].stop();

		// Uplift and Slope/Aspect Calcs
        timers["UpliftAndSlopeAspect"].start();
		for (i = 1; i <= lattice_size_x - 2; i++)
		{
			for (j = 1; j <= lattice_size_y - 2; j++)
			{
				topo[i][j] += U * ann_timestep;
				topoold[i][j] += U * ann_timestep;
                SlopeAspect(i, j);
			}
		}
        timers["UpliftAndSlopeAspect"].stop();

		// Update solar characteristics
        timers["SolarCharacteristics"].start();
		SunPosition();
		SolarInflux();
        timers["SolarCharacteristics"].stop();

		// Carry out melt on exposed pixels
        timers["Melt"].start();
		t = lattice_size_x * lattice_size_y;
		while (t > 0)
		{
			t--;
			topo_indexx.get_ij(t, i, j);
			if ((i > 0) && (i < (lattice_size_x - 1)) && (j > 0) && (j < (lattice_size_y - 1)))  // Do not alter boundary elements
			{
				MeltExposedIce(i, j);
			}
		}
        timers["Melt"].stop();

		//Channel erosion
        timers["ChannelErosion"].start();
		max = 0;
		for (i = 1; i <= lattice_size_x - 2; i++)
		{
			for (j = 1; j <= lattice_size_y - 2; j++)
			{				deltah = ann_timestep * K * sqrt( flow[i][j]/1e6 ) * deltax * slope[i][j];     // Fluvial erosion law; 
				topo[i][j] -= deltah;
				//std::cout << "ann_ts: " << ann_timestep << ", K: " << K << ", flow: " << flow[i][j] / 1e6 << ", slope: " << slope[i][j] << std::endl;

				if ( topo[i][j] < 0 ) { topo[i][j] = 0; }
				if ( K * sqrt( flow[i][j]/1e6 ) * deltax > max ) { max = K * sqrt( flow[i][j]/1e6 ) * deltax; }
			}
		}
        timers["ChannelErosion"].stop();

		// Update current time
        ct.increment(timestep);

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

        ct.print();

		// Write to file at intervals
		tstep += timestep;
		if (tstep >= printinterval) {
			char fname[100];
			sprintf(fname, "erosion_%i_%i_%.3f.asc", ct.get_day(), ct.get_hour(), r.altitude );
			PrintState(fname);
			tstep = 0;
		}
	}
    total_time.stop();

    // summarise timings (TODO: make this optional?)
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
	std::vector<std::vector<float>> mat = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));
	std::default_random_engine generator;
    if (fix_random_seed) {
        Util::Warning("Fixing random seed - this should only be used for testing/debugging!");
        generator.seed(12345);
    }
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
	    Util::Error("Well that didn't work ..!  Missing or invalid file: " + std::string(fname), 1);
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

	raster = std::vector<std::vector<float>>(lattice_size_x, std::vector<float>(lattice_size_y));

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
