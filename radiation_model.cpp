#include <cmath>
#include <algorithm>
#include "raster.h"
#include "grid_neighbours.h"
#include "global_defs.h"
#include "solar_geometry.h"
#include "model_time.h"
#include "utility.h"
#include "radiation_model.h"


RadiationModel::RadiationModel() : lattice_size_x(0), lattice_size_y(0) {}


/// Allocates all Rasters to the correct size. This should be called before running the radiation model.
void RadiationModel::initialise(Raster& topo, Parameters& params) {
    lattice_size_x = topo.get_size_x();
    lattice_size_y = topo.get_size_y();
    deltax = topo.get_deltax();
    deltax2 = deltax * deltax;

	solar_raster = Raster(lattice_size_x, lattice_size_y, 0.0);
	shade_raster = Raster(lattice_size_x, lattice_size_y);
    incoming_watts = Raster(lattice_size_x, lattice_size_y, 0.0);
	I_D = Raster(lattice_size_x, lattice_size_y);
	I_R = Raster(lattice_size_x, lattice_size_y);
	I_P = Raster(lattice_size_x, lattice_size_y);
	N_Ip = Raster(lattice_size_x, lattice_size_y);
	E_Ip = Raster(lattice_size_x, lattice_size_y);
	S_Ip = Raster(lattice_size_x, lattice_size_y);
	W_Ip = Raster(lattice_size_x, lattice_size_y);
	NE_Ip = Raster(lattice_size_x, lattice_size_y);
	SE_Ip = Raster(lattice_size_x, lattice_size_y);
	SW_Ip = Raster(lattice_size_x, lattice_size_y);
	NW_Ip = Raster(lattice_size_x, lattice_size_y);

    r = SolarGeometry(params);
}

/// This routine first updates the sun position and then computes the solar influx.
void RadiationModel::update_solar_characteristics(Raster& topo, ModelTime& ct) {
    if (lattice_size_x != topo.get_size_x() || lattice_size_y != topo.get_size_y()) {
        Util::Error("Must initialise RadiationModel", 1);
    }

    r.sun_position(ct);
    solar_influx(topo, ct);
}

void RadiationModel::solar_influx(Raster& topo, ModelTime& ct) {
    // Calculate shading from surrounding terrain
    real_type d80, I_o, M;
    real_type tau_b;
	std::vector<real_type> asp_4;

	// Solar position, in radians
	real_type azm = r.get_azimuth() * degrad;            //  Anything in the 'r' object uses degrees; converted here to radians
	real_type alt = r.get_altitude() * degrad;           //  Invert cos<>sin to obtain zenith angle. 0 is sun at zenith (flat terrain faces up); 90 sun is at the horizon (vertical terrain)
	real_type lat = r.get_latitude() * degrad;
	real_type dec = r.get_declination() * degrad;
	real_type sha = r.get_SHA() * degrad;

	d80 = (80 * degrad);
	asp_4 = { 0, 90, 180, 270, 45, 135, 225, 315 };

	I_o = 1367 * (1 + 0.0344 * cos(360 * ct.get_day() / 365 * degrad));
	M = sqrt(1229. + pow((614. * sin(alt)), 2.)) - 614 * sin(alt);               // Air mass ratio  (Keith and Kreider 1978)
	tau_b = 0.56 * (exp(-0.65 * M) + exp(-0.095 * M));                               // Atmospheric transmittance for beam radiation

    // Solar: Direct and Diffuse
    #pragma omp parallel for
	for (int i = 2; i <= lattice_size_x - 1; i++)
	{
		for (int j = 2; j <= lattice_size_y - 1; j++)
		{
            std::vector<real_type> cos_i80 = std::vector<real_type>(8);
            real_type m1, m2, m3, m4, m5, asp360;
            real_type cos_i;

			shade_raster(i, j) = ((sin(alt) * cos(topo.slope(i, j))) + (cos(alt) * sin(topo.slope(i, j)) * cos(azm - topo.aspect(i, j))));
			if (shade_raster(i, j) < 0) shade_raster(i, j) = 0;

			asp360 = topo.aspect(i, j);                 // Change aspect coordinates for flux estimates:
			if ( asp360 < PI) asp360 += PI / 2;    // N = 0; E = 1/2 pi; S = pi; W = 1.5 pi
			if (asp360 < 0) asp360 += 2 * PI;

			// Solar radiation striking a tilted surface
			m1 = sin(lat) * cos(topo.slope(i, j));
			m2 = cos(lat) * sin(topo.slope(i, j)) * cos( asp360 );
			m3 = cos(lat) * cos(topo.slope(i, j));
			m4 = sin(lat) * cos(topo.slope(i, j)) * cos( asp360 );
			m5 = cos(dec) * sin(topo.slope(i, j)) * sin( sha );
			// Incident angle of incoming beam radiation

			cos_i = sin(dec) * (m1 - m2) + cos(dec) * cos(sha) * (m3 + m4) + m5;

			I_P(i, j) = (I_o * tau_b) * cos_i;
			if (I_P(i, j) < 0) I_P(i, j) = 0;
			if (shade_raster(i, j) < 0) I_P(i, j) = 0;
			I_D(i, j) = I_o * ( 0.271 - 0.294 * tau_b ) * pow ( cos(topo.slope(i, j) / 2 ),  2 ) * sin(alt);  // Diffuse insolation
			I_R(i, j) = 0.2 * I_o * ( 0.271 + 0.706 * tau_b ) * pow ( sin(topo.slope(i, j) / 2 ),  2 ) * sin(alt);  // Reflected insolation
			if (I_R(i, j) < 0) I_R(i, j) = 0;
			I_P(i, j) += I_D(i, j);

			//  Maps of incident radiation on near-vertical (>80 deg) slopes in each cardinal direction
			//  8 Ip values and 'cos_i's; one for each direction

			m1 = sin(lat) * cos(d80);
			m3 = cos(lat) * cos(d80);
			m5 = cos(dec) * sin(d80) * sin(sha);

			for (int m = 0; m < 8; m++)
			{
				m2 = cos(lat) * sin(d80) * cos(asp_4[m] * degrad);
				m4 = sin(lat) * cos(d80) * cos(asp_4[m] * degrad);
				cos_i80[m] = sin(dec) * (m1 - m2) + cos(dec) * cos(sha) * (m3 + m4) + m5;
			}

			N_Ip(i, j) = (I_o * tau_b) * cos_i80[0] * shade_raster(i, j);
			E_Ip(i, j) = (I_o * tau_b) * cos_i80[1] * shade_raster(i, j);
			S_Ip(i, j) = (I_o * tau_b) * cos_i80[2] * shade_raster(i, j);
			W_Ip(i, j) = (I_o * tau_b) * cos_i80[3] * shade_raster(i, j);
			NE_Ip(i, j) = (I_o * tau_b) * cos_i80[4] * shade_raster(i, j);
			SE_Ip(i, j) = (I_o * tau_b) * cos_i80[5] * shade_raster(i, j);
			SW_Ip(i, j) = (I_o * tau_b) * cos_i80[6] * shade_raster(i, j);
			NW_Ip(i, j) = (I_o * tau_b) * cos_i80[7] * shade_raster(i, j);
		}
	}
}

/// This is a simple loop over all points in the DEM, i.e. there is no ordering from low to high
/// elevations.
void RadiationModel::melt_potential(Raster& topo, Raster& Sed_Track, Raster& flow, GridNeighbours& nebs) {
    if (lattice_size_x != topo.get_size_x() || lattice_size_y != topo.get_size_y()) {
        Util::Error("Must initialise RadiationModel", 1);
    }

    // first compute incoming watts at all pixels (except boundary?)
    incoming_watts.set_data(0.0);
    #pragma omp parallel for
    for (int i = 1; i < lattice_size_x - 1; i++) {
        for (int j = 1; j < lattice_size_y - 1; j++) {
            real_type N, E, S, W, NE, SE, SW, NW;
            real_type incoming = 0;
            real_type lowestpixel;         // Elevation of the lowest pixel in the 9-element neighbourhood.

            std::vector<real_type> neighb{ topo(nebs.idown(i), nebs.jup(j)), topo(i, nebs.jup(j)), topo(nebs.iup(i), nebs.jup(j)),    // Elevations within 9-element neighbourhood NW-N-NE-W-ctr-E-SW-S-SE
                topo(nebs.idown(i), j), topo(i, j), topo(nebs.iup(i), j),
                topo(nebs.idown(i), nebs.jdown(j)), topo(i, nebs.jdown(j)), topo(nebs.iup(i), nebs.jdown(j)) };

            // get the value of the lowest pixel within the 9-element neighbourhood
            lowestpixel = *min_element(neighb.begin(), neighb.end());

            if (topo(i, j) > lowestpixel)      // If any neighbouring pixels are higher than central pixel, then proceed with melt/avalanche algorithm
            {
                // Extent (m2) of exposed faces in each of 8 directions
                N = std::max<real_type>((topo(i, j) - Sed_Track(i, j) - topo(i, nebs.jup(j))), 0.0) * deltax * 0.8;  // If ice is exposed, positive value, otherwise zero
                E = std::max<real_type>((topo(i, j) - Sed_Track(i, j) - topo(nebs.iup(i), j)), 0.0) * deltax * 0.8;
                S = std::max<real_type>((topo(i, j) - Sed_Track(i, j) - topo(i, nebs.jdown(j))), 0.0) * deltax * 0.8;
                W = std::max<real_type>((topo(i, j) - Sed_Track(i, j) - topo(nebs.idown(i), j)), 0.0) * deltax * 0.8;
                NE = std::max<real_type>((topo(i, j) - Sed_Track(i, j) - topo(nebs.iup(i), nebs.jup(j))), 0.0) * deltax * 0.2;  //  Faces have 0.8 of deltax resolution; corners have 0.2
                SE = std::max<real_type>((topo(i, j) - Sed_Track(i, j) - topo(nebs.iup(i), nebs.jdown(j))), 0.0) * deltax * 0.2;
                SW = std::max<real_type>((topo(i, j) - Sed_Track(i, j) - topo(nebs.idown(i), nebs.jdown(j))), 0.0) * deltax * 0.2;
                NW = std::max<real_type>((topo(i, j) - Sed_Track(i, j) - topo(nebs.idown(i), nebs.jup(j))), 0.0) * deltax * 0.2;

                // Radiative flux (m2 * W·m-2 = W) to ice for each face and corner of the pixel block

                if (N > 0) {
                    incoming = N * I_P(i, j) * N_Ip(i, j);         //  Area exposed (m2) * direct+diffuse (W·m-2) * vertical faces (W·m-2)
                    if ((topo.aspect(i, nebs.jup(j)) < -7 * (PI / 8)) || (topo.aspect(i, nebs.jup(j)) > 7 * (PI / 8)))
                        incoming += N * I_R(i, nebs.jup(j));
                }           //  Add reflected radiation component (I_R), if applicable (e.g. pixel to the North is sloping Southward)
                if (E > 0) {
                    incoming += E * I_P(i, j) * E_Ip(i, j);
                    if ((topo.aspect(nebs.iup(i), j) < 7 * (PI / 8)) && (topo.aspect(nebs.iup(i), j) > 3 * (PI / 8)))
                        incoming += E * I_R(nebs.iup(i), j);
                }
                if (S > 0) {
                    incoming += S * I_P(i, j) * S_Ip(i, j);
                    if ((topo.aspect(i, nebs.jdown(j)) < (PI / 8)) && (topo.aspect(i, nebs.jdown(j)) > -1 * (PI / 8)))
                        incoming += S * I_R(i, nebs.jdown(j));
                }
                if (W > 0) {
                    incoming += W * I_P(i, j) * E_Ip(i, j);
                    if ((topo.aspect(nebs.idown(i), j) < -5 * (PI / 8)) && (topo.aspect(nebs.idown(i), j) > -3 * (PI / 8)))
                        incoming += W * I_R(nebs.idown(i), j);
                }
                if (NE > 0) {
                    incoming += NE * I_P(i, j) * NE_Ip(i, j);
                    if ((topo.aspect(nebs.iup(i), nebs.jup(j)) < 7 * (PI / 8)) && (topo.aspect(nebs.iup(i), nebs.jup(j)) > 5 * (PI / 8)))
                        incoming += NE * I_R(nebs.iup(i), nebs.jup(j));
                }
                if (SE > 0) {
                    incoming += SE * I_P(i, j) * SE_Ip(i, j);
                    if ((topo.aspect(nebs.iup(i), nebs.jdown(j)) < 3 * (PI / 8)) && (topo.aspect(nebs.iup(i), nebs.jdown(j)) > 1 * (PI / 8)))
                        incoming += SE * I_R(nebs.iup(i), nebs.jdown(j));
                }
                if (SW > 0) {
                    incoming += SW * I_P(i, j) * SW_Ip(i, j);
                    if ((topo.aspect(nebs.idown(i), nebs.jdown(j)) < -1 * (PI / 8)) && (topo.aspect(nebs.idown(i), nebs.jdown(j)) > -3 * (PI / 8)))
                        incoming += SW * I_R(nebs.idown(i), nebs.jdown(j));
                }
                if (NW > 0) {
                    incoming += NW * I_P(i, j) * NW_Ip(i, j);
                    if ((topo.aspect(nebs.idown(i), nebs.jup(j)) < -5 * (PI / 8)) && (topo.aspect(nebs.idown(i), nebs.jup(j)) > -7 * (PI / 8)))
                        incoming += NW * I_R(nebs.idown(i), nebs.jup(j));
                }

                // save incoming_watts to be applied during avalanche
                incoming_watts(i, j) = incoming;
            }
        }
    }
}

void RadiationModel::save_rasters(std::string prefix) {
    shade_raster.save(prefix + "_shade_raster.asc");
    I_P.save(prefix + "_I_P.asc");
    incoming_watts.save(prefix + "_incoming.asc");
}
