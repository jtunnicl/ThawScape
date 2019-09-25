#include <cmath>
#include <algorithm>
#include "raster.h"
#include "grid_neighbours.h"
#include "global_defs.h"
#include "solar_geometry.h"
#include "dem.h"
#include "model_time.h"
#include "radiation_model.h"


RadiationModel::RadiationModel(DEM& topo_, Raster& Sed_Track_, Raster& flow_, GridNeighbours& nebs_, Parameters& params_) :
        topo(topo_), Sed_Track(Sed_Track_), flow(flow_), nebs(nebs_), params(params_), r(params_), initialised(false) {}


/// Allocates all Rasters to the correct size. This should be called before running the radiation model.
void RadiationModel::initialise() {
    int lattice_size_x = topo.get_size_x();
    int lattice_size_y = topo.get_size_y();
	solar_raster = Raster(lattice_size_x, lattice_size_y, 0.0);
	shade_raster = Raster(lattice_size_x, lattice_size_y);
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

    initialised = true;
}

void RadiationModel::update_solar_characteristics(ModelTime& ct) {
    r.sun_position(ct);
    solar_influx(ct);
}

void RadiationModel::solar_influx(ModelTime& ct) {
    int lattice_size_x = topo.get_size_x();
    int lattice_size_y = topo.get_size_y();

    // Calculate shading from surrounding terrain
	int i, j, m;
	real_type m1, m2, m3, m4, m5, M, asp360, d80;
	real_type cos_i, I_o, tau_b;
	std::vector<real_type> cos_i80, asp_4;

	// Solar position, in radians
	real_type azm = r.get_azimuth() * degrad;            //  Anything in the 'r' object uses degrees; converted here to radians
	real_type alt = r.get_altitude() * degrad;           //  Invert cos<>sin to obtain zenith angle. 0 is sun at zenith (flat terrain faces up); 90 sun is at the horizon (vertical terrain)
	real_type lat = r.get_lattitude() * degrad;
	real_type dec = r.get_declination() * degrad;
	real_type sha = r.get_SHA() * degrad;

	cos_i80 = std::vector<real_type>(8);
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

			for (m = 0; m < 8; m++)
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

void RadiationModel::melt_exposed_ice() {
    int lattice_size_x = topo.get_size_x();
    int lattice_size_y = topo.get_size_y();
    real_type deltax = topo.get_deltax();
    real_type deltax2 = deltax * deltax;

    // sort by elevations
    topo.sort_data();

    int t = lattice_size_x * lattice_size_y;
    while (t > 0)
    {
        t--;
        int i, j;
        topo.get_sorted_ij(t, i, j);
        if ((i > 0) && (i < (lattice_size_x - 1)) && (j > 0) && (j < (lattice_size_y - 1)))  // Do not alter boundary elements
        {
            real_type N, E, S, W, NE, SE, SW, NW;
            real_type incoming = 0;
            real_type elev_drop = 0;       // Decrease in elevation at central pixel, following ice melt
            real_type accommodation = 0;   // Volume available to fill below central pixel, in the immediate neighbourhood
            real_type lowestpixel;         // Elevation of the lowest pixel in the 9-element neighbourhood.
            int m = 0;

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

                // Ice mass lost, based on ablation at each face
                // incoming watts / meltrate / pixel area

                for (m = 0; m < 8; m++) {
                    if (topo(i, j) - neighb[m] > 0) accommodation += deltax2 * (topo(i, j) - neighb[m]);   // sum up all the volume available on pixels below the central pixel
                }

                elev_drop = incoming / params.get_melt() / deltax2;
                if (elev_drop * deltax < accommodation)           // i.e. There is room to accommodate the failed mass in neighbouring cells
                    topo(i, j) -= elev_drop;
                else
                    topo(i, j) = lowestpixel;


                // Water lost in melt flows downstream (add to 'flow' raster)
            }
        }
    }
}
