#include <vector>
#include "raster.h"
#include "grid_neighbours.h"
#include "parameters.h"
#include "global_defs.h"
#include "hillslope_diffusion.h"

HillSlopeDiffusion::HillSlopeDiffusion(Raster& topo_, Raster& flow_, GridNeighbours& nebs_,
        Parameters& params_) : topo(topo_), flow(flow_), nebs(nebs_), params(params_) {}


void HillSlopeDiffusion::run() {
    int lattice_size_x = topo.get_size_x();
    int lattice_size_y = topo.get_size_y();

    // only reallocate if not correct size
    if (ax.size() != lattice_size_x) {
        ax = real_vector(lattice_size_x);
        bx = real_vector(lattice_size_x);
        cx = real_vector(lattice_size_x);
        ux = real_vector(lattice_size_x);
        rx = real_vector(lattice_size_x);
    }
    if (ay.size() != lattice_size_y) {
        ay = real_vector(lattice_size_y);
        by = real_vector(lattice_size_y);
        cy = real_vector(lattice_size_y);
        uy = real_vector(lattice_size_y);
        ry = real_vector(lattice_size_y);
    }
    topoold.resize(topo.get_size_x(), topo.get_size_y());

    real_type D = params.get_D();
    real_type deltax2 = topo.get_deltax() * topo.get_deltax();
	int count = 0;
	while (count < 5)
	{
		count++;
        #pragma omp parallel for
		for (int i = 0; i < lattice_size_x; i++)
			for (int j = 0; j < lattice_size_y; j++)
				topoold(i, j) = topo(i, j);
		for (int i = 0; i < lattice_size_x; i++)
		{
			for (int j = 0; j < lattice_size_y; j++)
			{
				real_type term1 = D * params.get_ann_timestep() / (deltax2);
				if (flow(i, j) < params.get_thresholdarea())
				{
					ay[j] = -term1;
					cy[j] = -term1;
					by[j] = 4 * term1 + 1;
					ry[j] = term1 * ( topo(nebs.iup(i), j) + topo(nebs.idown(i), j) ) + topoold(i, j);
				}
				else
				{
					by[j] = 1;
					ay[j] = 0;
					cy[j] = 0;
					ry[j] = topoold(i, j);
				}
				if (j == 0)
				{
					by[j] = 1;
					cy[j] = 0;
					ry[j] = topoold(i, j);
				}
				if (j == lattice_size_y-1)
				{
					by[j] = 1;
					ay[j] = 0;
					ry[j] = topoold(i, j);
				}
			}
			tridag(ay, by, cy, ry, uy, lattice_size_y);
            #pragma omp parallel for
			for (int j = 0; j < lattice_size_y; j++)
				topo(i, j) = uy[j];
		}
        #pragma omp parallel for
		for (int i = 0; i < lattice_size_x; i++)
			for (int j = 0; j < lattice_size_y; j++)
				topoold(i, j) = topo(i, j);
		for (int j = 0; j < lattice_size_y; j++)
		{
			for (int i = 0; i < lattice_size_x; i++)
			{
				real_type term1 = D * params.get_timestep() / ( deltax2 );
				if (flow(i, j) < params.get_thresholdarea())
				{
					ax[i] = -term1;
					cx[i] = -term1;
					bx[i] = 4 * term1 + 1;
					rx[i] = term1 * ( topo(i, nebs.jup(j)) + topo(i, nebs.jdown(j)) ) + topoold(i, j);
				}
				else
				{
					bx[i] = 1;
					ax[i] = 0;
					cx[i] = 0;
					rx[i] = topoold(i, j);
				}
				if (i == 0)
				{
					bx[i] = 1;
					cx[i] = 0;
					rx[i] = topoold(i, j);
				}
				if (i == lattice_size_x-1)
				{
					bx[i] = 1;
					ax[i] = 0;
					rx[i] = topoold(i, j);
				}
			}
			tridag(ax, bx, cx, rx, ux, lattice_size_x);
            #pragma omp parallel for
			for (int i = 0; i < lattice_size_x; i++)
				topo(i, j) = ux[i];
		}
	}
}

void HillSlopeDiffusion::tridag(real_vector& a, real_vector& b, real_vector& c, real_vector& r, real_vector& u, int n) {
	unsigned long j;
	real_type bet;
	real_vector gam(n);

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
