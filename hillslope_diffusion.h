#ifndef _HILL_SLOPE_DIFFUSION_H_
#define _HILL_SLOPE_DIFFUSION_H_

#include "raster.h"
#include "grid_neighbours.h"
#include "parameters.h"
#include "global_defs.h"

class HillSlopeDiffusion {
    private:
        Raster topoold;  ///< Raster containing old elevations
        int lattice_size_x;   ///< x dimension
        int lattice_size_y;   ///< y dimension
        real_type D;   ///< Diffusion rate
        real_type deltax2;   ///< Pixel size squared
        real_type ann_timestep;   ///< Time step in years
        real_type thresholdarea;
        real_vector ax;
        real_vector ay;
        real_vector bx;
        real_vector by;
        real_vector cx;
        real_vector cy;
        real_vector ux;
        real_vector uy;
        real_vector rx;
        real_vector ry;

        /// \brief Solve a tridiagonal system
        void tridag(real_vector& a, real_vector& b, real_vector& c, real_vector& r, real_vector& u, int n);

    public:
        /// \brief Create HillSlopeDiffusion object
        HillSlopeDiffusion();

        /// \brief Initialise the HillSlopeDiffusion object
        /// \param topo Flow accumulation Raster
        /// \param params Parameters object
        void initialise(Raster& topo, Parameters& params);

        /// \brief Run the HillSlopeDiffusion algorithm
        /// \param topo Elevations Raster
        /// \param flow Flow accumulation Raster
        /// \param nebs Neighbour indexing object
        void run(Raster& topo, Raster& flow, GridNeighbours& nebs);
};

#endif
