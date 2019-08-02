#ifndef _HILL_SLOPE_DIFFUSION_H_
#define _HILL_SLOPE_DIFFUSION_H_

#include "raster.h"
#include "grid_neighbours.h"
#include "parameters.h"
#include "global_defs.h"

class HillSlopeDiffusion {
    private:
        Raster& topo;  ///< Raster of elevations
        Raster& flow;  ///< Flow accumulation raster
        GridNeighbours& nebs;  ///< Grid neighbour indexing
        Parameters& params;  ///< Parameters object
        Raster topoold;  ///< Raster containing old elevations
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
        /// \param topo_ Elevations Raster, stored as a reference
        /// \param topo_ Flow accumulation Raster, stored as a reference
        /// \param nebs_ Neighbour indexing object, stored as a reference
        /// \param params_ Parameters object, stored as a reference
        HillSlopeDiffusion(Raster& topo_, Raster& flow_, GridNeighbours& nebs_, Parameters& params_);

        /// \brief Run the HillSlopeDiffusion algorithm
        void run();
};

#endif
