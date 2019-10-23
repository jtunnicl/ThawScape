#ifndef _FLOOD_H_
#define _FLOOD_H_

#include "global_defs.h"
#include "raster.h"
#include "grid_neighbours.h"
#include "parameters.h"
#include "Array2D.hpp"


class Flood {
    private:
        int size_x;  ///< Number of cells in the x dimension
        int size_y;  ///< Number of cells in the y dimension
        Array2D<real_type> elevation;  ///< Elevation array for passing to Barnes' routines
        int algorithm;  ///< Which algorithm to use

        /// \brief Run one of Barnes' flood algorithms
        void run_barnes_flood(Raster& topo);

        /// \brief Run Pelletier's algorithm
        void run_fillinpitsandflats(Raster& topo, GridNeighbours& nebs);

        /// \brief Do Pelletier's pit filling
        void fillinpitsandflats(int i, int j, Raster& topo, GridNeighbours& nebs);

    public:
        /// \brief Create a Flood object
        Flood();

        /// \brief Initialise the Flood object
        /// \param topo The Raster of elevations
        /// \param params Parameters object
        void initialise(Raster& topo, Parameters& params);
        
        /// \brief Run the algorithm
        /// \param topo The Raster of elevations
        /// \param nebs GridNeighbours instance for neighbour indexing
        void run(Raster& topo, GridNeighbours& nebs);
};

#endif
