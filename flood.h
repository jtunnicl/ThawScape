#ifndef _FLOOD_H_
#define _FLOOD_H_

#include "global_defs.h"
#include "dem.h"
#include "grid_neighbours.h"
#include "parameters.h"
#include "Array2D.hpp"


class Flood {
    private:
        int size_x;  ///< Number of cells in the x dimension
        int size_y;  ///< Number of cells in the y dimension
        DEM& topo;  ///< Raster of elevations
        GridNeighbours& nebs;  ///< Grid neighbour indexing
        Array2D<real_type> elevation;
        bool initialised;
        int algorithm;  ///< Which algorithm to use

        /// \brief Run one of Barnes' flood algorithms
        void run_barnes_flood();

        /// \brief Run Pelletier's algorithm
        void run_fillinpitsandflats();

        /// \brief Do Pelletier's pit filling
        void fillinpitsandflats(int i, int j);

    public:
        /// \brief Create a Flood object
        /// \param topo_ The Raster of elevations
        /// \param nebs_ GridNeighbours instance for neighbour indexing
        Flood(DEM& topo_, GridNeighbours& nebs_);

        /// \brief Initialise the MFDFlowRouter object
        /// \param params Parameters object
        void initialise(Parameters& params);
        
        /// \brief Run the algorithm
        void run();
};

#endif
