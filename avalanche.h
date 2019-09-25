#ifndef _AVALANCHE_H_
#define _AVALANCHE_H_

#include "grid_neighbours.h"
#include "dem.h"
#include "global_defs.h"

/// \brief Avalanching
class Avalanche {
    private:
        int size_x;  ///< Number of cells in the x dimension
        int size_y;  ///< Number of cells in the y dimension
        DEM& topo;  ///< Raster of elevations
        Raster& sed_track;  ///< Raster of sediment track depth
        GridNeighbours& nebs;  ///< Grid neighbour indexing
        real_type thresh;  ///< Critical height in m above neighbouring pixel
        real_type thresh_diag;  ///< Critical height in m above neighbouring pixel along a diagonal
        bool initialised;

    public:
        /// \brief Create an Avalanche object
        /// \param topo_ The Raster of elevations
        /// \param sed_track_ The Raster of sediment track depth
        /// \param nebs_ GridNeighbours instance for neighbour indexing
        Avalanche(DEM& topo_, Raster& sed_track_, GridNeighbours& nebs_);

        /// \brief Initialise the Avalanche object
        void initialise();
        
        /// \brief Run the avalanche code
        void run();
};

#endif
