#ifndef _AVALANCHE_H_
#define _AVALANCHE_H_

#include "grid_neighbours.h"
#include "raster.h"
#include "global_defs.h"

/// \brief Avalanching
class Avalanche {
    private:
        int size_x;  ///< Number of cells in the x dimension
        int size_y;  ///< Number of cells in the y dimension
        real_type thresh;  ///< Critical height in m above neighbouring pixel
        real_type thresh_diag;  ///< Critical height in m above neighbouring pixel along a diagonal

    public:
        /// \brief Create an Avalanche object
        Avalanche();

        /// \brief Initialise the Avalanche object
        /// \param topo The Raster of elevations
        void initialise(Raster& topo);
        
        /// \brief Run the avalanche code
        /// \param topo The Raster of elevations
        /// \param sed_track The Raster of sediment track depth
        /// \param nebs GridNeighbours instance for neighbour indexing
        void run(Raster& topo, Raster& sed_track, GridNeighbours& nebs);
};

#endif
