#ifndef _MFD_FLOW_ROUTE_
#define _MFD_FLOW_ROUTE_

#include <vector>
#include "raster.h"
#include "grid_neighbours.h"
#include "global_defs.h"

/// \brief Multiple flow direction flow routing
class MFDFlowRouter {
    private:
        int size_x;  ///< Number of cells in the x dimension
        int size_y;  ///< Number of cells in the y dimension
        Raster fa_bounds;  ///< Raster for flow coming in at the boundaries

    public:
        /// \brief Create an MFDFlowRouter object
        MFDFlowRouter();

        /// \brief Initialise the MFDFlowRouter object
        /// \param flow The flow accumulation Raster containing the initial (boundary) flow values
        void initialise(Raster& flow);
        
        /// \brief Do the flow routing
        /// \param topo The Raster of elevations
        /// \param flow The flow accumulation Raster that will contain the output flow values
        /// \param nebs GridNeighbours instance for neighbour indexing
        void run(Raster& topo, Raster& flow, GridNeighbours& nebs);
};

#endif
