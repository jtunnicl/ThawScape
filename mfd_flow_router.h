#ifndef _MFD_FLOW_ROUTE_
#define _MFD_FLOW_ROUTE_

#include <vector>
#include "raster.h"
#include "global_defs.h"

/// \brief Multiple flow direction flow routing
class MFDFlowRouter {
    private:
        int size_x;  ///< Number of cells in the x dimension
        int size_y;  ///< Number of cells in the y dimension
        Raster& topo;  ///< Raster of elevations
        Raster& flow;  ///< Flow accumulation Raster
        std::vector<int>& iup; ///< Grid neighbour index vector
        std::vector<int>& idown; ///< Grid neighbour index vector
        std::vector<int>& jup; ///< Grid neighbour index vector
        std::vector<int>& jdown; ///< Grid neighbour index vector
        Raster flow1;  // these can probably be removed...
        Raster flow2;
        Raster flow3;
        Raster flow4;
        Raster flow5;
        Raster flow6;
        Raster flow7;
        Raster flow8;
        Raster fa_bounds;  ///< Raster for flow coming in at the boundaries

    public:
        /// \brief Create an MFDFlowRouter object and initialise it (this interface may change)
        /// \param topo_ The Raster of elevations
        /// \param flow_ The flow accumulation Raster that will contain the initial flow values
        /// \param iup_ Grid neighbour index vector
        /// \param idown_ Grid neighbour index vector
        /// \param jup_ Grid neighbour index vector
        /// \param jdown_ Grid neighbour index vector
        MFDFlowRouter(Raster& topo_, Raster& flow_, std::vector<int>& iup_, std::vector<int>& idown_,
                std::vector<int>& jup_, std::vector<int>& jdown_);  // interface may change

        /// \brief Initialise the MFDFlowRouter object
        void initialise();
        
        /// \brief Do the flow routing
        void run();
};

#endif
