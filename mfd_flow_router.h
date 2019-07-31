#ifndef _MFD_FLOW_ROUTE_
#define _MFD_FLOW_ROUTE_

#include <vector>
#include "raster.h"
#include "global_defs.h"

class MFDFlowRouter {
    private:
        int size_x;
        int size_y;
        Raster flow1;
        Raster flow2;
        Raster flow3;
        Raster flow4;
        Raster flow5;
        Raster flow6;
        Raster flow7;
        Raster flow8;
        Raster fa_bounds;

    public:
        MFDFlowRouter();

        void initialise(Raster& initial_flow);
        
        void run(Raster& topo, Raster& flow, std::vector<int>& iup, std::vector<int>& idown,
                std::vector<int>& jup, std::vector<int>& jdown);  // expect this to change (having to pass lots of Rasters in)
};

#endif
