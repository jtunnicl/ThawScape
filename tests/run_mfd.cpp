#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include "raster.h"
#include "dem.h"
#include "flood.h"
#include "mfd_flow_router.h"
#include "grid_neighbours.h"
#include "parameters.h"


int main(int argc, char** argv) {
    std::cout << "Running MFD flow routing" << std::endl;

    // arg parsing, could be done better, may change
    // usage: run_mfd <parameter_file>
    if (argc != 2) {
        std::cerr << "Must specify a parameter file" << std::endl;
        return 1;
    }
    // load parameters
    Parameters params(argv[1]);
    std::string topo_file = params.get_topo_file();
    std::string flow_file = params.get_fa_file();

    // load topo
    std::cout << "Topo file = " << topo_file << std::endl;
    DEM topo(topo_file);

    // flow accumulation raster
    std::cout << "Flow accumulation file = " << flow_file << std::endl;
    Raster flow(flow_file);

    // setup grid neighbour indexing
    GridNeighbours nebs;
    nebs.setup(topo.get_size_x(), topo.get_size_y());

    // flood object
    Flood flood;
    flood.initialise(topo, params);

    // flow routing object
    MFDFlowRouter flow_router;
    flow_router.initialise(flow);

    // save initial flow
    flow.save("flow0input.asc");

    // topo must be sorted before calling mfd and, since it won't change, we sort it once now
    topo.sort_data();

    // run the flood algorithm
    flood.run(topo, nebs);

    // do the flow routing (TODO: can it just be run by itself with nothing else?)
    flow_router.run(topo, flow, nebs);

    // save the flow
    flow.save("flow1output.asc");

    return 0;
}
