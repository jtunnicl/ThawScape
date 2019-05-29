# ThawScape

ThawScape is experimental code to simulate the influence of thermokarst dynamics upon ice-cored terrain.
The code uses some well-established landscape evolution algorithms (cf. Pelletier, 2008) to simulate water and sediment flux over a landscape that is underlain by many (10-100s) meters of ground ice (sediment/ice ~= 15/85). As the cover soil layer is eroded to expose the ice, melt proceeds according to sun conditions and precipitation. Meltwater contributes to the overland and channelized flow component. If the thermokarst is covered in sediment, the melt/ablation process will cease until the ice is exposed again.

Sediment transport occurs both by hillslope diffusion ('D') and advective channel transport ('K'). Uplift ('U') is included, to maintain erosive potential, particularly over long timescales. Diffusive transport is higher within saturated sediment that runs off from melting thermokarst erosion sites. The ice melt rate ('melt') is calibrated to observed rates of thermokarst retreat, and requires some further experimentation to link values to physical (field) variables. A vegetation layer has been included to modify the erosion threshold and diffusion rate; vegetation will gradually re-establish after 3 years of no erosion activity.

Input requirements are (1) an ascii, arc-type DEM (e.g. topo.asc) with 6-line header indicating row and columns, lower left corner, cell resolution, and NoData values. (2) a flow accumulation raster (e.g. FA.asc), built with a multi-flow direction (MFD) algorithm. The code uses the boundary elements of this raster to determine flow contributions from outside the model domain. (3) Finally, an input file (e.g. ThawScape.ini) is used to specify model parameters.

Output is written as a sequence of ascii rasters that can be read with any GIS package.

## Building and Running ThawScape

ThawScape uses the CMake build system and requires a C++11 compiler.

```
mkdir build && cd build     # create a build directory and switch to it
cmake ..                    # run CMake to configure the build
make                        # build the code
ctest                       # run the tests
```

After the build has completed there should be a `ThawScape` executable in the
build directory.

The `ThawScape.ini` input file must be present in the directory you run the
executable from.
