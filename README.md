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

After the build has completed there should be a *ThawScape* executable in the
build directory.

The `ThawScape.ini` input file must be present in the directory you run the
executable from.

The default is to build a double precision version of ThawScape. If you would
like to use single precision then add the option `-DDOUBLE_PRECISION=OFF`.

## Code layout

The code is driven from *main.cpp*, which creates a `StreamPower` object (from
*streampower.cpp*):

- Inputs are loaded and components initialised during `StreamPower::Init()`
  - Parameters are loaded from the input file, *ThawScape.ini* (`Parameters()` from
    *parameters.cpp*)
  - DEM / `topo` Raster is loaded (`StreamPower::SetTopo()` from *streampower.cpp*)
  - Flow accumulation / `flow` Raster is loaded (`StreamPower::SetFA()` from
    *streampower.cpp*)
  - Most components need to be initialised with some combination of the DEM and
    flow `Raster`s and `Parameters` object before they can be used
  - Diffusion is initialised (`StreamPower::InitDiffusion()` from *streampower.cpp*)
- The main simulation loop occurs during `StreamPower::Start()`

The components are mainly split into their own files/classes as discussed below:

- Radiation model (`RadiationModel` from *radiation_model.cpp*)
  1. Compute slope and aspect of the DEM, required for radiation model calculations
     (`Raster::compute_slope_and_aspect()` from *raster.cpp*)
  2. Update solar characteristics (`RadiationModel::update_solar_characteristics()`
     from *radiation_model.cpp*)
     1. Update sun position (`SolarGeometry::sun_position()` from *solar_geometry.cpp*)
     2. Compute solar influx (`RadiationModel::solar_influx` from *radiation_model.cpp*)
  3. Compute melt potential (`RadiationModel::melt_potential()` from *radiation_model.cpp*)
     - **Note**: application of the melt potential occurs later during avalanching
- Hydro processes
  1. Flood / fill in pits and flats in the DEM, required before sorting points by elevation
     (`Flood::run()` from *flood.cpp*)
  2. Sort points in the elevation model by height, required as flow routing must proceed
     from high to low (`Raster::sort_data()` from *raster.cpp*)
  3. Calculate flow routing / flow accumulation (`MFDFlowRouter::run()` from
     *mfd_flow_router.cpp*)
  4. Compute slope and aspect of the DEM, required for channel erosion
     (`Raster::compute_slope_and_aspect()` from *raster.cpp*)
  5. Compute channel erosion based on flow accumulation (`StreamPower::channel_erosion()` 
     from *streampower.cpp*)
- Erosion
  1. Flood / fill in pits and flats in the DEM, required before sorting points by elevation
     (`Flood::run()` from *flood.cpp*)
  2. Sort points in the elevation model by height, required as avalanching must proceed
     from low to high (`Raster::sort_data()` from *Raster.cpp*)
  3. Landsliding / avalanching (`Avalanche::run()` from *avalanche.cpp*)
     - **Note**: The melt potential, computed earlier, is applied during this routine
  4. Diffusive hill slope erosion (`HillslopeDiffusion::run()` from *hillslope_diffusion.cpp*)
- Apply uplift (`StreamPower::uplift()` from *streampower.cpp*)

More developer documentation here: https://jtunnicl.github.io/ThawScape
