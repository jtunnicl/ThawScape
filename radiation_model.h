#ifndef _RADIATION_MODEL_H
#define _RADIATION_MODEL_H

#include "raster.h"
#include "grid_neighbours.h"
#include "global_defs.h"
#include "solar_geometry.h"
#include "dem.h"
#include "parameters.h"
#include "model_time.h"

/// \brief RadiationModel class for carrying out melt
class RadiationModel {
    private:
        int lattice_size_x;
        int lattice_size_y;
        real_type deltax;
        real_type deltax2;
        real_type melt;   ///< Melt parameter
        SolarGeometry r;  ///< Solar geometry
        Raster solar_raster;
        Raster shade_raster; 
        Raster I_D;
        Raster I_R;
        Raster I_P;
        Raster N_Ip;
        Raster E_Ip;
        Raster S_Ip;
        Raster W_Ip;
        Raster NE_Ip;
        Raster SE_Ip;
        Raster SW_Ip;
        Raster NW_Ip;
        Raster Ip_D8;  ///< Map of incoming solar flux, 8 directions

        /// \brief Compute solar influx
        /// \param topo The DEM
        /// \param ct ModelTime object
        void solar_influx(DEM& topo, ModelTime& ct);

    public:
        /// \brief Create RadiationModel
        RadiationModel();

        /// \brief Initialise the RadiationModel object
        /// \param topo The DEM
        /// \param params Parameters object
        void initialise(DEM& topo, Parameters& params);

        /// \brief Update the solar characteristics
        /// \param topo The DEM
        /// \param ct Current ModelTime
        void update_solar_characteristics(DEM& topo, ModelTime& ct);

        /// \brief Melt exposed ice
        /// \param topo The DEM
        /// \param Sed_Track Sediment track depth Raster
        /// \param flow The flow accumulation Raster
        /// \param nebs Grid neighbour indexing
        void melt_exposed_ice(DEM& topo, Raster& Sed_Track, Raster& flow, GridNeighbours& nebs);

        /// \brief Get the current solar altitude
        /// \returns altitude The current solar altitude
        real_type get_solar_altitude() { return r.get_altitude(); }
};

#endif
