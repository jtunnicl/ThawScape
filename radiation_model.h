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
        DEM& topo;  ///< Reference to topo DEM
        Raster& Sed_Track;  /// Raster of sediment track depth
        Raster& flow;  ///< Reference to flow Raster
        GridNeighbours& nebs;  ///< Reference to grid neighbour indexing
        Parameters& params;  ///< References to parameters object
        SolarGeometry r;  ///< Solar geometry
        bool initialised;  ///< Whether the RadiationModel has been initialised
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
        void solar_influx(ModelTime& ct);

    public:
        /// \brief Create RadiationModel
        /// \param topo_ The DEM
        /// \param flow_ The flow accumulation Raster
        /// \param nebs_ Grid neighbour indexing
        /// \param params Parameters object
        RadiationModel(DEM& topo_, Raster& Sed_Track_, Raster& flow, GridNeighbours& nebs_, Parameters& params);

        /// \brief Initialise the RadiationModel object
        void initialise();

        /// \brief Update the solar characteristics
        void update_solar_characteristics(ModelTime& ct);

        /// \brief Melt exposed ice
        void melt_exposed_ice();

        /// \brief Get the current solar altitude
        /// \returns altitude The current solar altitude
        real_type get_solar_altitude() { return r.get_altitude(); }
};

#endif
