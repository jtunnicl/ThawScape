#ifndef _DEM_H_
#define _DEM_H_

#include "global_defs.h"
#include "grid_neighbours.h"
#include "raster.h"

/// \brief The DEM class adds methods for computing the slope and aspect of a Raster
class DEM : public Raster {
    private:
        Raster slope_;  ///< The slope at each point in the Raster of elevations
        Raster aspect_;  ///< The aspect of each pixel in the Raster of elevations

    public:
        /// \brief Create an empty DEM object
        DEM();

        /// \brief Create a DEM object with the given dimensions and all elements initialised to the same value
        /// \param size_x_ x dimension of the DEM
        /// \param size_y_ y dimension of the DEM
        /// \param value The value to initialise the DEM with
        DEM(int size_x_, int size_y_, real_type value);

        /// \brief Create a DEM object by loading it from a file
        /// \param filename The name of the file to load the DEM from
        DEM(const std::string &filename);

        /// \brief Compute the slope and aspect of the DEM
        void compute_slope_and_aspect(const GridNeighbours& nebs);

        /// \brief Get the slope at a point in the DEM
        /// \param i x index of the pixel to return the slope at
        /// \param j y index of the pixel to return the slope at
        /// \returns slope the slope at the given pixel
        real_type slope(const int i, const int j) const;

        /// \brief Get the aspect at a point in the DEM
        /// \param i x index of the pixel to return the aspect at
        /// \param j y index of the pixel to return the aspect at
        /// \returns aspect the aspect at the given pixel
        real_type aspect(const int i, const int j) const;
};

#endif
