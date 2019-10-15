#ifndef _SOLAR_GEOMETRY_
#define _SOLAR_GEOMETRY_

#include "global_defs.h"
#include "model_time.h"
#include "parameters.h"

class SolarGeometry {
    private:
        real_type latitude;
        real_type longitude;
        real_type stdmed;          ///< LSTM = (UTC - 7H * 15 deg)
        real_type declination;     ///< Declination of sun from equatorial plane
        real_type altitude;        ///< Sun altitude in the sky
        real_type azimuth;         ///< Compass angle of sun
        real_type incidence;       ///< Angle of sun's incidence
        real_type SHA;             ///< \brief Solar Hour Angle
                                   ///<
                                   ///< Solar Hour Angle is 0° at solar noon.  Since the Earth rotates 15° per hour,
                                   ///< each hour away from solar noon corresponds to an angular motion of the sun in the sky of 15°.
                                   ///< In the morning the hour angle is negative, in the afternoon the hour angle is positive.

    public:
        /// \brief Create SolarGeometry object
        SolarGeometry();

        /// \brief Create SolarGeometry object initialised with Parameters
        SolarGeometry(Parameters& params);

        /// \brief Compute the sun position
        void sun_position(const ModelTime& ct);

        /// \brief Returns the latitude
        real_type get_latitude() const { return latitude; }

        /// \brief Returns the longitude
        real_type get_longitude() const { return longitude; }

        /// \brief Returns stdmed
        real_type get_stdmed() const { return stdmed; }

        /// \brief Returns the declination of sun from equatorial plane
        real_type get_declination() const { return declination; }

        /// \brief Returns the sun's altitude in the sky
        real_type get_altitude() const { return altitude; }

        /// \brief Returns compass angle of the sun
        real_type get_azimuth() const { return azimuth; }

        /// \brief Returns the angle of sun's incidence
        real_type get_incidence() const { return incidence; }

        /// \brief Returns the solar hour angle
        real_type get_SHA() const { return SHA; }
};

#endif
