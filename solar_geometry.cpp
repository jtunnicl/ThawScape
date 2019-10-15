#include "global_defs.h"
#include "model_time.h"
#include "solar_geometry.h"


SolarGeometry::SolarGeometry() : latitude(0), longitude(0), stdmed(0),
        declination(0), altitude(0), azimuth(0), incidence(0), SHA(0) {}


SolarGeometry::SolarGeometry(Parameters& params) : SolarGeometry() {
	latitude = params.get_latitude();
	longitude = params.get_longitude();
	stdmed = params.get_stdmed();
	declination = params.get_declination();
	altitude = params.get_altitude();
	azimuth = params.get_azimuth();
}


/// Documentation for Sun Position and Solar influx come from the following sources:
/// ME 4131 THERMAL ENVIRONMENTAL ENGINEERING LABORATORY MANUAL, Appendix D
/// http://www.me.umn.edu/courses/me4131/LabManual/AppDSolarRadiation.pdf
///
/// Kumar, L., Skidmore, A. and Knowles, E. 1997. Modelling topographic variation
/// in solar radiation in a GIS environment. Int J. Information Science, 11(5), p.475-497
///
/// Angus, R., Muneer, T. 1993. Sun position for daylight models: Precise algorithms for
/// determination. Lighting Research and Technology 25(2) 81-83.
void SolarGeometry::sun_position(const ModelTime& ct) {
	real_type m, n, B, EOT, LST;
    int day = ct.get_day();
    int hour = ct.get_hour();
    int minute = ct.get_minute();

	// Solar declination: angle of suns rays relative to equato Max pos at summer equinox = +23.5
	declination = 23.45 * sin( 360. / 365. * (284. + day) * degrad);

	B = ( 360. / 364. ) * ( day - 81. ) * degrad;              // result in radians
	EOT = 0.165 * sin(2. * B ) - 0.126 * cos( B ) - 0.025 * sin( B );
	                           // Equation of Time [hr]
	LST = hour + ( minute/60 ) + ( stdmed - longitude ) / 15 + EOT - 0;     //  Last term is Daylight Saving (e.g. +1)
	                           // Local Solar Time, correcting for distance from nearest time zone meridian
	SHA = 15 * (LST - 12);   // Local Solar Hour  (negative before solar noon, positive after)

	m = sin(latitude * degrad) * sin(declination * degrad);
	n = cos(latitude * degrad) * cos(declination * degrad) * cos(SHA * degrad);
	altitude = asin(m + n) / degrad;

	m = sin(latitude * degrad) * cos(declination * degrad) * cos(SHA * degrad);
	n = cos(latitude * degrad) * sin(declination * degrad);
	azimuth = acos ( (m - n) / cos(altitude * degrad) ) / degrad - 180;   // in degrees

	if (SHA > 0) azimuth = -azimuth;              // N = 0; E = -90; S = -180/+180; W = +90
}
