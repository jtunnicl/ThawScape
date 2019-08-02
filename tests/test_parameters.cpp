#include <string>

#include "catch2/catch.hpp"
#include "parameters.h"


TEST_CASE("Parameters", "[parameters]") {
	Parameters p("ThawScapeTestInit.ini");

    // check parameters are correct
    REQUIRE(p.get_topo_file().compare(std::string("topotest.asc")) == 0);
    REQUIRE(p.get_fa_file().compare(std::string("FAtest.asc")) == 0);
    REQUIRE(p.get_U() == Approx(0.014));
    REQUIRE(p.get_K() == Approx(0.055));
    REQUIRE(p.get_D() == Approx(1.505));
    REQUIRE(p.get_melt() == Approx(255));
    REQUIRE(p.get_thresholdarea() == Approx(0.2));
    REQUIRE(p.get_init_exposure_age() == Approx(1));
    REQUIRE(p.get_init_sed_track() == Approx(4));
    REQUIRE(p.get_init_veg() == Approx(2));
    REQUIRE(p.get_timestep() == Approx(4));
    REQUIRE(p.get_ann_timestep() == Approx(4.0 / 8760.0));
    REQUIRE(p.get_printinterval() == 98);
    REQUIRE(p.get_year() == 2012);
    REQUIRE(p.get_day() == 74);
    REQUIRE(p.get_hour() == 11);
    REQUIRE(p.get_end_year() == 2014);
    REQUIRE(p.get_end_day() == 1);
    REQUIRE(p.get_lattitude() == Approx(1));
    REQUIRE(p.get_longitude() == Approx(3.5));
    REQUIRE(p.get_stdmed() == Approx(135));
    REQUIRE(p.get_declination() == Approx(1.0));
    REQUIRE(p.get_altitude() == Approx(0.3));
    REQUIRE(p.get_azimuth() == Approx(0.4));
}
