#include <string>

#include "catch2/catch.hpp"
#include "streampower.h"


TEST_CASE("StreamPower Init", "[streampower]") {
	StreamPower sp = StreamPower(10, 10);
	sp.Init("ThawScapeTestInit.ini");

    // check parameters are correct
    REQUIRE(sp.topo_file.compare(std::string("topotest.asc")) == 0);
    REQUIRE(sp.fa_file.compare(std::string("FAtest.asc")) == 0);
    REQUIRE(sp.U == Approx(0.014));
    REQUIRE(sp.K == Approx(0.055));
    REQUIRE(sp.D == Approx(1.505));
    REQUIRE(sp.melt == Approx(255));
    REQUIRE(sp.deltax == Approx(10.5));
    REQUIRE(sp.nodata == Approx(-9999.5));
    REQUIRE(sp.xllcorner == Approx(5));
    REQUIRE(sp.yllcorner == Approx(2));
    REQUIRE(sp.thresholdarea == Approx(0.2));
    REQUIRE(sp.init_exposure_age == Approx(1));
    REQUIRE(sp.init_sed_track == Approx(4));
    REQUIRE(sp.init_veg == Approx(2));
    REQUIRE(sp.timestep == Approx(4));
    REQUIRE(sp.printinterval == 98);
    REQUIRE(sp.ct.year == 2012);
    REQUIRE(sp.ct.day == 74);
    REQUIRE(sp.ct.hour == 11);
    REQUIRE(sp.ct.end_year == 2014);
    REQUIRE(sp.r.lattitude == Approx(1));
    REQUIRE(sp.r.longitude == Approx(3.5));
    REQUIRE(sp.r.stdmed == Approx(135));
    REQUIRE(sp.r.declination == Approx(1.0));
    REQUIRE(sp.r.altitude == Approx(0.3));
    REQUIRE(sp.r.azimuth == Approx(0.4));
}
