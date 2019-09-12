#include <iostream>
#include <string>
#include "inih/INIReader.h"
#include "utility.h"
#include "parameters.h"

// set the default values here
Parameters::Parameters() : U(0.01), K(0.001), D(1.5), melt(250), timestep(1), printinterval(1),
        ann_timestep(timestep / 8760), thresholdarea(1e35), init_exposure_age(0),
        init_sed_track(2), init_veg(8), year(2010), day(145), hour(12), minute(0),
        end_year(2015), end_day(1), lattitude(0), longitude(0), stdmed(0), declination(0),
        altitude(0), azimuth(0), topo_file("topo.asc"), fa_file("FA.asc"),
        sed_file("SedThickness.asc"), fix_random_seed(false), save_topo(true), save_flow(false) {}


/// Load parameter values from a .INI file. Parameters can be omitted from the file,
/// in which case they will take their default value.
Parameters::Parameters(const std::string& parameter_file) : Parameters() {
    // load parameter file if it exists
    std::cout << "Loading parameter file: " << parameter_file << std::endl;
    INIReader reader(parameter_file);
    if (reader.ParseError() != 0) {
        Util::Error(std::string("Cannot load parameter file: ") + parameter_file, 1);
    }

    // Setup Key Model Variables
    set_U(reader.GetReal("model", "U", U));        // 'Uplift', m yr^-1
    set_K(reader.GetReal("model", "K", K));        // Stream Power, yr^-1
    set_D(reader.GetReal("model", "D", D));        // Diffusion, yr^-1
    set_melt(reader.GetReal("model", "melt", melt));    // Reciprocal melt rate, for a given radiation input

    set_timestep(reader.GetReal("time", "timestep", timestep));   // Time step in hours
    set_printinterval(reader.GetInteger("output", "printinterval", printinterval)); // Output timestep, in hours
    set_save_topo(reader.GetBoolean("output", "save_topo", save_topo));
    set_save_flow(reader.GetBoolean("output", "save_flow", save_flow));

//  thresh(0.577 * deltax;   // Critical height in m above neighbouring pixel, at 30 deg  (TAN(RADIANS(33deg))*deltax
//  thresh_diag(thresh * sqrt2;
    set_thresholdarea(reader.GetReal("model", "thresholdarea", thresholdarea));  // Threshold for diffusion domain - to prevent diffusion in channels, etc. (m2)

    set_init_exposure_age(reader.GetReal("model", "init_exposure_age", init_exposure_age));    // Variables used to initiate exposure time, sed depth and veg age rasters
    set_init_sed_track(reader.GetReal("model", "init_sed_track", init_sed_track));
    set_init_veg(reader.GetReal("model", "init_veg", init_veg));

    set_year(reader.GetInteger("time", "year", year));
    set_day(reader.GetInteger("time", "day", day));  // 144;              // May 25th is the start of melt/rain season
    set_hour(reader.GetInteger("time", "hour", hour));    // 24-hr clock
    set_minute(reader.GetInteger("time", "hour", minute));   // 0 in most cases
    set_end_year(reader.GetInteger("time", "end_year", end_year));  // Model execution ends on the first day of this year
    set_end_day(reader.GetInteger("time", "end_day", end_day));
//  ct(ModelTime(year, day, hour, minute, end_year, end_day));
//  duration(end_year - year;   // Model execution time, in years, keeping in mind melt season is 138 days

    set_lattitude(reader.GetReal("solar_geom", "lattitude", lattitude)); // 67.3;
    set_longitude(reader.GetReal("solar_geom", "longitude", longitude)); // 134.9;         // Dempster Coordinates
    set_stdmed(reader.GetReal("solar_geom", "stdmed", stdmed)); //9 * 15;         // Standard meridian of nearest time zone. LSTM((UTC - 9H * 15 deg) n.b. Alaska Time Meridian(135 deg W
    set_declination(reader.GetReal("solar_geom", "declination", declination));
    set_altitude(reader.GetReal("solar_geom", "altitude", altitude));
    set_azimuth(reader.GetReal("solar_geom", "azimuth", azimuth));

    // input file names
    set_topo_file(reader.Get("input", "topo", topo_file));
    set_fa_file(reader.Get("input", "FA", fa_file));
    set_sed_file(reader.Get("input", "sed", sed_file));

    // should we fix the random number seed (for testing)
    set_fix_random_seed(reader.GetBoolean("random", "fix_seed", fix_random_seed));
}

void Parameters::set_U(real_type U_) {
    U = U_;
}

void Parameters::set_altitude(real_type altitude_) {
    altitude = altitude_;
}

void Parameters::set_azimuth(real_type azimuth_) {
    azimuth = azimuth_;
}

void Parameters::set_K(real_type K_) {
    K = K_;
}

void Parameters::set_D(real_type D_) {
    D = D_;
}

void Parameters::set_day(int day_) {
    if (day_ < 1 || day_ > 365) {
        Util::Error("Trying to set day to out of range value", 1);
    }
    else {
        day = day_;
    }
}

void Parameters::set_year(int year_) {
    year = year_;
}

void Parameters::set_hour(int hour_) {
    if (hour_ < 1 || hour_ > 24) {
        Util::Error("Trying to set hour to out of range value", 1);
    }
    else {
        hour = hour_;
    }
}

void Parameters::set_minute(int minute_) {
    if (minute_ < 1 || minute_ > 60) {
        Util::Error("Trying to set minute to out of range value", 1);
    }
    else {
        minute = minute_;
    }
}

void Parameters::set_end_year(int end_year_) {
    if (end_year_ <  year) {
        Util::Error("Trying to set end_year to less than current year", 1);
    }
    else {
        end_year = end_year_;
    }
}

void Parameters::set_end_day(int end_day_) {
    if (end_day_ < 1 || end_day_ > 365) {
        Util::Error("Trying to set end_day to out of range value", 1);
    }
    else {
        end_day = end_day_;
    }
}

void Parameters::set_declination(real_type declination_) {
    declination = declination_;
}

void Parameters::set_fa_file(const std::string& fa_file_) {
    fa_file = fa_file_;
}

void Parameters::set_sed_file(const std::string& sed_file_) {
    sed_file = sed_file_;
}

void Parameters::set_topo_file(const std::string& topo_file_) {
    topo_file = topo_file_;
}

void Parameters::set_fix_random_seed(bool fix_random_seed_) {
    fix_random_seed = fix_random_seed_;
}

void Parameters::set_init_exposure_age(real_type init_exposure_age_) {
    init_exposure_age = init_exposure_age_;
}

void Parameters::set_init_sed_track(real_type init_sed_track_) {
    init_sed_track = init_sed_track_;
}

void Parameters::set_init_veg(real_type init_veg_) {
    init_veg = init_veg_;
}

void Parameters::set_lattitude(real_type lattitude_) {
    lattitude = lattitude_;
}

void Parameters::set_longitude(real_type longitude_) {
    longitude = longitude_;
}

void Parameters::set_melt(real_type melt_) {
    melt = melt_;
}

void Parameters::set_stdmed(real_type stdmed_) {
    stdmed = stdmed_;
}

void Parameters::set_printinterval(int printinterval_) {
    if (printinterval < 1) {
        Util::Error("Print interval must be greater than 0", 1);
    }
    else {
        printinterval = printinterval_;
    }
}

void Parameters::set_thresholdarea(real_type thresholdarea_) {
    thresholdarea = thresholdarea_;
}

void Parameters::set_timestep(real_type timestep_) {
    if (timestep <= 0) {
        Util::Error("Timestep must be greater than 0", 1);
    }
    else {
        timestep = timestep_;
        ann_timestep = timestep / 8760;
    }
}

void Parameters::set_save_topo(bool save_topo_) {
    save_topo = save_topo_;
}

void Parameters::set_save_flow(bool save_flow_) {
    save_flow = save_flow_;
}
