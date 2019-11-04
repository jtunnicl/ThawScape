#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <string>
#include "global_defs.h"

/// \brief Class for loading and accessing the parameters
///
/// Each parameter has its own private variable defined within this class and public "getter" and "setter"
/// methods that can be used to access the variable. Usually the setter is only used while reading
/// parameters from file and then the getter is used wherever the parameter is required in the code.
///
/// In order to add a new parameter one should:
///   -# add a new private variable for the parameter, for example `real_type my_new_param`
///   -# add a public method for setting the parameter, including any checks that the parameter is "sane",
///      for example `void Parameters::set_my_new_param(real_type my_new_param_) { my_new_param = my_new_param_; }`
///   -# add a public method for getting the parameter, for example
///      `real_type Parameters::get_my_new_param() const { return my_new_param; }`
///   -# add the new parameter with default value to the initialiser list in `Parameters::Parameters()` from
///      parameters.cpp
///   -# add a line to `Parameters::Parameters(const std::string& parameter_file)` to load the new parameter
///      from the input file, for example `set_my_new_param(reader.GetReal("example", "my_new_param", my_new_param));`
///      (assuming the new parameter is in the "example" section in the input file)
///   -# add the parameter into your input file (ThawScape.ini), for example
///      \code
///      [example]
///      my_new_param = 1.0
///      \endcode
class Parameters {
    private:
        real_type U;  ///< Uplift, m yr^-1
        real_type K;  ///< Stream Power, yr^-1
        real_type D;  ///< Diffusion, yr^-1
        real_type melt;  ///< Reciprocal melt rate, for a given radiation input
        real_type timestep;  ///< Time step in hours
        int printinterval;  ///< Output timestep, in hours
        real_type ann_timestep;  ///< Used in formula based on annual rates (e.g. 2 hours over 8760 hrs in 365 days)
        real_type thresholdarea;  ///< Threshold for diffusion domain - to prevent diffusion in channels, etc. (m2)
        real_type init_exposure_age;  ///< Used to initialise exposure time
        real_type init_sed_track;  ///< Used to initialise sed depth
        real_type init_veg;  ///< Used to to initialise vegetation layer
        int year;  ///< Initial year
        int day;  ///< Initial day
        int hour;  ///< Initial hour
        int minute;  ///< Initial minute
        int end_year;  ///< End year of the simulation
        int end_day;  ///< End day of the simulation (the simulation ends at the beginning of this day)
        real_type latitude;
        real_type longitude;
        real_type stdmed;  ///< Standard meridian of nearest time zone
        real_type declination;
        real_type altitude;
        real_type azimuth;
        std::string topo_file;  ///< Input file name for topo (elevations) raster
        std::string fa_file;  ///< Input file name for the flow accumulation raster
        std::string sed_file;  ///< Input file name for the sediment layer
        bool fix_random_seed;  ///< Fix the random seed
        bool save_topo;  ///< Save the topo (elevations) raster
        bool save_flow;  ///< Save the flow accumulation raster
        int flood_algorithm;  ///< Choose the algorithm for flood/pit-filling
        bool avalanche;  ///< Enable the avalanche component
        bool flood;  ///< Enable the flood component
        bool flow_routing;  ///< Enable the flow routing component
        bool diffusive_erosion;  ///< Enable the diffusive hillslope erosion component
        bool uplift;  ///< Enable the uplift component
        bool melt_component;  ///< Enable the melt component
        bool channel_erosion;  ///< Enable the channel erosion component
        bool debug_melt;   ///< Write out additional Rasters for debugging the RadiationModel

    public:
        /// \brief Default Parameters object
        Parameters();

        /// \brief Create Parameters object with values loaded from file
        Parameters(const std::string& parameter_file);

        void set_U(real_type U_) { U = U_; }
        real_type get_U() const { return U; }

        void set_K(real_type K_) { K = K_; }
        real_type get_K() const { return K; }

        void set_D(real_type D_) { D = D_; }
        real_type get_D() const { return D; }

        void set_melt(real_type melt_) { melt = melt_; }
        real_type get_melt() const { return melt; }

        void set_timestep(real_type timestep_);
        real_type get_timestep() const { return timestep; }
        real_type get_ann_timestep() const { return ann_timestep; }

        void set_printinterval(int printinterval_);
        int get_printinterval() const { return printinterval; }

        void set_thresholdarea(real_type thresholdarea_) { thresholdarea = thresholdarea_; }
        real_type get_thresholdarea() const { return thresholdarea; }

        void set_init_exposure_age(real_type init_exposure_age_) { init_exposure_age = init_exposure_age_; }
        real_type get_init_exposure_age() const { return init_exposure_age; }

        void set_init_sed_track(real_type init_sed_track_) { init_sed_track = init_sed_track_; }
        real_type get_init_sed_track() const { return init_sed_track; }

        void set_init_veg(real_type init_veg_) { init_veg = init_veg_; }
        real_type get_init_veg() const { return init_veg; }

        void set_year(int year_) { year = year_; }
        int get_year() const { return year; }

        void set_day(int day_);
        int get_day() const { return day; }

        void set_hour(int hour_);
        int get_hour() const { return hour; }

        void set_minute(int minute_);
        int get_minute() const { return minute; }

        void set_end_year(int end_year_);
        int get_end_year() const { return end_year; }

        void set_end_day(int end_day_);
        int get_end_day() const { return end_day; }

        void set_latitude(real_type latitude_) { latitude = latitude_; }
        real_type get_latitude() const { return latitude; }

        void set_longitude(real_type longitude_) { longitude = longitude_; }
        real_type get_longitude() const { return longitude; }

        void set_stdmed(real_type stdmed_) { stdmed = stdmed_; }
        real_type get_stdmed() const { return stdmed; }

        void set_declination(real_type declination_) { declination = declination_; }
        real_type get_declination() const { return declination; }

        void set_altitude(real_type altitude_) { altitude = altitude_; }
        real_type get_altitude() const { return altitude; }

        void set_azimuth(real_type azimuth_) { azimuth = azimuth_; }
        real_type get_azimuth() const { return azimuth; }

        void set_fa_file(const std::string& fa_file_) { fa_file = fa_file_; }
        std::string get_fa_file() const { return fa_file; }

        void set_topo_file(const std::string& topo_file_) { topo_file = topo_file_; }
        std::string get_topo_file() const { return topo_file; }

        void set_sed_file(const std::string& sed_file_) { sed_file = sed_file_; }
        std::string get_sed_file() const { return sed_file; }

        void set_fix_random_seed(bool fix_random_seed_) { fix_random_seed = fix_random_seed_; }
        bool get_fix_random_seed() const { return fix_random_seed; }

        void set_save_topo(bool save_topo_) { save_topo = save_topo_; }
        bool get_save_topo() const { return save_topo; }

        void set_save_flow(bool save_flow_) { save_flow = save_flow_; }
        bool get_save_flow() const { return save_flow; }

        void set_flood_algorithm(int flood_algorithm_);
        int get_flood_algorithm() const { return flood_algorithm; }

        void set_avalanche(bool avalanche_) { avalanche = avalanche_; }
        bool get_avalanche() const { return avalanche; }

        void set_flood(bool flood_) { flood = flood_; }
        bool get_flood() const { return flood; }

        void set_flow_routing(bool flow_routing_) { flow_routing = flow_routing_; }
        bool get_flow_routing() const { return flow_routing; }

        void set_diffusive_erosion(bool diffusive_erosion_) { diffusive_erosion = diffusive_erosion_; }
        bool get_diffusive_erosion() const { return diffusive_erosion; }

        void set_uplift(bool uplift_) { uplift = uplift_; }
        bool get_uplift() const { return uplift; }

        void set_melt_component(bool melt_component_) { melt_component = melt_component_; }
        bool get_melt_component() const { return melt_component; }

        void set_channel_erosion(bool channel_erosion_) { channel_erosion = channel_erosion_; }
        bool get_channel_erosion() const { return channel_erosion; }

        void set_debug_melt(bool debug_melt_) { debug_melt = debug_melt_; }
        bool get_debug_melt() const { return debug_melt; }
};

#endif
