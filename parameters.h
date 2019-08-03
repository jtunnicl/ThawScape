#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <string>
#include "global_defs.h"

/// \brief Input Parameters
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
        real_type lattitude;
        real_type longitude;
        real_type stdmed;  ///< Standard meridian of nearest time zone
        real_type declination;
        real_type altitude;
        real_type azimuth;
        std::string topo_file;  ///< Input file name for topo (elevations) raster
        std::string fa_file;  ///< Input file name for the flow accumulation raster
        std::string sed_file;  ///< Input file name for the sediment layer
        bool fix_random_seed;  ///< Fix the random seed

    public:
        /// \brief Default Parameters object
        Parameters();

        /// \brief Create Parameters object with values loaded from file
        Parameters(const std::string& parameter_file);

        /// \brief Set parameter value
        /// \param U New parameter value
        void set_U(real_type U_);

        /// \brief Get parameter value
        /// \returns U Parameter value
        real_type get_U() const { return U; }

        /// \brief Set parameter value
        /// \param K New parameter value
        void set_K(real_type K_);

        /// \brief Get parameter value
        /// \returns K Parameter value
        real_type get_K() const { return K; }

        /// \brief Set parameter value
        /// \param D New parameter value
        void set_D(real_type D_);

        /// \brief Get parameter value
        /// \returns D Parameter value
        real_type get_D() const { return D; }

        /// \brief Set parameter value
        /// \param melt New parameter value
        void set_melt(real_type melt_);

        /// \brief Get parameter value
        /// \returns melt Parameter value
        real_type get_melt() const { return melt; }

        /// \brief Set parameter value
        /// \param timestep New parameter value
        void set_timestep(real_type timestep_);

        /// \brief Get parameter value
        /// \returns timestep Parameter value
        real_type get_timestep() const { return timestep; }

        /// \brief Get parameter value
        /// \returns ann_timestep Parameter value
        real_type get_ann_timestep() const { return ann_timestep; }

        /// \brief Set parameter value
        /// \param printinterval New parameter value
        void set_printinterval(int printinterval_);

        /// \brief Get parameter value
        /// \returns printinterval Parameter value
        int get_printinterval() const { return printinterval; }

        /// \brief Set parameter value
        /// \param thresholdarea New parameter value
        void set_thresholdarea(real_type thresholdarea_);

        /// \brief Get parameter value
        /// \returns thresholdarea Parameter value
        real_type get_thresholdarea() const { return thresholdarea; }

        /// \brief Set parameter value
        /// \param init_exposure_age New parameter value
        void set_init_exposure_age(real_type init_exposure_age_);

        /// \brief Get parameter value
        /// \returns init_exposure_age Parameter value
        real_type get_init_exposure_age() const { return init_exposure_age; }

        /// \brief Set parameter value
        /// \param init_sed_track New parameter value
        void set_init_sed_track(real_type init_sed_track_);

        /// \brief Get parameter value
        /// \returns init_sed_track Parameter value
        real_type get_init_sed_track() const { return init_sed_track; }

        /// \brief Set parameter value
        /// \param init_veg New parameter value
        void set_init_veg(real_type init_veg_);

        /// \brief Get parameter value
        /// \returns init_veg Parameter value
        real_type get_init_veg() const { return init_veg; }

        /// \brief Set parameter value
        /// \param year New parameter value
        void set_year(int year_);

        /// \brief Get parameter value
        /// \returns year Parameter value
        int get_year() const { return year; }

        /// \brief Set parameter value
        /// \param day New parameter value
        void set_day(int day_);

        /// \brief Get parameter value
        /// \returns day Parameter value
        int get_day() const { return day; }

        /// \brief Set parameter value
        /// \param hour New parameter value
        void set_hour(int hour_);

        /// \brief Get parameter value
        /// \returns hour Parameter value
        int get_hour() const { return hour; }

        /// \brief Set parameter value
        /// \param minute New parameter value
        void set_minute(int minute_);

        /// \brief Get parameter value
        /// \returns minute Parameter value
        int get_minute() const { return minute; }

        /// \brief Set parameter value
        /// \param end_year New parameter value
        void set_end_year(int end_year_);

        /// \brief Get parameter value
        /// \returns end_year Parameter value
        int get_end_year() const { return end_year; }

        /// \brief Set parameter value
        /// \param end_day New parameter value
        void set_end_day(int end_day_);

        /// \brief Get parameter value
        /// \returns end_day Parameter value
        int get_end_day() const { return end_day; }

        /// \brief Set parameter value
        /// \param lattitude New parameter value
        void set_lattitude(real_type lattitude_);

        /// \brief Get parameter value
        /// \returns lattitude Parameter value
        real_type get_lattitude() const { return lattitude; }

        /// \brief Set parameter value
        /// \param longitude New parameter value
        void set_longitude(real_type longitude_);

        /// \brief Get parameter value
        /// \returns longitude Parameter value
        real_type get_longitude() const { return longitude; }

        /// \brief Set parameter value
        /// \param stdmed New parameter value
        void set_stdmed(real_type stdmed_);

        /// \brief Get parameter value
        /// \returns stdmed Parameter value
        real_type get_stdmed() const { return stdmed; }

        /// \brief Set parameter value
        /// \param declination New parameter value
        void set_declination(real_type declination_);

        /// \brief Get parameter value
        /// \returns declination Parameter value
        real_type get_declination() const { return declination; }

        /// \brief Set parameter value
        /// \param altitude New parameter value
        void set_altitude(real_type altitude_);

        /// \brief Get parameter value
        /// \returns altitude Parameter value
        real_type get_altitude() const { return altitude; }

        /// \brief Set parameter value
        /// \param azimuth New parameter value
        void set_azimuth(real_type azimuth_);

        /// \brief Get parameter value
        /// \returns azimuth Parameter value
        real_type get_azimuth() const { return azimuth; }

        /// \brief Set parameter value
        /// \param fa_file_ New parameter value
        void set_fa_file(const std::string& fa_file_);

        /// \brief Get parameter value
        /// \returns fa_file Parameter value
        std::string get_fa_file() const { return fa_file; }

        /// \brief Set parameter value
        /// \param topo_file_ New parameter value
        void set_topo_file(const std::string& topo_file_);

        /// \brief Get parameter value
        /// \returns topo_file Parameter value
        std::string get_topo_file() const { return topo_file; }

        /// \brief Set parameter value
        /// \param sed_file_ New parameter value
        void set_sed_file(const std::string& sed_file_);

        /// \brief Get parameter value
        /// \returns sed_file Parameter value
        std::string get_sed_file() const { return sed_file; }

        /// \brief Set parameter value
        /// \param fix_random_seed_ New parameter value
        void set_fix_random_seed(bool fix_random_seed_);

        /// \brief Get parameter value
        /// \returns fix_random_seed Parameter value
        bool get_fix_random_seed() const { return fix_random_seed; }
};

#endif
