#ifndef _RASTER_HPP_
#define _RASTER_HPP_

#include <vector>
#include <string>
#include "global_defs.h"

class Raster {
    private:
        int size_x, size_y;
        std::vector<real_type> data;
        real_type xllcorner, yllcorner;
        real_type deltax;
        real_type nodata;
        std::vector<int> idx;

    public:
        Raster();
        Raster(int size_x_, int size_y_);
        Raster(int size_x_, int size_y_, real_type value);
        Raster(const std::string &filename);
        const real_type& operator()(int i, int j) const;
        real_type& operator()(int i, int j);
        void resize(int size_x_, int size_y_);
        void load(const std::string &filename);
        void save(const std::string &filename);
        void set_data(real_type value);
        bool is_nodata(int i, int j);
        void set_nodata(int i, int j);
        void sort_data();
        void get_sorted_ij(int t, int &i, int &j);
        real_type get_size_x() { return size_x; }
        real_type get_size_y() { return size_y; }
        real_type get_xllcorner() { return xllcorner; }
        real_type get_yllcorner() { return yllcorner; }
        real_type get_deltax() { return deltax; }
        real_type get_nodata() { return nodata; }
};

#endif
