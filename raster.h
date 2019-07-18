#ifndef _RASTER_HPP_
#define _RASTER_HPP_

#include <vector>
#include <string>
#include "global_defs.h"

class Raster {
    private:
        int size_x, size_y;
        std::vector<calcs_t> data;
        calcs_t xllcorner, yllcorner;
        calcs_t deltax;
        calcs_t nodata;
        std::vector<int> idx;

    public:
        Raster();
        Raster(int size_x_, int size_y_);
        Raster(int size_x_, int size_y_, calcs_t value);
        Raster(const std::string &filename);
        const calcs_t& operator()(int i, int j) const { return data.at(i * size_y + j); }
        calcs_t& operator()(int i, int j) { return data.at(i * size_y + j); }
        void resize(int size_x_, int size_y_);
        void load(const std::string &filename);
        void save(const std::string &filename);
        void set_data(calcs_t value);
        bool is_nodata(int i, int j);
        void set_nodata(int i, int j);
        void sort_data();
        void get_sorted_ij(int t, int &i, int &j);
        calcs_t get_size_x() { return size_x; }
        calcs_t get_size_y() { return size_y; }
        calcs_t get_xllcorner() { return xllcorner; }
        calcs_t get_yllcorner() { return yllcorner; }
        calcs_t get_deltax() { return deltax; }
        calcs_t get_nodata() { return nodata; }
};

#endif
