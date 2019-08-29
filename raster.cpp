#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include "global_defs.h"
#include "utility.h"
#include "raster.h"



// create empty Raster
Raster::Raster() : size_x(0), size_y(0), data(), xllcorner(0), yllcorner(0), deltax(1), nodata(-99999), save_prec(-1) {}

// create Raster of given size with no data
Raster::Raster(int size_x_, int size_y_) : Raster() {
    resize(size_x_, size_y_);
}

// create Raster of given size with data initialised to specified value
Raster::Raster(int size_x_, int size_y_, real_type value) : Raster(size_x_, size_y_) {
    set_data(value);
}

// create Raster by loading from file
Raster::Raster(const std::string &filename) : Raster() {
    load(filename);
}

// operators for accessing underlying data
const real_type& Raster::operator()(int i, int j) const {
#ifdef NDEBUG
    return data[i * size_y + j];
#else
    return data.at(i * size_y + j);
#endif
}

real_type& Raster::operator()(int i, int j) {
#ifdef NDEBUG
    return data[i * size_y + j];
#else
    return data.at(i * size_y + j);
#endif
}

// destructively resize the data
void Raster::resize(int size_x_, int size_y_) {
    if (size_x_ != size_x || size_y_ != size_y) {
        size_x = size_x_;
        size_y = size_y_;
        data = std::vector<real_type>(size_x * size_y);
        idx = std::vector<int>();
    }
}

// load Raster from file
void Raster::load(const std::string &filename) {
    std::ifstream fin(filename);

    if (!fin) {
        Util::Error("Well that didn't work ..!  Missing or invalid file: " + std::string(filename), 1);
    }
    else {
        Util::Warning("Reading raster without any checks or guarantees ...");
    }

    // read 6 lines of metadata
    std::string key;
    fin >> key; fin >> size_y; // ncols //NOTE: Pelltier's code was originally written for [x][y] indexing; Saga uses [y][x].
    fin >> key; fin >> size_x; // nrows
    fin >> key; fin >> xllcorner;
    fin >> key; fin >> yllcorner;
    fin >> key; fin >> deltax;
    fin >> key; fin >> nodata;

    // create vector
    data = std::vector<real_type>(size_x * size_y);

    // read data
    for (int x = 0; x < size_x; x++)
    {
        for (int y = 0; y < size_y; y++)
        {
            fin >> this->operator()(x, y);
        }
    }

    Util::Info("Done reading raster");
}

// save Raster to file
void Raster::save(const std::string &filename) {
    std::ofstream fout(filename);
    if (!fout) {
        Util::Error("Error opening file to save raster: " + filename, 1);
    }

    // write arcgrid format
    fout << "ncols " << size_y << std::endl;
    fout << "nrows " << size_x << std::endl;
    fout << "xllcorner " << xllcorner << std::endl;
    fout << "yllcorner " << yllcorner << std::endl;
    fout << "cellsize " << deltax << std::endl;
    fout << "NODATA_value " << nodata << std::endl;

    // setting precision of output
    if (save_prec > 0) {
        fout << std::fixed << std::setprecision(save_prec);
    }

    // write the data
    for (int i = 0; i < size_x; i++)
    {
        for (int j = 0; j < size_y; j++)
        {
            fout << this->operator()(i, j) << " ";
        }
        fout << std::endl;
    }
    fout.close();
}

// set all elements of data to a value
void Raster::set_data(real_type value) {
    std::fill(data.begin(), data.end(), value);
}

// check if an element is set to nodata
bool Raster::is_nodata(int i, int j) {
    return (this->operator()(i, j) == nodata);
}

// set a pixel to nodata
void Raster::set_pixel_nodata(int i, int j) {
    this->operator()(i, j) = nodata;
}

// make a list of indices sorted by data values
void Raster::sort_data() {
    idx = std::vector<int>(size_x * size_y);

    // initialise index array
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values (lowest to highest) in data
    std::sort(idx.begin(), idx.end(), [this](int i1, int i2) {return data[i1] < data[i2]; });
}

// return specified element from the ordered data values
void Raster::get_sorted_ij(int t, int &i, int &j) {
    if (t >= idx.size() || t < 0) {
        std::cerr << "Warning: out of bounds in Raster::get_sorted_ij()" << std::endl;
        i = -1;
        j = -1;
    }
    else {
        // return i, j indices in original array corresponding to t'th value in ordered array
        int index = idx[t];
        i = index / size_y;
        j = index % size_y;
    }
}

void Raster::set_save_precision(int prec) {
    save_prec = prec;
}

void Raster::set_deltax(const real_type deltax_) {
    deltax = deltax_;
}
