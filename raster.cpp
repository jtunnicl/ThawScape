#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include "global_defs.h"
#include "utility.h"
#include "raster.h"



// create empty Raster
Raster::Raster() : size_x(0), size_y(0), data(), xllcorner(0), yllcorner(0), deltax(1), nodata(-99999) {}

// create Raster of given size with no data
Raster::Raster(int size_x_, int size_y_) : Raster() {
    resize(size_x_, size_y_);
}

// create Raster of given size with data initialised to specified value
Raster::Raster(int size_x_, int size_y_, calcs_t value) : Raster(size_x_, size_y_) {
    set_data(value);
}

// load Raster from file
Raster::Raster(std::string filename) {
    load(filename);
}

// destructively resize the data
void Raster::resize(int size_x_, int size_y_) {
    size_x = size_x_;
    size_y = size_y_;
    data = std::vector<calcs_t>(size_x * size_y);
    idx = std::vector<int>();
}

void Raster::load(std::string filename) {
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
    data = std::vector<calcs_t>(size_x * size_y);

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
void Raster::save(std::string filename) {
    std::ofstream fout(filename);
    if (!fout) {
        Util::Error("Error opening file to save raster: " + filename, 1);
    }

    // setting precision of output
//    fout << std::fixed << std::setprecision(12);

    // write arcgrid format
    fout << "ncols " << size_y << std::endl;
    fout << "nrows " << size_x << std::endl;
    fout << "xllcorner " << xllcorner << std::endl;
    fout << "yllcorner " << yllcorner << std::endl;
    fout << "cellsize " << deltax << std::endl;
    fout << "NODATA_value " << nodata << std::endl;

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
void Raster::set_data(calcs_t value) {
    std::fill(data.begin(), data.end(), value);
}

// check if an element is set to nodata
bool Raster::is_nodata(int i, int j) {
    return (this->operator()(i, j) == nodata);
}

// set an element to nodata
void Raster::set_nodata(int i, int j) {
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
    }

    // return i, j indices in original array corresponding to t'th value in ordered array
    int index = idx[t];
    i = index / size_y;
    j = index % size_y;
}
