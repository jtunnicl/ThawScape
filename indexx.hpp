
#ifndef _INDEXX_H_
#define _INDEXX_H_

#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>


template <typename T>
class Indexx {
    /*
     * Takes a 2d array, converts to a flattened array, orders by value and returns
     * i, j indices from the input array
     * Note: there is no bounds checking
     */
    private:
        std::vector<int> idx;
        std::vector<T> flat_array;  // instead of allocating a new vector every time update is called
        int nx, ny, nxy;

    public:
        Indexx() : nx(0), ny(0), nxy(0) {}
        Indexx(int nx_, int ny_) : nx(nx_), ny(ny_), nxy(nx_ * ny_), flat_array(nx_ * ny_), idx(nx_ * ny_) {}
        
        void update_array(std::vector< std::vector<T> > &array) {
            // put 2d input array into flat array
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    flat_array[i * ny + j] = array[i][j];
                }
            }
            
            // initialise index array
            std::iota(idx.begin(), idx.end(), 0);

            // sort indexes based on comparing values in flat_array
            std::sort(idx.begin(), idx.end(), [this](int i1, int i2) {return flat_array[i1] < flat_array[i2]; });
        }

        void get_ij(int t, int &i, int &j) {
            if (t >= nxy || t < 0) {
                std::cerr << "Warning: out of bounds in Indexx::get_ij()" << std::endl;
            }

            // return i, j indices in original array corresponding to t'th highest value
            int index = idx[t];
            i = index / ny;
            j = index % ny;
        }
};

#endif
