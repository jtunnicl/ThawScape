
#ifndef _INDEXX_H_
#define _INDEXX_H_

#include <vector>
#include <numeric>
#include <algorithm>


template <typename T>
class Indexx {
    /*
     * Takes a 2d array, converts to a flattened array, orders by value and returns
     * i, j indices from the input array
     */
    private:
        std::vector<int> idx;
        int nx;
        int ny;
        int nxy;

    public:
        Indexx(int nx_, int ny_, std::vector< std::vector<T> > &array) : nx(nx_), ny(ny_), nxy(nx * ny) {
            update(array);
        }
        
        void update(std::vector< std::vector<T> > &array) {
            std::vector<T> flat_array(nxy);
            idx.resize(nxy);

            // put 2d input array into flat array
            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    flat_array[i * ny + j] = array[i][j];
                }
            }
            
            // initialise index array
            std::iota(idx.begin(), idx.end(), 0);

            // sort indexes based on comparing values in flat_array
            std::sort(idx.begin(), idx.end(), [&flat_array](int i1, int i2) {return flat_array[i1] < flat_array[i2]; });
        }

        void get_ij(int t, int &i, int &j) {
            // return i, j indices in original array corresponding to t'th highest value
            int index = idx[t];
            i = index / ny;
            j = index % ny;
        }
};

#endif
