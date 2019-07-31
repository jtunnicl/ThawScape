#ifndef _GRID_NEIGHBOURS_H_
#define _GRID_NEIGHBOURS_H_

#include <vector>
#include "global_defs.h"

class GridNeighbours {
    private:
        std::vector<int> iup_;
        std::vector<int> idown_;
        std::vector<int> jup_;
        std::vector<int> jdown_;

    public:
        GridNeighbours();

        void setup(int size_x, int size_y);

        int iup(int i);

        int idown(int i);
        
        int jup(int j);
        
        int jdown(int j);
};

#endif
