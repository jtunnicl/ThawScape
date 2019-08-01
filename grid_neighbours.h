#ifndef _GRID_NEIGHBOURS_H_
#define _GRID_NEIGHBOURS_H_

#include <vector>


/// \brief Indexing of neighbouring cells in the grid
class GridNeighbours {
    private:
        std::vector<int> iup_;  ///< Indexes of up neighbours in x direction
        std::vector<int> idown_;  ///< Indexes of down neighbours in x direction
        std::vector<int> jup_;  ///< Indexes of up neighbours in y direction
        std::vector<int> jdown_;  ///< Indexes of down neighbours in y direction

    public:
        /// \brief Create empty GridNeighbours object
        GridNeighbours();

        /// \brief Initialise GridNeighbours for a given size grid
        /// \param size_x Number of cells in x direction
        /// \param size_y Number of cells in y direction
        void setup(int size_x, int size_y);

        /// \brief Get index of up neighbour in x direction for given cell i
        /// \param i Index of cell to get the neighbour of
        int iup(int i);

        /// \brief Get index of down neighbour in x direction for given cell i
        /// \param i Index of cell to get the neighbour of
        int idown(int i);
        
        /// \brief Get index of up neighbour in y direction for given cell j
        /// \param j Index of cell to get the neighbour of
        int jup(int j);
        
        /// \brief Get index of down neighbour in y direction for given cell j
        /// \param j Index of cell to get the neighbour of
        int jdown(int j);
};

#endif
