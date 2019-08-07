#include <iostream>
#include <vector>
#include "grid_neighbours.h"

GridNeighbours::GridNeighbours() {}

GridNeighbours::GridNeighbours(const int size_x, const int size_y) {
    setup(size_x, size_y);
}

void GridNeighbours::setup(const int size_x, const int size_y) {
	idown_ = std::vector<int>(size_x);
	iup_ = std::vector<int>(size_x);
	jup_ = std::vector<int>(size_y);
	jdown_ = std::vector<int>(size_y);

    #pragma omp parallel for
	for (int i = 0; i <= size_x - 1; i++)
	{
		idown_[i] = i - 1;
		iup_[i] = i + 1;
	}
	idown_[0] = 0;
	iup_[size_x - 1] = size_x - 1;

    #pragma omp parallel for
	for (int j = 0; j <= size_y - 1; j++)
	{
		jdown_[j] = j - 1;
		jup_[j] = j + 1;
	}
	jdown_[0] = 0;
	jup_[size_y - 1] = size_y - 1;
}

int GridNeighbours::iup(const int i) {
    return iup_[i];
}

int GridNeighbours::idown(const int i) {
    return idown_[i];
}

int GridNeighbours::jup(const int j) {
    return jup_[j];
}

int GridNeighbours::jdown(const int j) {
    return jdown_[j];
}
