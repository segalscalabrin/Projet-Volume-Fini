#ifndef GRID_HPP
#define GRID_HPP

#include "include.hpp"

void initGrid(Grid *&grid, int dim, int level, int maxLevel);

std::vector<double> getLevelSet(Grid *grid, ASphere *geo);

#endif
