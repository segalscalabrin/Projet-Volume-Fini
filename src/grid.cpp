#include "grid.hpp"

void initGrid(Grid *&grid, int dim, int level, int maxLevel) 
{
    // Validation des paramètres
    if (level <= 0 || dim <= 0 || maxLevel < level) {
        throw std::invalid_argument("Invalid grid parameters.");
    }

    // Crée la grille avec les dimensions spécifiées
    grid = new Grid(-1.0, -1.0, 0.0, 2.0, 2.0 / std::pow(2, level), dim);
}


std::vector<double> getLevelSet(Grid *grid, ASphere *geo)
{
    std::vector<double> phi(grid->nbCells(), 0.0);
    for (auto & cell : grid->getCells()) {
        const long &id = cell.getId();
        phi[id] = geo->getLevelSet(id);
    }
    return phi;
}