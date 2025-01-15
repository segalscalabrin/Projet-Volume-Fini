#include "fonctions.hpp"

double computeError(Grid *grid, PiercedVector<double> phi, PiercedVector<double> phiExact)
{
    double error(0.0);
    for (int i=0; i<grid->nbCells(); i++) {
        error += (phi[i] - phiExact[i])*(phi[i] - phiExact[i]);
    }   

    return sqrt(error/grid->nbCells());
}