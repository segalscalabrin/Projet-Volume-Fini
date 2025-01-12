#include "fonctions.hpp"

double vitesseX(double x, double y)
{
    return 1.0;
}

double vitesseY(double x, double y)
{
    return 1.0;
}


double computeError(Grid *grid, PiercedVector<double> phi, PiercedVector<double> phiExact)
{
    double error(0.0);
    for (int i=0; i<phi.size(); i++) {
        error += (phi[i] - phiExact[i])*(phi[i] - phiExact[i]);
    }   

    return sqrt(error/phi.size());
}