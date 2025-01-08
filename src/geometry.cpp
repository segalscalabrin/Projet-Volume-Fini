#include "geometry.hpp"

PiercedVector<double> initializePhiWithGeometry(Grid *grid) 
{
    ASphere *geo  = new ASphere(0, 0, 0, 1.5, 2);
    PiercedVector<double> phi;

    geo->computeLevelSet(grid);
    phi = geo->getPhi();

    delete geo;
    return phi;
}