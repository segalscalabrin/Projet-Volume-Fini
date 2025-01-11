#ifndef EXACTE_HPP
#define EXACTE_HPP

#include "include.hpp"
#include "fonctions.hpp"

double solExacte(double x, double y, double t, int cas);

PiercedVector<double> computeSolExacte(Grid *grid, AGeometry *geo, Data *data, double t, int cas);

#endif
