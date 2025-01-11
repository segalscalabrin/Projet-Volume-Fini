#ifndef FONCTIONS_HPP
#define FONCTIONS_HPP

#include "include.hpp"

double vitesseX(double x, double y, int cas);
double vitesseY(double x, double y, int cas);
double solInit(double x, double y, int cas);
PiercedVector<double> computeSolInit(Grid *grid, AGeometry *geo, int cas);

#endif
