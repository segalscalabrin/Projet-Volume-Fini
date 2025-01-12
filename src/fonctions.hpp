#ifndef FONCTIONS_HPP
#define FONCTIONS_HPP

#include "include.hpp"
double vitesseX(double x, double y);
double vitesseY(double x, double y);

double computeError(Grid *grid, PiercedVector<double> phi, PiercedVector<double> phiExact);

#endif
