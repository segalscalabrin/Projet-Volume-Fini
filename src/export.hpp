#ifndef EXPORT_HPP
#define EXPORT_HPP

#include "include.hpp"

void exportResults(Grid *grid, const PiercedVector<double>& levelSetValues, std::string fileName, int iteration);

double computeError(Grid *grid, PiercedVector<double> phi, PiercedVector<double> phiExact);

#endif
