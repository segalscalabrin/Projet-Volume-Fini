#ifndef INCLUDE_HPP
#define INCLUDE_HPP

#include "Neos.hpp"

using namespace neos;

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>

struct Data {
    int dim = 2;
    int level = 3;

    int ordre = 1;

    double l = 1.0;

    double xmin = 0.0;
    double ymin = 0.0;
    double xmax = 1.0;
    double ymax = 1.0;

    double tmax = 1.0;

    std::vector<double> u;
    std::string solName;
    std::string solNameExacte;
    std::vector<double> phi;

    // Constructor to initialize members with a given level
    Data(int ord, int lev)
        : level(lev),
          ordre(ord),
          solName("solutions/sol_numerique/levelset_"),
          solNameExacte("solutions/sol_exacte/levelset_exact_")
    {
        u.push_back(1.0);
        u.push_back(0.0);
        u.push_back(0.0);
    }
};

#endif
