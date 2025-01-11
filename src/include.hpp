#ifndef INCLUDE_HPP
#define INCLUDE_HPP

#include "Neos.hpp"

using namespace neos;

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm> 

struct Data {
    int dim = 2;
    int level = 3;
    int maxLevel = 7;

    int nbSteps = 1000;
    double dt = 0.001;

    int ordre = 1;
    double h;

    double l = 2.0;

    double xmin = -1.0;
    double ymin = -1.0;
    double xmax = 1.0;
    double ymax = 1.0;

    int cas;
    std::vector<double> u;
    std::string solName;
    std::string solNameExacte;
    std::vector<double> phi;

    // Constructor to initialize members with a given level
    Data(int c, int ord, int lev)
        : level(lev),
          ordre(ord),
          cas(c),
          h(2.0 / std::pow(2, lev)),
          solName("solutions/sol_numerique/levelset_"),
          solNameExacte("solutions/sol_exacte/levelset_exact_")
    {
        u.push_back(1.0);
        u.push_back(0.0);
        u.push_back(0.0);
    }
};

#endif
