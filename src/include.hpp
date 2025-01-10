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
    int level = 6;
    int maxLevel = 7;

    int nbSteps = 1000;
    double dt = 0.001;

    int ordre = 2;
    double h;
    std::vector<double> u;
    std::string solName;
    std::string solNameExacte;
    std::vector<double> phi;

    // Constructor to initialize members
    Data()
        : h(2.0 / std::pow(2, level)), solName("solutions/sol_numerique/levelset_"), solNameExacte("solutions/sol_exacte/levelset_") {
        u.push_back(1.0);
        u.push_back(0.0);
        u.push_back(0.0);
    }
};

#endif
