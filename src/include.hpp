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

    int nbSteps = 300;
    double dt = 0.1;
    double h;
    std::vector<double> u;
    std::string solName;
    std::vector<double> phi;

    // Constructor to initialize members
    Data()
        : h(2.0 / std::pow(2, level)), solName("solutions/levelset_") {
        u.push_back(0.1);
        u.push_back(0.1);
        u.push_back(0.0);
    }
};

#endif
