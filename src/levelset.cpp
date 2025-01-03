/* -------------------------------------------------------------------------*\
 *
 *  NEOS
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of Neos.
 *
 *  Neos is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  Neos is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Neos. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

/* sudo docker run -it -v $(pwd):/builds/workspace neos */

#include "Neos.hpp"

#include "main_loop.h"

#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <array>

#define TYPEGEOM 1

using namespace neos;

void initializeProgram(int &ac, char **av);
Grid *createGrid(int level, int dimension);
void finalizeProgram(Grid *g);

int main(int ac, char **av) {
    initializeProgram(ac, av);

    int level = 5;
    int dimension = 2;
    double final_time = 1.0;
    int time_steps = 5;
    double dt = final_time / time_steps;

    Grid *g = createGrid(level, dimension);
    timeLoop(g, final_time, time_steps, dt);

    finalizeProgram(g);
    return 0;
}


void initializeProgram(int &ac, char **av) {
    Neos_Init(&ac, &av);
}

Grid *createGrid(int level, int dimension) {
    Grid *g = new Grid(-1.0, -1.0, 0.0, 2.0, 2.0 / pow(2, level), dimension);
    std::cerr << "Initial cell count: " << g->getCellCount() << std::endl;
    return g;
}



void finalizeProgram(Grid *g) {
    delete g;
    Neos_Finalize();
}
