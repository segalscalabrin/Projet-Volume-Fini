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
/* -------------------------------------------------------------------------*\
sudo docker run -it -v $(pwd):/builds/workspace neos
\*---------------------------------------------------------------------------*/


#include "include.hpp"

#include "export.hpp"
#include "transport.hpp"


int main(int argc, char **argv) {
    Neos_Init(&argc, &argv);
    int i(0);

    // Create struct and class pointer
    Data            data(atoi(argv[argc - 2]), atoi(argv[argc - 1]));
    Grid            *grid = new Grid(data.xmin, data.ymin, 0.0, data.l, data.l / std::pow(2, data.level), data.dim);
    ASphere         *geo  = new ASphere(0.0, 0.0, 0, 0.25, 2);
    TransportScheme *trpt = new TransportScheme(&data, grid, geo);

    // Initial pos
    trpt->initializePhi();
    exportResults(grid, trpt->getPhi(), data.solName, i);
    exportResults(grid, trpt->getPhiExact(), data.solNameExacte, i);

    // Loop across the scheme
    while (trpt->getT()<data.tmax) {
        i++;
        trpt->computePhi();
        exportResults(grid, trpt->getPhi(), data.solName, i);
        exportResults(grid, trpt->getPhiExact(), data.solNameExacte, i);
    }

    std::cout << "Erreur: " << computeError(grid, trpt->getPhi(), trpt->getPhiExact()) << std::endl;
    std::cout << "Dt: " << trpt->getDT() << std::endl;
    
    delete grid;
    delete geo;
    delete trpt;
    
    Neos_Finalize();

    return 0;
}


