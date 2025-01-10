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
    Data            data;
    Grid            *grid = new Grid(-1.0, -1.0, 0.0, 2.0, 2.0 / std::pow(2, data.level), data.dim);
    ASphere         *geo  = new ASphere(-0.5, 0, 0, 0.5, 2);
    TransportScheme *trpt = new TransportScheme(&data, grid, geo);

    // Initial pos
    trpt->initializePhi();
    exportResults(grid, trpt->getPhi(), data.solName, i);

    // Loop across the scheme
    for (i = 1; i <= data.nbSteps; i++) {
        trpt->computePhi();
        if (i%10==0) {
            exportResults(grid, trpt->getPhi(), data.solName, i/10);
            exportResults(grid, trpt->getPhiExact(), data.solNameExacte, i/10);
        }
    }
    
    delete grid;
    delete trpt;
    delete geo;
    
    Neos_Finalize();

    return 0;
}


