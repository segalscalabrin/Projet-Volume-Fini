#include "fonctions.hpp"


double vitesseX(double x, double y, int cas)
{
    switch (cas) {
        case 1:
            return 1.0;
        case 2:
            return 0.5;
        case 3:
            return -y;
        case 4:
            return 1.0;
        default:
            std::cout << "Erreur dans le numéro du cas" << std::endl;
            exit(0);
    }
}

double vitesseY(double x, double y, int cas)
{
    switch (cas) {
        case 1:
            return 1.0;
        case 2:
            return 0.5;
        case 3:
            return -y;
        case 4:
            return 0.0;
        default:
            std::cout << "Erreur dans le numéro du cas" << std::endl;
            exit(0);
    }
}

double solInit(double x, double y, int cas)
{
    switch (cas) {
        case 1:
            return exp(-x*x/0.0075 - y*y/0.0075);
        case 2:
            if (sqrt(x*x + y*y)<=0.5) {
                return std::cos(x*std::acos(-1)*2.) + std::cos(y*std::acos(-1)*2.);
            }
            else {
                return 0;
            }
        case 3:
            if (sqrt((x-0.25)*(x-0.25) + (y-0.25)*(y-0.25))<=0.4) {
                return 1;
            }
            else {
                return 0;
            }
        default:
            std::cout << "Erreur dans le numéro du cas" << std::endl;
            exit(0);
    }
}

PiercedVector<double> computeSolInit(Grid *grid, AGeometry *geo, int cas)
{
    PiercedVector<double> vectSolInit;
    if (cas < 4) {
        std::vector<double> solInit0(grid->nbCells(), 0.0);
        vectSolInit = VtoPV(solInit0, grid);
        for(auto cell : grid->getCells()) {
            long cellId = cell.getId();
            NPoint center;
            center = grid->evalCellCentroid(cellId);
            vectSolInit[cellId] = solInit(center[NPX], center[NPY], cas);
        }

        return vectSolInit;
    }
    else {
        geo->computeLevelSet(grid);
        vectSolInit = geo->getPhi();

        return vectSolInit;
    }
}