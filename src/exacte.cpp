#include "exacte.hpp"

double solExacte(double x, double y, double t, int cas)
{
    double uExacte;
    double Vx(vitesseX(x, y, cas)), Vy(vitesseY(x, y, cas));

    uExacte = solInit(x - Vx * t, y - Vy * t, cas);

    return uExacte;
}

PiercedVector<double> computeSolExacte(Grid *grid, AGeometry *geo, Data *data, double t, int cas)
{
    PiercedVector<double> vectSolExacte;
    if (cas < 4) {
        std::vector<double> solExacte0(grid->nbCells(), 0.0);
        vectSolExacte = VtoPV(solExacte0, grid);
        for(auto cell : grid->getCells()) {
            long cellId = cell.getId();
            NPoint center;
            center = grid->evalCellCentroid(cellId);
            vectSolExacte[cellId] = solExacte(center[NPX], center[NPY], t, cas);
        }

        return vectSolExacte;
    }
    else {
        geo->updateCenter(data->dt, data->u);
        geo->computeLevelSet(grid);
        vectSolExacte = geo->getPhi();

        return vectSolExacte;
    }

}

