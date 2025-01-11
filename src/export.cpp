#include "export.hpp"

void exportResults(Grid *grid, const PiercedVector<double>& levelSetValues, std::string fileName, int iteration) 
{
    std::string file;
    file = fileName + std::to_string(iteration);
    grid->setExportName(file);

    std::vector<double> LSV = PVtoV(levelSetValues, grid);
    grid->addData("LevelSet", LSV);

    grid->write();
}


double computeError(Grid *grid, PiercedVector<double> phi, PiercedVector<double> phiExact)
{
    double error(0.0);
    for (int i=0; i<phi.size(); i++) {
        error += (phi[i] - phiExact[i])*(phi[i] - phiExact[i]);
    }   

    return sqrt(error/phi.size());
}