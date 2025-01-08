#include "export.hpp"

void exportResults(Grid *grid, std::vector<double> levelSetValues, std::string fileName, int iteration) 
{
    std::string file;

    file = fileName + std::to_string(iteration);

    grid->setExportName(file);

    grid->addData("LevelSet", levelSetValues);

    grid->write();
}
