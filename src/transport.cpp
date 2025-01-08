#include "transport.hpp"


void TransportScheme::initializePhi()
{
    _geo->computeLevelSet(_grid);
}

std::vector<double> TransportScheme::getPhi()  
{
    std::vector<double> phi(_grid->nbCells(), 0.0);

    for (auto & cell : _grid->getCells()) {
        const long &id = cell.getId();
        phi[id] = _geo->getLevelSet(id);
    }

    return phi;   
}

void TransportScheme::compute()
{
    trpt->compute(_phi, data->u, data->dt);
    geo->updateCenter(data.dt, data.u);
}

