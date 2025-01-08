#include "transport.hpp"


void TransportScheme::initializePhi(PiercedVector<double> phi)
{
    _phi = phi;
}

const PiercedVector<double>& TransportScheme::getPhi()
{
    return _phi;
}

void TransportScheme::computePhi()
{
    std::vector<double> new_phi0(_grid->nbCells(), 0.0);
    PiercedVector<double> new_phi = VtoPV(new_phi0, _grid);

    for (auto &cell : _grid->getCells()) {
        const long &id = cell.getId();
        if (_grid->getCell(id).isInterior()) {
            std::vector<double> flux, delta;
            std::vector<long> neighborsId;
            std::vector<double> neighborsPhi;

            neighborsId = getCellNeighbour(id);
            neighborsPhi = getPhiNeighbour(neighborsId);

            delta = getDelta(id, neighborsId);
            flux = getFlux(id, neighborsPhi);

            new_phi[id] = _phi[id] - _data->dt/delta[0] * (flux[0] - flux[1]) - _data->dt/delta[0] * (flux[2] - flux[3]);
        }
    }

    _phi = new_phi;
}

std::vector<long> TransportScheme::getCellNeighbour(long cellId)
{
    std::vector<long> neighborsIdUnordered, neighborsId(4);

    /*
    
        A CHANGER CAR PAS POSSIBLE DE TROUVER LES PTN DE VOISINS DE SES MORTS

    */

    _grid->findCellNeighs(cellId);

    NPoint centerTarget = _grid->evalCellCentroid(cellId); 

    for (long neighborId : neighborsIdUnordered) {
        NPoint centerNeighbor = _grid->evalCellCentroid(neighborId); 

        if (centerNeighbor[NPX] > centerTarget[NPX] && centerNeighbor[NPZ] == centerTarget[NPZ]) {
            neighborsId[1] = neighborId;
        } 
        else if (centerNeighbor[NPX] < centerTarget[NPX] && centerNeighbor[NPZ] == centerTarget[NPZ]) {
            neighborsId[0] = neighborId;
        } 
        else if (centerNeighbor[NPY] > centerTarget[NPY] && centerNeighbor[NPZ] == centerTarget[NPZ]) {
            neighborsId[3] = neighborId;
        } 
        else if (centerNeighbor[NPY] < centerTarget[NPY] && centerNeighbor[NPZ] == centerTarget[NPZ]) {
            neighborsId[2] = neighborId;
        }
    }

    return neighborsId;
}

std::vector<double> TransportScheme::getPhiNeighbour(std::vector<long> neighborsId)
{
    std::vector<double> neighborsPhi;
    for(auto neighId : neighborsId) {
        neighborsPhi.push_back(_phi[neighId]);
    }

    return neighborsPhi;
}


std::vector<double> TransportScheme::getFlux(long cellId, std::vector<double> neighborsPhi)
{
    std::vector<double> flux;   

    flux.push_back(F_ordre1(neighborsPhi[1], _phi[cellId]));
    flux.push_back(F_ordre1(_phi[cellId], neighborsPhi[0]));
    flux.push_back(F_ordre1(neighborsPhi[3], _phi[cellId]));
    flux.push_back(F_ordre1(_phi[cellId], neighborsPhi[2]));

    return flux;
}


std::vector<double> TransportScheme::getDelta(long cellId, std::vector<long> neighborsId)
{
    std::vector<double> delta;
    std::array<double, 3> ic, xc, yc;
    ic = _grid->evalCellCentroid(cellId);
    xc = _grid->evalCellCentroid(neighborsId[0]);
    yc = _grid->evalCellCentroid(neighborsId[2]);

    delta.push_back(abs(ic[0]-xc[0]));
    delta.push_back(abs(ic[1]-yc[1]));

    return delta;
}


double TransportScheme::F_ordre1(double up, double um)
{   
    double ux(_data->u[0]);
    if (ux>=0) {
        return ux*um;
    }
    else {
        return ux*up;
    }
}

double TransportScheme::G_ordre1(double up, double um)
{   
    double uy(_data->u[1]);
    if (uy>=0) {
        return uy*um;
    }
    else {
        return uy*up;
    }
}