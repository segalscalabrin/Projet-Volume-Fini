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
    // Initialiser les vecteurs auxiliaires
    std::vector<double> new_phi0(_grid->nbCells(), 0.0);
    std::vector<double> fp0(_grid->nbCells(), 0.0), fm0(_grid->nbCells(), 0.0);
    std::vector<double> gp0(_grid->nbCells(), 0.0), gm0(_grid->nbCells(), 0.0);

    // Convertir en PiercedVector
    PiercedVector<double> new_phi = VtoPV(new_phi0, _grid);
    PiercedVector<double> fp = VtoPV(fp0, _grid);
    PiercedVector<double> fm = VtoPV(fm0, _grid);
    PiercedVector<double> gp = VtoPV(gp0, _grid);
    PiercedVector<double> gm = VtoPV(gm0, _grid);

    // Regrouper les flux dans un std::array
    std::array<PiercedVector<double>, 4> flux = {fp, fm, gp, gm};

    computeFlux(flux);

    for (auto &cell : _grid->getCells()) {
        const long &id = cell.getId();
        double delta;
        std::vector<long> neighborsId;
        std::vector<double> neighborsPhi;

        delta = _data->h;

        new_phi[id] = _phi[id] - _data->dt/delta * (flux[1][id] - flux[0][id]) - _data->dt/delta * (flux[3][id] - flux[2][id]);
    }

    _phi = new_phi;
}



void TransportScheme::computeFlux(std::array<PiercedVector<double>, 4>& flux)
{
    std::array<long, 2> ownersId;
    NPoint center1, center2;
    auto& interfaces = _grid->getInterfaces();

    double f(0), g(0);

    for(auto inter : interfaces) {
        ownersId = inter.getOwnerNeigh();
        center1 = _grid->evalCellCentroid(ownersId[0]);
        if (!inter.isBorder()) {
            center2 = _grid->evalCellCentroid(ownersId[1]);
            if (center1[NPX] < center2[NPX] && center1[NPZ] == center2[NPZ]) {
                f = F_ordre1(_phi[ownersId[1]], _phi[ownersId[0]]);
                flux[1][ownersId[0]] = f;
                flux[0][ownersId[1]] = f;
            } 
            else if (center1[NPX] > center2[NPX] && center1[NPZ] == center2[NPZ]) {
                f = F_ordre1(_phi[ownersId[0]], _phi[ownersId[1]]);
                flux[0][ownersId[0]] = f;
                flux[1][ownersId[1]] = f;
            } 
            else if (center1[NPY] < center2[NPY] && center1[NPZ] == center2[NPZ]) {
                f = G_ordre1(_phi[ownersId[1]], _phi[ownersId[0]]);
                flux[3][ownersId[0]] = f;
                flux[2][ownersId[1]] = f;
            } 
            else if (center1[NPY] > center2[NPY] && center1[NPZ] == center2[NPZ]) {
                f = G_ordre1(_phi[ownersId[0]], _phi[ownersId[1]]);
                flux[2][ownersId[0]] = f;
                flux[3][ownersId[1]] = f;
            }
        }
        else {
            center2 = _grid->evalInterfaceCentroid(inter.getId());
            f = F_ordre1(_phi[ownersId[0]], _phi[ownersId[0]]);
            g = G_ordre1(_phi[ownersId[0]], _phi[ownersId[0]]);
            if (center1[NPX] < center2[NPX] && center1[NPZ] == center2[NPZ]) {
                flux[1][ownersId[0]] = f;
            } 
            else if (center1[NPX] > center2[NPX] && center1[NPZ] == center2[NPZ]) {
                flux[0][ownersId[0]] = f;
            } 
            else if (center1[NPY] < center2[NPY] && center1[NPZ] == center2[NPZ]) {
                flux[3][ownersId[0]] = g;
            } 
            else if (center1[NPY] > center2[NPY] && center1[NPZ] == center2[NPZ]) {
                flux[2][ownersId[0]] = g;
            }
        }
    }
}

/*
std::vector<long> TransportScheme::getCellNeighbour(long cellId)
{
    int k(0);
    std::vector<long> neighboursId(6, -1);
    long neighbourId(-1);
    NPoint centerTarget;
    NPoint centerNeighbor;

    centerTarget = _grid->evalCellCentroid(cellId);

    auto& interfaces(_grid->getInterfaces());







    std::cout << "x: " << centerTarget[0] << " y: " << centerTarget[1] << " z: " << centerTarget[2] << std::endl;
    std::cout << std::endl;

    for (auto inter : interfaces) {
        NPoint centerInter = _grid->evalInterfaceCentroid(inter.getId());
        std::cout << "x: " << centerInter[0] << " y: " << centerInter[1] << " z: " << centerInter[2] << std::endl;
        std::cout << inter.isBorder() << std::endl;
        std::cout << k << std::endl;      
        k++;
    }
    
    for (auto inter : interfaces) {
        if(!inter.isBorder()) {
            neighbourId = inter.getOwnerNeigh()[1];
            NPoint centerNeighbor = _grid->evalCellCentroid(neighbourId); 

            if (centerNeighbor[NPX] > centerTarget[NPX] && centerNeighbor[NPZ] == centerTarget[NPZ]) {
                neighboursId[1] = neighbourId;
            } 
            else if (centerNeighbor[NPX] < centerTarget[NPX] && centerNeighbor[NPZ] == centerTarget[NPZ]) {
                neighboursId[0] = neighbourId;
            } 
            else if (centerNeighbor[NPY] > centerTarget[NPY] && centerNeighbor[NPZ] == centerTarget[NPZ]) {
                neighboursId[3] = neighbourId;
            } 
            else if (centerNeighbor[NPY] < centerTarget[NPY] && centerNeighbor[NPZ] == centerTarget[NPZ]) {
                neighboursId[2] = neighbourId;
            }
        }
    }

    for (int i=0; i<4; i++) {
        std::cout << neighboursId[i] << std::endl;
        if(neighboursId[i] == -1) {
            if (i==0) {
                neighboursId[i] = neighboursId[1];
            }
            else if (i==1) {
                neighboursId[i] = neighboursId[0];
            }
            else if (i==2) {
                neighboursId[i] = neighboursId[3];
            }
            else if (i==3) {
                neighboursId[i] = neighboursId[2];
            }
        }
    }    
    for (int i=0; i<4; i++) {
        std::cout << neighboursId[i] << std::endl;
    }

    return neighboursId;
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
*/

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