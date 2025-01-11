#include "transport.hpp"


void TransportScheme::initializePhi()
{    
    _phi = computeSolInit(_grid, _geo, _data->cas);
    _phiExact = computeSolExacte(_grid, _geo, _data, _t, _data->cas);
}

const PiercedVector<double>& TransportScheme::getPhi()
{
    return _phi;
}

const PiercedVector<double>& TransportScheme::getPhiExact()
{
    return _phiExact;
}

const double& TransportScheme::getT()
{
    return _t;
}


void TransportScheme::computePhi()
{
    // Initialiser les vecteurs auxiliaires
    std::vector<double> new_phi0(_grid->nbCells(), 0.0);
    PiercedVector<double> new_phi = VtoPV(new_phi0, _grid);
    _t += _data->dt;

    for (auto &cell : _grid->getCells()) {
        const long &cellId = cell.getId();
        double delta;
        
        std::array<double, 4> flux;
        int ordreMood(_data->ordre);

        flux = computeCellFlux(cellId, ordreMood);
        new_phi[cellId] = _phi[cellId]  - (_data->dt / _data->h) * (flux[1] - flux[0])
                                        - (_data->dt / _data->h) * (flux[3] - flux[2]);

        if (cellId == 0) {
            std::array<long, 4> neigh(getCellNeighs(cellId));
            std::cout << "Phi i: " << _phi[cellId] << std::endl;
            std::cout << "Phi imx: " << _phi[neigh[0]] << std::endl;
            std::cout << "Phi ipx: " << _phi[neigh[1]] << std::endl;
            std::cout << std::endl;
            std::cout << "Phi new: " << new_phi[cellId] << std::endl;

            std::cout << "Flux mx: " << flux[0] << std::endl;
            std::cout << "Flux px: " << flux[1] << std::endl;
            std::cout << "Flux my: " << flux[2] << std::endl;
            std::cout << "Flux py: " << flux[3] << std::endl;
            std::cout << std::endl;
            std::cout << std::endl;

        }

        while(!critereMood(cellId, new_phi[cellId]) && ordreMood > 1) {
            ordreMood--;
            flux = computeCellFlux(cellId, ordreMood);
            new_phi[cellId] = _phi[cellId]  - (_data->dt / _data->h) * (flux[1] - flux[0])
                                            - (_data->dt / _data->h) * (flux[3] - flux[2]);
        }
    }

    _phi = new_phi;
    _phiExact = computeSolExacte(_grid, _geo, _data, _t, _data->cas);
}


std::array<long, 4> TransportScheme::getCellNeighs(long cellId)
{
    std::vector<long> unorderedNeighs;
    std::array<long, 4> neighs{-1, -1, -1, -1};
    unorderedNeighs = _grid->findCellNeighs(cellId, 1, false);
    NPoint center1, center2;
    center1 = _grid->evalCellCentroid(cellId);

    for(auto neighId : unorderedNeighs) {
        
        center2 = _grid->evalCellCentroid(neighId);
        if (center1[NPX] > center2[NPX] && center1[NPZ] == center2[NPZ]) {
            neighs[0] = neighId;
        } 
        else if (center1[NPX] < center2[NPX] && center1[NPZ] == center2[NPZ]) {
            neighs[1] = neighId;
        } 
        else if (center1[NPY] > center2[NPY] && center1[NPZ] == center2[NPZ]) {
            neighs[2] = neighId;
        } 
        else if (center1[NPY] < center2[NPY] && center1[NPZ] == center2[NPZ]) {
            neighs[3] = neighId;
        }
    }

    return neighs;
}

bool TransportScheme::critereMood(long cellId, double new_u) {
    std::array<long, 4> neighsId;
    double umin(_phi[cellId]), umax(_phi[cellId]);

    neighsId = getCellNeighs(cellId);

    for (int i=0; i<4; i++) {
        if (neighsId[i] != -1) {
            umin = std::min(umin, _phi[neighsId[i]]);
            umax = std::max(umax, _phi[neighsId[i]]);
        }
    }

    return (new_u > umin && new_u < umax);
}

double TransportScheme::boundaryCondition(std::array<double, 3> coord)
{
    if (_data->cas < 4) {
        std::cout << solExacte(coord[0], coord[1], _t, _data->cas) << std::endl;
        return solExacte(coord[0], coord[1], _t, _data->cas);
    }
    else {
        return _geo->getLevelSet(coord);
    }
}

double TransportScheme::cellBorderCheck(long cellId, long cellToCheckId, std::string signe, std::string direction) {
    double cellSize(0), borderPhi(0);
    if(cellToCheckId != -1) {
        borderPhi = _phi[cellToCheckId];
        return borderPhi;
    }
    else {
        std::array<double, 3> borderCoord(_grid->evalCellCentroid(cellId));
        cellSize = _grid->evalCellSize(cellId);
        if (direction==std::string("horizontal")) {
            if (signe==std::string("moins")) {
                borderCoord[0] -= cellSize;
                borderPhi = boundaryCondition(borderCoord);
                return borderPhi;
            }
            else {
                borderCoord[0] += cellSize;
                borderPhi = boundaryCondition(borderCoord);
                return borderPhi;  
            }
        }
        else if (direction==std::string("vertical")) {
            if (signe==std::string("moins")) {
                borderCoord[1] -= cellSize;
                borderPhi = boundaryCondition(borderCoord);
                return borderPhi;
            }
            else {
                borderCoord[1] += cellSize;
                borderPhi = boundaryCondition(borderCoord);
                return borderPhi;  
            }
        }
        else {
            std::cerr << "Erreur border check" << std::endl;
            exit(1);
        }
    }
}

std::array<double, 4> TransportScheme::computeCellFlux(int cellId, int ordre)
{
    std::array<long, 4> neighsId;
    std::array<double, 4> flux;

    neighsId = getCellNeighs(cellId);

    flux[0] = computeFlux(neighsId[0], cellId, neighsId[1], ordre, std::string("moins"), std::string("horizontal"));
    flux[1] = computeFlux(neighsId[0], cellId, neighsId[1], ordre, std::string("plus"), std::string("horizontal"));
    flux[2] = computeFlux(neighsId[2], cellId, neighsId[3], ordre, std::string("moins"), std::string("vertical"));
    flux[3] = computeFlux(neighsId[2], cellId, neighsId[3], ordre, std::string("plus"), std::string("vertical"));

    return flux;
}


double TransportScheme::computeFlux(long cellmId, long cellId, long cellpId, int ordre, std::string signe, std::string direction)
{   
    double phiM, phi, phiP, vx, vy;

    phi = _phi[cellId];
    phiM = cellBorderCheck(cellId, cellmId, std::string("moins"), direction);
    phiP = cellBorderCheck(cellId, cellpId, std::string("plus"), direction);

    std::array<double, 3> cellCoord(_grid->evalCellCentroid(cellId));
    vx = vitesseX(cellCoord[0], cellCoord[1], _data->cas);
    vy = vitesseY(cellCoord[0], cellCoord[1], _data->cas);

    if (direction==std::string("horizontal")) {
        return Flux_F(phiM, phi, phiP, vx, vy, ordre, signe);
    }
    else if (direction==std::string("vertical")) {
        return Flux_G(phiM, phi, phiP, vx, vy, ordre, signe);
    }
    else {
        std::cerr << "Direction flux non conforme" << std::endl;
        exit(1);
    }
}

double TransportScheme::Flux_F(double um, double u, double up, double vx, double vy, int ordre, std::string signe)
{   
    if (ordre == 3) {
        return F_ordre3(um, u, up, vx, vy, signe);
    }
    else if(ordre == 2) {
        if (signe==std::string("plus")) {
            return F_ordre2(u, up, vx, vy);
        }
        else if (signe==std::string("moins")) {
            return F_ordre2(um, u, vx, vy);
        }
        else {
            std::cerr << "Signe flux non conforme" << std::endl;
            exit(1);
        }
    }
    else if(ordre == 1) {
        if (signe==std::string("plus")) {
            return F_ordre1(u, up, vx, vy);
        }
        else if (signe==std::string("moins")) {
            return F_ordre1(um, u, vx, vy);
        }
        else {
            std::cerr << "Signe flux non conforme" << std::endl;
            exit(1);
        }
    }
    else {
        std::cerr << "Ordre schema non conforme" << std::endl;
        exit(1);
    }
}

double TransportScheme::Flux_G(double um, double u, double up, double vx, double vy, int ordre, std::string signe)
{   
    if (ordre == 3) {
        return G_ordre3(um, u, up, vx, vy, signe);
    }
    else if(ordre == 2) {
        if (signe==std::string("plus")) {
            return G_ordre2(u, up, vx, vy);
        }
        else if (signe==std::string("moins")) {
            return G_ordre2(um, u, vx, vy);
        }
        else {
            std::cerr << "Signe flux non conforme" << std::endl;
            exit(1);
        }
    }
    else if(ordre == 1) {
        if (signe==std::string("plus")) {
            return G_ordre1(u, up, vx, vy);
        }
        else if (signe==std::string("moins")) {
            return G_ordre1(um, u, vx, vy);
        }
        else {
            std::cerr << "Signe flux non conforme" << std::endl;
            exit(1);
        }
    }
    else {
        std::cerr << "Ordre schema non conforme" << std::endl;
        exit(1);
    }
}


double TransportScheme::F_ordre3(double um, double u, double up, double vx, double vy, std::string signe)
{   
    double ratio1(1./8.), ratio2(3./8.); 
    if (signe==std::string("plus")) {
        return vx*(ratio2*u + ratio2*up - ratio1*um);
    }
    else {
        return vx*(ratio2*u + ratio2*um - ratio1*up);
    }
}

double TransportScheme::G_ordre3(double um, double u, double up, double vx, double vy, std::string signe)
{   
    double ratio1(1./8.), ratio2(3./8.); 
    if (signe==std::string("plus")) {
        return vy*(ratio2*u + ratio2*up - ratio1*um);
    }
    else {
        return vy*(ratio2*u + ratio2*um - ratio1*up);
    }
}


double TransportScheme::F_ordre2(double um, double up, double vx, double vy)
{   
    double c = vx * _data->dt / _data->h;

    return 0.5 * vx * (up + um - c * (up - um));
}

double TransportScheme::G_ordre2(double um, double up, double vx, double vy)
{   
    double c = vy * _data->dt / _data->h; 

    return 0.5 * vy * (up + um - c * (up - um));
}

double TransportScheme::F_ordre1(double um, double up, double vx, double vy)
{   
    if (vx>=0) {
        return vx*um;
    }
    else {
        return vx*up;
    }
}

double TransportScheme::G_ordre1(double um, double up, double vx, double vy)
{   
    if (vy>=0) {
        return vy*um;
    }
    else {
        return vy*up;
    }
}