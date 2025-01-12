#include "transport.hpp"


void TransportScheme::initializePhi()
{    
    _data->u[0] = vitesseX(0., 0.);
    _data->u[1] = vitesseY(0., 0.);
    _geo->computeLevelSet(_grid);
    _phi = _geo->getPhi();
    std::vector<double> _phiPrime0(_grid->nbCells(), 0.0);
    std::vector<double> _phiSecond0(_grid->nbCells(), 0.0);
    _phiPrime = VtoPV(_phiPrime0, _grid);
    _phiSecond = VtoPV(_phiSecond0, _grid);
    _phiExact = _geo->getPhi();
    _t = 0;
    _dt = 0.5*_grid->evalCellSize(0)/std::abs(std::max(_data->u[0], _data->u[1]));
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

const double& TransportScheme::getDT()
{
    return _dt;
}

void TransportScheme::computePhi()
{
    _t += _dt;
    // Initialiser les vecteurs auxiliaires
    std::vector<double> new_phi0(_grid->nbCells(), 0.0);
    PiercedVector<double> new_phi_ord3 = VtoPV(new_phi0, _grid);
    PiercedVector<double> new_phi_ord2 = VtoPV(new_phi0, _grid);
    PiercedVector<double> new_phi_ord1 = VtoPV(new_phi0, _grid);
    PiercedVector<double> new_phi = VtoPV(new_phi0, _grid);
    int compteurMood(0);

    if (_data->ordre == 3) {
        new_phi_ord3 = computePhiRK(3);
    }
    if (_data->ordre == 3 || _data->ordre == 2) {
        new_phi_ord2 = computePhiRK(2);
    }
    new_phi_ord1 = computePhiRK(1);

    if (_data->ordre == 3) {
        new_phi = new_phi_ord3;
    }
    else if (_data->ordre == 2) {
        new_phi = new_phi_ord2;
    }
    else {
        new_phi = new_phi_ord1;
    }
    /*
    if (_data->ordre > 1) {
        for (auto &cell : _grid->getCells()) {
            long cellId = cell.getId();
            if (!critereMood(cellId, new_phi[cellId]) && _data->ordre==3) {  
                new_phi[cellId] = new_phi_ord2[cellId];
            }
            if (!critereMood(cellId, new_phi[cellId])) {  
                new_phi[cellId] = new_phi_ord1[cellId];
                compteurMood++;
            }
        }
    }
    */

   

    _phi = new_phi;

    _geo->updateCenter(_dt, _data->u);
    _geo->computeLevelSet(_grid);
    _phiExact = _geo->getPhi();
    
}


bool TransportScheme::critereMood(long cellId, double new_u) {
    std::array<long, 4> neighsId;
    double umin(_phi[cellId]), umax(_phi[cellId]);

    neighsId = orderedNeigh(cellId);

    for (int i=0; i<4; i++) {
        if (neighsId[i] != -1) {
            umin = std::min(umin, _phi[neighsId[i]]);
            umax = std::max(umax, _phi[neighsId[i]]);
        }
    }

    return (new_u > umin && new_u < umax);
}


PiercedVector<double> TransportScheme::computePhiRK(int ordre)
{
    std::vector<double> new_phi0(_grid->nbCells(), 0.0);
    PiercedVector<double> new_phi = VtoPV(new_phi0, _grid);
    std::array<double, 8> borderPhi;
    std::array<double, 8> borderPhiPrime;
    double vx(_data->u[0]), vy(_data->u[1]);
    double cx( _dt * vx / _grid->evalCellSize(0)), cy( _dt * vy / _grid->evalCellSize(0));

    for (auto &cell : _grid->getCells()) {
        const long &cellId = cell.getId();
        borderPhi = computeNeighValue(cellId, 0);
        _phiPrime[cellId] = _phi[cellId] - cx * computeFlux(_phi[cellId], borderPhi, ordre);
    }
    for (auto &cell : _grid->getCells()) {
        const long &cellId = cell.getId();
        borderPhi = computeNeighValue(cellId, 0);
        borderPhiPrime = computeNeighValue(cellId, 1);
        new_phi[cellId] = _phi[cellId] - cx/2.*(computeFlux(_phiPrime[cellId], borderPhiPrime, ordre) + computeFlux(_phi[cellId], borderPhi, ordre));
    }
    /*
    for (auto &cell : _grid->getCells()) {
        const long &cellId = cell.getId();
        borderPhi = computeNeighValue(cellId, 2);
        new_phi[cellId] = (1./3.)*_phi[cellId] + (2./3.) * computeNewPhi(_phiSecond[cellId], borderPhi, ordre);
    }
*/
    return new_phi;
}


double TransportScheme::computeNewPhi(double phi, std::array<double, 8> borderPhi, int ordre)
{
    double vx(_data->u[0]), vy(_data->u[1]);
    double cx( _dt * vx / _grid->evalCellSize(0)), cy( _dt * vy / _grid->evalCellSize(0));
    if (ordre == 3) {
        return phi
                -(cx / 2.) * (borderPhi[1]  - borderPhi[0]) 
                +(cx * cx / 2.) * (borderPhi[1] - 2*phi + borderPhi[0]) 
                -(cx * cx * cx / 12.) * (borderPhi[5] - 2.*borderPhi[1] + 2.*borderPhi[0] - borderPhi[4])
                -(cy / 2.) * (borderPhi[3]  - borderPhi[2]) 
                +(cy * cy / 2.) * (borderPhi[3] - 2*phi + borderPhi[2]) 
                -(cy * cy * cy / 12.) * (borderPhi[7] - 2.*borderPhi[3] + 2.*borderPhi[2] - borderPhi[6]);
    }
    else if (ordre == 2) {
        return phi - cx * computeFlux(phi, borderPhi, ordre);
    }
    else {
        if (cx >= 0) {
            if (cy >= 0) {
                return phi 
                        - cx * (phi - borderPhi[0])  // gauche
                        - cy * (phi - borderPhi[2]); // bas
            } else {
                return phi 
                        - cx * (phi - borderPhi[0])  // gauche
                        - cy * (borderPhi[3] - phi); // haut
            }
        } else {
            if (cy >= 0) {
                return phi 
                        - cx * (borderPhi[1] - phi)  // droite
                        - cy * (phi - borderPhi[2]); // bas
            } else {
                return phi 
                        - cx * (borderPhi[1] - phi)  // droite
                        - cy * (borderPhi[3] - phi); // haut
            }
        }
    }
}

double TransportScheme::computeFlux(double phi, std::array<double, 8> borderPhi, int ordre)
{
    double vx(_data->u[0]), vy(_data->u[1]);
    double cx( _dt * vx / _grid->evalCellSize(0)), cy( _dt * vy / _grid->evalCellSize(0));
    if (ordre == 3) {
        return phi
                -(cx / 2.) * (borderPhi[1]  - borderPhi[0]) 
                +(cx * cx / 2.) * (borderPhi[1] - 2*phi + borderPhi[0]) 
                -(cx * cx * cx / 12.) * (borderPhi[5] - 2.*borderPhi[1] + 2.*borderPhi[0] - borderPhi[4])
                -(cy / 2.) * (borderPhi[3]  - borderPhi[2]) 
                +(cy * cy / 2.) * (borderPhi[3] - 2*phi + borderPhi[2]) 
                -(cy * cy * cy / 12.) * (borderPhi[7] - 2.*borderPhi[3] + 2.*borderPhi[2] - borderPhi[6]);
    }
    else if (ordre == 2) {
        return  - (vx / 2.) * (borderPhi[1] - borderPhi[0]) + (vx * cx / 2.) * (borderPhi[1] - 2*phi + borderPhi[0])
                - (vy / 2.) * (borderPhi[3] - borderPhi[2]) + (vy * cy / 2.) * (borderPhi[3] - 2*phi + borderPhi[2]);
    }
    else {
        if (cx >= 0) {
            if (cy >= 0) {
                return phi 
                        - vx * (phi - borderPhi[0])  // gauche
                        - vy * (phi - borderPhi[2]); // bas
            } else {
                return phi 
                        - vx * (phi - borderPhi[0])  // gauche
                        - vy * (borderPhi[3] - phi); // haut
            }
        } else {
            if (cy >= 0) {
                return phi 
                        - cx * (borderPhi[1] - phi)  // droite
                        - cy * (phi - borderPhi[2]); // bas
            } else {
                return phi 
                        - cx * (borderPhi[1] - phi)  // droite
                        - cy * (borderPhi[3] - phi); // haut
            }
        }
    }
}


std::array<double, 8> TransportScheme::computeNeighValue(int cellId, int RK)
{
    std::array<long, 8> neighsId;
    std::array<double, 8> neighValue;

    neighsId = getCellNeighs(cellId);

    neighValue[0] = cellBorderCheck(cellId, neighsId[0], RK, std::string("moins"), std::string("horizontal"));
    neighValue[1] = cellBorderCheck(cellId, neighsId[1], RK, std::string("plus"), std::string("horizontal"));
    neighValue[2] = cellBorderCheck(cellId, neighsId[2], RK, std::string("moins"), std::string("vertical"));
    neighValue[3] = cellBorderCheck(cellId, neighsId[3], RK, std::string("plus"), std::string("vertical"));
    neighValue[4] = cellBorderCheck(cellId, neighsId[4], RK, std::string("moinsmoins"), std::string("horizontal"));
    neighValue[5] = cellBorderCheck(cellId, neighsId[5], RK, std::string("plusplus"), std::string("horizontal"));
    neighValue[6] = cellBorderCheck(cellId, neighsId[6], RK, std::string("moinsmoins"), std::string("vertical"));
    neighValue[7] = cellBorderCheck(cellId, neighsId[7], RK, std::string("plusplus"), std::string("vertical"));

    return neighValue;
}



std::array<long, 8> TransportScheme::getCellNeighs(long cellId)
{
    std::array<long, 4> nearNeighs{-1, -1, -1, -1};
    std::array<long, 8> neighs{-1, -1, -1, -1, -1, -1, -1, -1};
    nearNeighs = orderedNeigh(cellId);

    for(int i=0; i<4; i++){
        neighs[i] = nearNeighs[i];
    }

    if (neighs[0] != -1) {
        neighs[4] = orderedNeigh(neighs[0])[0];
    }
    if (neighs[1] != -1) {
        neighs[5] = orderedNeigh(neighs[1])[1];
    }
    if (neighs[2] != -1) {
        neighs[6] = orderedNeigh(neighs[2])[2];
    }
    if (neighs[3] != -1) {
        neighs[7] = orderedNeigh(neighs[3])[3];
    }

    return neighs;
}

std::array<long, 4> TransportScheme::orderedNeigh(long cellId)
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



double TransportScheme::cellBorderCheck(long cellId, long cellToCheckId, int RK, std::string signe, std::string direction) {
    double cellSize(0), borderPhi(0);
    if(cellToCheckId != -1) {
        if (RK==1) {
            borderPhi = _phiPrime[cellToCheckId];
        }
        else if (RK==2) {
            borderPhi = _phiSecond[cellToCheckId];
        }
        else {
            borderPhi = _phi[cellToCheckId];
        }
        return borderPhi;
    }
    else {
        std::array<double, 3> borderCoord(_grid->evalCellCentroid(cellId));
        cellSize = _grid->evalCellSize(cellId);
        if (direction==std::string("horizontal")) {
            if (signe==std::string("moins")) {
                borderCoord[0] -= cellSize;
                borderPhi = _geo->getLevelSet(borderCoord);
                return borderPhi;
            }
            else if (signe==std::string("plus")){
                borderCoord[0] += cellSize;
                borderPhi = _geo->getLevelSet(borderCoord);
                return borderPhi;  
            }
            else if (signe==std::string("moinsmoins")){
                borderCoord[0] -= 2*cellSize;
                borderPhi = _geo->getLevelSet(borderCoord);
                return borderPhi;  
            }
            else if (signe==std::string("plusplus")){
                borderCoord[0] += 2*cellSize;
                borderPhi = _geo->getLevelSet(borderCoord);
                return borderPhi;  
            }
        }
        else if (direction==std::string("vertical")) {
            if (signe==std::string("moins")) {
                borderCoord[1] -= cellSize;
                borderPhi = _geo->getLevelSet(borderCoord);
                return borderPhi;
            }
            else if (signe==std::string("plus")){
                borderCoord[1] += cellSize;
                borderPhi = _geo->getLevelSet(borderCoord);
                return borderPhi;  
            }
            else if (signe==std::string("moinsmoins")){
                borderCoord[1] -= 2*cellSize;
                borderPhi = _geo->getLevelSet(borderCoord);
                return borderPhi;  
            }
            else if (signe==std::string("plusplus")){
                borderCoord[1] += 2*cellSize;
                borderPhi = _geo->getLevelSet(borderCoord);
                return borderPhi;  
            }
        }
        else {
            std::cerr << "Erreur border check" << std::endl;
            exit(1);
        }
    }
    return 0.0;
}