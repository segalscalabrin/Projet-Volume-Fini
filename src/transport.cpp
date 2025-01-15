#include "transport.hpp"


void TransportScheme::initializeScheme()
{    
    _geo->computeLevelSet(_grid);
    _phi = _geo->getPhi();
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
    // Initialiser les vecteurs auxiliaires
    std::vector<double> new_phi0(_grid->nbCells(), 0.0);
    PiercedVector<double> new_phi = VtoPV(new_phi0, _grid);
    _t += _dt;

    for (auto &cell : _grid->getCells()) {
        const long &cellId = cell.getId();
        int ordreMood(_data->ordre);

        std::array<double, 8> borderPhi;
        borderPhi = computeNeighValue(cellId);

        new_phi[cellId] = computeNewPhi(_phi[cellId], borderPhi, ordreMood);

        while(!critereMood(cellId, new_phi[cellId]) && ordreMood > 1) {  
            ordreMood--;
            new_phi[cellId] = computeNewPhi(_phi[cellId], borderPhi, ordreMood);
        }
    }

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


double TransportScheme::computeNewPhi(double phi, std::array<double, 8> borderPhi, int ordre)
{
    double vx(_data->u[0]), vy(_data->u[1]);
    double cx( _dt * vx / _grid->evalCellSize(0)), cy( _dt * vy / _grid->evalCellSize(0));
    if (ordre == 3) {
        return  phi 
                + phi * (- cx - cx*cx + cx*cx*cx/2.) 
                + borderPhi[1] * (cx*cx/2. + cx*cx*cx/6.) 
                + borderPhi[0] * (cx + cx*cx/2. + cx*cx*cx/3.) 
                + borderPhi[4] * (cx*cx*cx/6.) 
                + phi * (- cy - cy*cy + cy*cy*cy/2.) 
                + borderPhi[3] * (cy*cy/2. + cy*cy*cy/6.) 
                + borderPhi[2] * (cy + cy*cy/2. + cy*cy*cy/3.) 
                + borderPhi[6] * (cy*cy*cy/6.);
    }
    else if (ordre == 2) {
        return phi
                - (cx / 2.) * (borderPhi[1] - borderPhi[0]) + (cx * cx / 2.) * (borderPhi[1] - 2.*phi + borderPhi[0])
                - (cy / 2.) * (borderPhi[3] - borderPhi[2]) + (cy * cy / 2.) * (borderPhi[3] - 2.*phi + borderPhi[2]);
    }
    else {
        return phi
                - std::max(cx, 0.0) * (phi - borderPhi[0]) - std::min(cx, 0.0) * (borderPhi[1] - phi)
                - std::max(cy, 0.0) * (phi - borderPhi[2]) - std::min(cy, 0.0) * (borderPhi[3] - phi);
    }
}


std::array<double, 8> TransportScheme::computeNeighValue(int cellId)
{
    std::array<long, 8> neighsId;
    std::array<double, 8> neighValue;

    neighsId = getCellNeighs(cellId);

    neighValue[0] = cellBorderCheck(cellId, neighsId[0], std::string("moins"), std::string("horizontal"));
    neighValue[1] = cellBorderCheck(cellId, neighsId[1], std::string("plus"), std::string("horizontal"));
    neighValue[2] = cellBorderCheck(cellId, neighsId[2], std::string("moins"), std::string("vertical"));
    neighValue[3] = cellBorderCheck(cellId, neighsId[3], std::string("plus"), std::string("vertical"));
    neighValue[4] = cellBorderCheck(cellId, neighsId[4], std::string("moinsmoins"), std::string("horizontal"));
    neighValue[5] = cellBorderCheck(cellId, neighsId[5], std::string("plusplus"), std::string("horizontal"));
    neighValue[6] = cellBorderCheck(cellId, neighsId[6], std::string("moinsmoins"), std::string("vertical"));
    neighValue[7] = cellBorderCheck(cellId, neighsId[7], std::string("plusplus"), std::string("vertical"));

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
    NPoint centerCell, centerBorder;
    centerCell = _grid->evalCellCentroid(cellId);

    for(auto neighId : unorderedNeighs) {
        
        centerBorder = _grid->evalCellCentroid(neighId);
        if (centerCell[NPX] > centerBorder[NPX] && centerCell[NPY] == centerBorder[NPY] && centerCell[NPZ] == centerBorder[NPZ]) {
            neighs[0] = neighId;
        } 
        else if (centerCell[NPX] < centerBorder[NPX] && centerCell[NPY] == centerBorder[NPY] && centerCell[NPZ] == centerBorder[NPZ]) {
            neighs[1] = neighId;
        } 
        else if (centerCell[NPX] == centerBorder[NPX] && centerCell[NPY] > centerBorder[NPY] && centerCell[NPZ] == centerBorder[NPZ]) {
            neighs[2] = neighId;
        } 
        else if (centerCell[NPX] == centerBorder[NPX] && centerCell[NPY] < centerBorder[NPY] && centerCell[NPZ] == centerBorder[NPZ]) {
            neighs[3] = neighId;
        }
    }

    return neighs;
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