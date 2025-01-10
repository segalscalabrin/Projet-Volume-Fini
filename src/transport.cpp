#include "transport.hpp"


void TransportScheme::initializePhi()
{    
    _geo->computeLevelSet(_grid);
    _phi = _geo->getPhi();
    _phiExact = _geo->getPhi();
}

const PiercedVector<double>& TransportScheme::getPhi()
{
    return _phi;
}

const PiercedVector<double>& TransportScheme::getPhiExact()
{
    return _phiExact;
}

void TransportScheme::computePhi()
{
    // Initialiser les vecteurs auxiliaires
    std::vector<double> new_phi0(_grid->nbCells(), 0.0);
    PiercedVector<double> new_phi = VtoPV(new_phi0, _grid);

    for (auto &cell : _grid->getCells()) {
        const long &cellId = cell.getId();
        double delta;
        
        std::array<double, 4> flux;
        int ordreMood(_data->ordre);

        flux = computeCellFlux(cellId, ordreMood);
        new_phi[cellId] = _phi[cellId] - (_data->dt / _data->h) * (flux[1] - flux[0])
                            - (_data->dt / _data->h) * (flux[3] - flux[2]);

        while(!critereMood(cellId, new_phi[cellId]) && ordreMood > 1) {  
            ordreMood--;
            flux = computeCellFlux(cellId, ordreMood);
            new_phi[cellId] = _phi[cellId] - (_data->dt / _data->h) * (flux[1] - flux[0])
                                - (_data->dt / _data->h) * (flux[3] - flux[2]);
        }
        
    }

    _phi = new_phi;

    _geo->updateCenter(_data->dt, _data->u);
    _geo->computeLevelSet(_grid);
    _phiExact = _geo->getPhi();
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

long cellBorderCheck(long cellId, long cellToCheckId) {
    if(cellToCheckId != -1) {
        return cellToCheckId;
    }
    else {
        return cellId;
    }
}

std::array<double, 4> TransportScheme::computeCellFlux(int cellId, int ordre)
{
    std::array<long, 4> neighsId;
    std::array<double, 4> flux;

    neighsId = getCellNeighs(cellId);

    flux[0] = computeFlux(cellBorderCheck(cellId, neighsId[0]), cellId, cellBorderCheck(cellId, neighsId[1]), ordre, std::string("moins"), std::string("horizontal"));
    flux[1] = computeFlux(cellBorderCheck(cellId, neighsId[0]), cellId, cellBorderCheck(cellId, neighsId[1]), ordre, std::string("plus"), std::string("horizontal"));

    flux[2] = computeFlux(cellBorderCheck(cellId, neighsId[2]), cellId, cellBorderCheck(cellId, neighsId[3]), ordre, std::string("moins"), std::string("vertical"));
    flux[3] = computeFlux(cellBorderCheck(cellId, neighsId[2]), cellId, cellBorderCheck(cellId, neighsId[3]), ordre, std::string("plus"), std::string("vertical"));

    return flux;
}


double TransportScheme::computeFlux(long cellmId, long cellId, long cellpId, int ordre, std::string signe, std::string direction)
{   
    if (direction==std::string("horizontal")) {
        return Flux_F(_phi[cellmId], _phi[cellId], _phi[cellpId], ordre, signe);
    }
    else if (direction==std::string("vertical")) {
        return Flux_G(_phi[cellmId], _phi[cellId], _phi[cellpId], ordre, signe);
    }
    else {
        std::cerr << "Direction flux non conforme" << std::endl;
        exit(1);
    }
}

double TransportScheme::Flux_F(double um, double u, double up, int ordre, std::string signe)
{   
    if (ordre == 3) {
        std::cerr << "pas implémenté" << std::endl;
        exit(1);
    }
    else if(ordre == 2) {
        if (signe==std::string("plus")) {
            return F_ordre2(u, up);
        }
        else if (signe==std::string("moins")) {
            return F_ordre2(um, u);
        }
        else {
            std::cerr << "Signe flux non conforme" << std::endl;
            exit(1);
        }
    }
    else if(ordre == 1) {
        if (signe==std::string("plus")) {
            return F_ordre1(u, up);
        }
        else if (signe==std::string("moins")) {
            return F_ordre1(um, u);
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

double TransportScheme::Flux_G(double um, double u, double up, int ordre, std::string signe)
{   
    if (ordre == 3) {
        std::cerr << "pas implémenté" << std::endl;
        exit(1);
    }
    else if(ordre == 2) {
        if (signe==std::string("plus")) {
            return G_ordre2(u, up);
        }
        else if (signe==std::string("moins")) {
            return G_ordre2(um, u);
        }
        else {
            std::cerr << "Signe flux non conforme" << std::endl;
            exit(1);
        }
    }
    else if(ordre == 1) {
        if (signe==std::string("plus")) {
            return G_ordre1(u, up);
        }
        else if (signe==std::string("moins")) {
            return G_ordre1(um, u);
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


double TransportScheme::F_ordre2(double um, double up)
{   
    double ux = _data->u[0]; 
    double c = ux * _data->dt / _data->h;

    return 0.5 * ux * (up + um - c * (up - um));
}

double TransportScheme::G_ordre2(double um, double up)
{   
    double uy = _data->u[1]; 
    double c = uy * _data->dt / _data->h; 

    return 0.5 * uy * (up + um - c * (up - um));
}

double TransportScheme::F_ordre1(double um, double up)
{   
    double ux(_data->u[0]);
    if (ux>=0) {
        return ux*um;
    }
    else {
        return ux*up;
    }
}

double TransportScheme::G_ordre1(double um, double up)
{   
    double uy(_data->u[1]);
    if (uy>=0) {
        return uy*um;
    }
    else {
        return uy*up;
    }
}