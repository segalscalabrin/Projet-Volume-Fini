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

    computeFlux(flux, 1);

    for (auto &cell : _grid->getCells()) {
        const long &id = cell.getId();
        double delta;
        std::vector<long> neighborsId;
        std::vector<double> neighborsPhi;

        new_phi[id] = _phi[id] - (_data->dt / _data->h) * (flux[1][id] - flux[0][id])
                               - (_data->dt / _data->h) * (flux[3][id] - flux[2][id]);
    }

    _phi = new_phi;
}



void TransportScheme::computeFlux(std::array<PiercedVector<double>, 4>& flux, int ordre)
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
                f = Flux_F(ordre, _phi[ownersId[1]], _phi[ownersId[0]]);
                std::cout << f << std::endl;
                flux[1][ownersId[0]] = f;
                flux[0][ownersId[1]] = f;
            } 
            else if (center1[NPX] > center2[NPX] && center1[NPZ] == center2[NPZ]) {
                f = Flux_F(ordre, _phi[ownersId[0]], _phi[ownersId[1]]);
                flux[0][ownersId[0]] = f;
                flux[1][ownersId[1]] = f;
            } 
            else if (center1[NPY] < center2[NPY] && center1[NPZ] == center2[NPZ]) {
                g = Flux_G(ordre, _phi[ownersId[1]], _phi[ownersId[0]]);
                flux[3][ownersId[0]] = g;
                flux[2][ownersId[1]] = g;
            } 
            else if (center1[NPY] > center2[NPY] && center1[NPZ] == center2[NPZ]) {
                g = Flux_G(ordre, _phi[ownersId[0]], _phi[ownersId[1]]);
                flux[2][ownersId[0]] = g;
                flux[3][ownersId[1]] = g;
            }
        }
        else {
            center2 = _grid->evalInterfaceCentroid(inter.getId());
            f = Flux_F(ordre, _phi[ownersId[0]], _phi[ownersId[0]]);
            g = Flux_G(ordre, _phi[ownersId[0]], _phi[ownersId[0]]);
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


double TransportScheme::Flux_F(int ordre, double up, double um)
{   
    if (ordre == 3) {
        std::cerr << "pas implémenté" << std::endl;
        exit(1);
    }
    else if(ordre == 2) {
        return F_ordre2(up, um);
    }
    else {
        return F_ordre1(up, um);
    }
}

double TransportScheme::Flux_G(int ordre, double up, double um)
{   
    if (ordre == 3) {
        std::cerr << "pas implémenté" << std::endl;
        exit(1);
    }
    else if(ordre == 2) {
        return F_ordre2(up, um);
    }
    else {
        return F_ordre1(up, um);
    }
}


double TransportScheme::F_ordre2(double up, double um)
{   
    double ux = _data->u[0]; 
    double c = ux * _data->dt / _data->h;

    return 0.5 * ux * (up + um - c * (up - um));
}

double TransportScheme::G_ordre2(double up, double um)
{   
    double uy = _data->u[1]; 
    double c = uy * _data->dt / _data->h; 

    return 0.5 * uy * (up + um - c * (up - um));
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