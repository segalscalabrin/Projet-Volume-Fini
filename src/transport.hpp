#ifndef TRANSPORT_HPP
#define TRANSPORT_HPP

#include "include.hpp"

class TransportScheme
{
    private:
        Data *_data;
        Grid *_grid;
        AGeometry  *_geo;

        PiercedVector<double> _phi;
        PiercedVector<double> _phiExact;

    public:
        TransportScheme(Data *data, Grid *grid, AGeometry *geo) : _data(data), _grid(grid), _geo(geo) {}
        ~TransportScheme() {}

        void initializePhi();
        const PiercedVector<double>& getPhi();
        const PiercedVector<double>& getPhiExact();

        void computePhi();

        std::array<long, 4> getCellNeighs(long cellId);

        bool critereMood(long cellId, double new_u);

        double cellBorderCheck(long cellId, long cellToCheckId, std::string signe, std::string direction);

        std::array<double, 4> computeCellFlux(int cellId, int ordre);

        double computeFlux(long cellmId, long cellId, long cellpId, int ordre, std::string signe, std::string direction);

        double Flux_F(double um, double u, double up, int ordre, std::string signe);
        double Flux_G(double um, double u, double up, int ordre, std::string signe);

        double F_ordre2(double um, double up);
        double G_ordre2(double um, double up);
        double F_ordre1(double um, double up);
        double G_ordre1(double um, double up);

};

#endif
