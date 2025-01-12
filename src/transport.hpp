#ifndef TRANSPORT_HPP
#define TRANSPORT_HPP

#include "include.hpp"

#include "fonctions.hpp"

class TransportScheme
{
    private:
        Data *_data;
        Grid *_grid;
        AGeometry  *_geo;

        PiercedVector<double> _phi;
        PiercedVector<double> _phiPrime;
        PiercedVector<double> _phiSecond;
        PiercedVector<double> _phiExact;
        double _t = 0; 
        double _dt = 0.1;

    public:
        TransportScheme(Data *data, Grid *grid, AGeometry *geo) : _data(data), _grid(grid), _geo(geo) {}
        ~TransportScheme() {}

        void initializePhi();
        const PiercedVector<double>& getPhi();
        const PiercedVector<double>& getPhiExact();
        const double& getT();
        const double& getDT();

        void computePhi();

        bool critereMood(long cellId, double new_u);

        PiercedVector<double> computePhiRK(int ordre);
        double computeNewPhi(double phi, std::array<double, 8> borderPhi, int ordre);
        double computeFlux(double phi, std::array<double, 8> borderPhi, int ordre);
        std::array<double, 8> computeNeighValue(int cellId, int RK);
        std::array<long, 8> getCellNeighs(long cellId);
        std::array<long, 4> orderedNeigh(long cellId);
        double cellBorderCheck(long cellId, long cellToCheckId, int RK, std::string signe, std::string direction);
};

#endif
