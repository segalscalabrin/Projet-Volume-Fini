#ifndef TRANSPORT_HPP
#define TRANSPORT_HPP

#include "include.hpp"

class TransportScheme
{
    private:
        Data *_data;
        Grid *_grid;

        PiercedVector<double> _phi;

    public:
        TransportScheme(Data *data, Grid *grid) : _data(data), _grid(grid) {}
        ~TransportScheme() {}

        void initializePhi(PiercedVector<double> phi);
        const PiercedVector<double>& getPhi();

        void computePhi();

        void computeFlux(std::array<PiercedVector<double>, 4>& flux);

        std::vector<long> getCellNeighbour(long cellId);
        std::vector<double> getPhiNeighbour(std::vector<long> neighborsId);


        std::vector<double> getFlux(long id, std::vector<double> neighborsPhi);
        std::vector<double> getDelta(long id, std::vector<long> neighborsId);

        double F_ordre1(double up, double um);
        double G_ordre1(double up, double um);

};

#endif
