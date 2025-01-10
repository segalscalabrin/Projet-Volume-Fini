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

        void computeFlux(std::array<PiercedVector<double>, 4>& flux, int ordre);

        std::array<long, 4> getCellNeighs(long cellId);

        double Flux_F(int ordre, double up, double um);
        double Flux_G(int ordre, double up, double um);

        double F_ordre2(double up, double um);
        double G_ordre2(double up, double um);
        double F_ordre1(double up, double um);
        double G_ordre1(double up, double um);

};

#endif
