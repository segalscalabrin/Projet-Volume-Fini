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

        void compute();

};

#endif
