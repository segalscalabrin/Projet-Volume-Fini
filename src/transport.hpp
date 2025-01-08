#ifndef TRANSPORT_HPP
#define TRANSPORT_HPP

#include "include.hpp"

class TransportScheme
{
    private:
        Data *_data;
        Grid *_grid;
        Geometry *_geo;

        std::vector<double> _phi;

    public:
        TransportScheme(Data *data, Grid *grid, Geometry *geo) : _data(data), _grid(grid), _geo(geo) {}
        ~TransportScheme() {}

        void initializePhi();
        std::vector<double> getPhi();

        void compute();

};

#endif
