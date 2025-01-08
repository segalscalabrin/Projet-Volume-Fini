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
    return;
}

