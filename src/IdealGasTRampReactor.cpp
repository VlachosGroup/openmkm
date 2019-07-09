//! @file IdealGasTRampReactor.cpp A zero-dimensional reactor

// This file is part of OpenMKM. See License.txt in the top-level directory 
// for license and copyright information.

#include <iostream>

#include "cantera/zeroD/IdealGasReactor.h"
#include "IdealGasTRampReactor.h"

using namespace Cantera;
using namespace std;

namespace OpenMKM
{

void IdealGasTRampReactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    IdealGasReactor::evalEqs(time, y, ydot, params);

    if (!m_energy) {
        ydot[2] = m_beta;
    }
}

}
