//! @file IdealGasReactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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
    //cout << "Mass " << y[0] <<  " Volume " << y[1] << " Temperature " << y[2] << endl;

    IdealGasReactor::evalEqs(time, y, ydot, params);

    if (!m_energy) {
        ydot[2] = m_beta;
    }

}


}
