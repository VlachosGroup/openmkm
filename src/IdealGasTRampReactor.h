//! @file IdealGasTRampReactor.h

// This file is part of openMKM. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef HCT_IDEALGASTRAMPREACTOR_H
#define CT_IDEALGASTRAMPREACTOR_H

#include "cantera/zeroD/IdealGasReactor.h"

namespace HeteroCt
{

/**
 * Class IdealGasTRampReactor is a class for stirred reactors that is specifically
 * optimized for ideal gases and supports linear temperature scaling if chosen. 
 */
class IdealGasTRampReactor : public Cantera::IdealGasReactor
{
public:
    IdealGasTRampReactor(double beta = 0.0) : m_beta(beta) {}

    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    void setBeta(double beta) { //! Set the temperature ramp
        m_beta = beta;
    }

protected:
    double m_beta;  //!< T_ramp value
};

}

#endif
