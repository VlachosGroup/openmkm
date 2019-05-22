//! @file ODE_integrators.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/ct_defs.h"
#include "NonLinearSolver.h"
#include "KIN_Solver.h"

namespace Cantera
{

NonLinearSolver* newNonLinearSolver(const std::string& itype)
{
    if (itype == "KINSOL") {
        return new KIN_Solver();
    } else {
        throw CanteraError("newNonLinearSolver",
                           "unknown nonlinear solver: "+itype);
    }
    return 0;
}

}

