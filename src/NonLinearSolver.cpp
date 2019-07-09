//! @file NonLinearSolver.cpp

// This file is part of OpenMKM. See License.txt in the top-level directory 
// for license and copyright information.

#include "cantera/base/ct_defs.h"
#include "NonLinearSolver.h"
#include "KIN_Solver.h"

namespace Cantera   // Cantera namespace for convenience
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
