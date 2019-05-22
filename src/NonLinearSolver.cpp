//! @file NonLinearsolver.cpp Factory routine for picking the NonLinear solver package

// This file is part of OpenMKM. See License.txt in the top-level directory.

#include "cantera/base/ct_defs.h"
#include "cantera/numerics/NonLinearSolver.h"
#include "cantera/numerics/KIN_Solver.h"

// DAE_DEVEL is turned off at the current time
#define DAE_DEVEL
#ifdef DAE_DEVEL

namespace Cantera
{
NonLinearSolver* newNonLinearSolver(const std::string& itype, ResidJacEval& f)
{
    if (itype == "KIN") {
        return new KIN_Solver(f);
    } else {
        throw CanteraError("newNonLinearSolver",
                           "unknown Nonlinear solver: "+itype);
    }
}
}

#endif
