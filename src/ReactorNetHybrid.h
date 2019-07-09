//! @file ReactorNetHybrid.h

// This file is part of OpenMKM. See License.txt in the top-level directory 
// for license and copyright information.

#ifndef OMKM_REACTORNETHYBRID_H
#define OMKM_REACTORNETHYBRID_H

#include "cantera/zeroD/ReactorNet.h"
#include "KIN_Solver.h"
#include "NonLinearSolver.h"

namespace Cantera   // Cantera namespace used for convenience
{

//! A class representing a network of connected reactors.
/*!
 *  This class is used to solve the steady state of governing equations for
 *  a network of reactors (Reactor, ConstPressureReactor) connected by various
 *  means, e.g. Wall, MassFlowController, Valve, PressureController.
 */
class ReactorNetHybrid : public ReactorNet
{
public:
    ReactorNetHybrid();
    virtual ~ReactorNetHybrid() {};
    ReactorNetHybrid(const ReactorNetHybrid&) = delete;
    ReactorNetHybrid& operator=(const ReactorNetHybrid&) = delete;

    void solve();
    void nonlinSolverInitialize();
    void setIntegratorEndTime(double tend) 
    {
        m_final_time = tend;
    }

    void setTolerances(double rtol, double atol);

    //void reinitialize();

protected:
    
    //! Initialize the reactor network. Called automatically the first time
    //! advance or step or solve is called.
    void initialize();

    std::unique_ptr<NonLinearSolver> m_nonlin_sol;
    bool m_nonlin_sol_init; //!< True if integrator initialization is current
    bool m_steady_state;    //!< True if steady state solver has to be used.
    double m_final_time;    //!< Time to stop the integrator when trying steady state
    double m_ftol;

};
}

#endif
