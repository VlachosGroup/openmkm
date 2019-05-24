//! @file ReactorNet.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORNETHYBRID_H
#define CT_REACTORNETHYBRID_H

#include "cantera/zeroD/ReactorNet.h"
#include "KIN_Solver.h"
#include "NonLinearSolver.h"

namespace Cantera
{

//! A class representing a network of connected reactors.
/*!
 *  This class is used to integrate the time-dependent governing equations for
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


protected:
    
    //! Initialize the reactor network. Called automatically the first time
    //! advance or step is called.
    void initialize();

    std::unique_ptr<NonLinearSolver> m_nonlin_sol;
    bool m_nonlin_sol_init; //!< True if integrator initialization is current
    bool m_steady_state;    //!< True if steady state solver has to be used.
    double m_final_time;    //!< Time to stop the integrator when trying steady state

};
}

#endif
