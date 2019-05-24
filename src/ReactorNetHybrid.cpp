//! @file ReactorNetHybrid.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "ReactorNetHybrid.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

ReactorNetHybrid::ReactorNetHybrid() :
    ReactorNet(), 
    m_nonlin_sol(newNonLinearSolver("KINSOL")), 
    m_nonlin_sol_init(false), 
    m_steady_state(false)
{
}



void ReactorNetHybrid::initialize()
{
    ReactorNet::initialize();
    m_nonlin_sol_init = false;
}

void ReactorNetHybrid::reinitialize()
{
    if (m_init) {
        debuglog("Re-initializing reactor network.\n", m_verbose);
        m_integ->reinitialize(m_time, *this);
        m_integrator_init = true;
        m_nonlin_sol_init = false;
    } else {
        initialize();
    }
}

void ReactorNetHybrid::nonlinSolverInitialize()
{
    debuglog("Initializing nonlinear solver.\n", m_verbose);
    
    m_nonlin_sol->setMaxStepSize(0.5); //TODO: Find a good value
    m_nonlin_sol->setMaxSteps(300);    //TODO: Find a good value
    m_nonlin_sol->initialize(*this);
    m_nonlin_sol_init = true;
    m_steady_state = false;
}

void ReactorNetHybrid::solve()
{
    if (!m_steady_state) {
        // Advance the reactor with ode integrator first.
        auto new_time = m_time ? m_time * 2 : 1e-6;
        advance(new_time);
    }

    // Try nonlinear solver
    // If nonlinear solver fails, advance the time with integrator
    // till the nonlinear solver works
    auto nls_conv_flag = false;
    if (!m_nonlin_sol_init)
        nonlinSolverInitialize();
    //else
    //    nonlinsol_reinitialize();
    //vector<double> tmp_y(neq());
    while (!nls_conv_flag) {
        //getState(tmp_y.data());
        int nonconv_status = m_nonlin_sol->solve();
        if (nonconv_status) {
            auto new_time = m_time ? m_time * 2 : 1e-6;
            if (new_time > m_final_time) {
                break;
            }
            advance(new_time);
        } else {
            nls_conv_flag = true;
        }
    }
    m_steady_state = true;

}

/*
void ReactorNetHybrid::eval(doublereal t, doublereal* y,
                      doublereal* ydot, doublereal* p)
{
    updateState(y);
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->evalEqs(t, y + m_start[n], ydot + m_start[n], p);
    }
    checkFinite("ydot", ydot, m_nv);
}
*/


/*
void ReactorNetHybrid::evalJacobian(doublereal t, doublereal* y,
                              doublereal* ydot, doublereal* p, Array2D* j)
{
    //evaluate the unperturbed ydot
    eval(t, y, ydot, p);
    for (size_t n = 0; n < m_nv; n++) {
        // perturb x(n)
        double ysave = y[n];
        double dy = m_atol[n] + fabs(ysave)*m_rtol;
        y[n] = ysave + dy;
        dy = y[n] - ysave;

        // calculate perturbed residual
        eval(t, y, m_ydot.data(), p);

        // compute nth column of Jacobian
        for (size_t m = 0; m < m_nv; m++) {
            j->value(m,n) = (m_ydot[m] - ydot[m])/dy;
        }
        y[n] = ysave;
    }
}
*/




}
