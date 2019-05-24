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



void ReactorNetHybrid::setMaxErrTestFails(int nmax)
{
    m_maxErrTestFails = nmax;
    m_init = false;
}



void ReactorNetHybrid::initialize()
{
    m_nv = 0;
    debuglog("Initializing reactor network.\n", m_verbose);
    if (m_reactors.empty()) {
        throw CanteraError("ReactorNetHybrid::initialize",
                           "no reactors in network!");
    }
    m_start.assign(1, 0);
    for (size_t n = 0; n < m_reactors.size(); n++) {
        Reactor& r = *m_reactors[n];
        r.initialize(m_time);
        size_t nv = r.neq();
        m_nv += nv;
        m_start.push_back(m_nv);

        if (m_verbose) {
            writelog("Reactor {:d}: {:d} variables.\n", n, nv);
            writelog("              {:d} sensitivity params.\n", r.nSensParams());
        }
        if (r.type() == FlowReactorType && m_reactors.size() > 1) {
            throw CanteraError("ReactorNetHybrid::initialize",
                               "FlowReactors must be used alone.");
        }
    }

    m_ydot.resize(m_nv,0.0);
    m_atol.resize(neq());
    fill(m_atol.begin(), m_atol.end(), m_atols);
    m_integ->setTolerances(m_rtol, neq(), m_atol.data());
    m_integ->setSensitivityTolerances(m_rtolsens, m_atolsens);
    m_integ->setMaxStepSize(m_maxstep);
    m_integ->setMaxErrTestFails(m_maxErrTestFails);
    if (m_verbose) {
        writelog("Number of equations: {:d}\n", neq());
        writelog("Maximum time step:   {:14.6g}\n", m_maxstep);
    }
    m_integ->initialize(m_time, *this);
    m_integrator_init = true;
    m_init = true;
}

void ReactorNetHybrid::reinitialize()
{
    if (m_init) {
        debuglog("Re-initializing reactor network.\n", m_verbose);
        m_integ->reinitialize(m_time, *this);
        m_integrator_init = true;
    } else {
        initialize();
    }
}

void ReactorHybrid::nonlinSolverInitialize()
{
    debuglog("Initializing nonlinear solver.\n", m_verbose);
    
    m_nonlinsol->setMaxStepSize(0.5); //TODO: Find a good value
    m_nonlinsol->setMaxSteps(300);    //TODO: Find a good value
    m_nonlin_sol->initialize(*this);
    m_nonlin_sol_init = true;
    m_steady_sate = false;
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
    nls_conv_flag = false;
    if (!m_nonlin_sol_init)
        nonlinsol_initialize();
    else
        nonlinsol_reinitialize();
    vector<double> tmp_y(neq());
    while (!nls_conv_flag) {
        getState(tmp_y.data());
        int nonconv_status = m_nonlin_sol->solve(tmp_y);
        if (nonconv_status) {
            auto new_time = m_time ? m_time * 2 : 1e-6;
            if (new_time > m_final_time) {
                break;
            }
            advance(new_time);
        }
    }
    m_steady_state = true;

}

void ReactorNetHybrid::eval(doublereal t, doublereal* y,
                      doublereal* ydot, doublereal* p)
{
    updateState(y);
    for (size_t n = 0; n < m_reactors.size(); n++) {
        m_reactors[n]->evalEqs(t, y + m_start[n], ydot + m_start[n], p);
    }
    checkFinite("ydot", ydot, m_nv);
}


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
