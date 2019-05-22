/**
 *  @file NonLinearSolver.h
 *  Header file for class NonLinearSolver
 */

// This file is part of OpenMKM. See License.txt in the top-level directory.

#ifndef CT_NONLINEAR_SOLVER_H
#define CT_NONLINEAR_SOLVER_H

#include "ResidJacEval.h"
#include "cantera/base/global.h"

namespace Cantera
{



/**
 * Wrapper for Nonlinear solvers
 *
 */
class NonLinearSolver
{
public:
    NonLinearSolver(FuncEval& f) :
        m_func(f),
        m_neq(f.nEquations()) {
    }

    virtual ~NonLinearSolver() {}

    /**
     * Set error tolerances. This version specifies a scalar
     * relative tolerance, and a vector absolute tolerance.
     */
    virtual void setTolerances(doublereal reltol,
                               doublereal* abstol) {
        warn("setTolerances");
    }

    /**
     * Set error tolerances. This version specifies a scalar
     * relative tolerance, and a scalar absolute tolerance.
     */
    virtual void setTolerances(doublereal reltol, doublereal abstol) {
        warn("setTolerances");
    }

    /**
     * Specify a Jacobian evaluator. If this method is not called,
     * the Jacobian will be computed by finite difference.
     */
    void setJacobian(Jacobian& jac) {
        warn("setJacobian");
    }

    virtual void setLinearSolverType(int solverType) {
        warn("setLinearSolverType");
    }

    virtual void setDenseLinearSolver() {
        warn("setDenseLinearSolver");
    }

    virtual void setBandedLinearSolver(int m_upper, int m_lower) {
        warn("setBandedLinearSolver");
    }
    virtual void setMaxStepSize(doublereal dtmax) {
        warn("setMaxStepSize");
    }
    virtual void setMaxOrder(int n) {
        warn("setMaxOrder");
    }
    virtual void setMaxNumSteps(int n) {
        warn("setMaxNumSteps");
    }
    virtual void setInitialStepSize(doublereal h0) {
        warn("setInitialStepSize");
    }
    virtual void setStopTime(doublereal tstop) {
        warn("setStopTime");
    }
    virtual void setMaxErrTestFailures(int n) {
        warn("setMaxErrTestFailures");
    }
    virtual void setMaxNonlinIterations(int n) {
        warn("setMaxNonlinIterations");
    }
    virtual void setMaxNonlinConvFailures(int n) {
        warn("setMaxNonlinConvFailures");
    }
    virtual void inclAlgebraicInErrorTest(bool yesno) {
        warn("inclAlgebraicInErrorTest");
    }

    /**
     * Solve the system of equations.
     */
    virtual int solve() {
        warn("solve");
        return 0;
    }


    /// Number of equations.
    int nEquations() const {
        return m_func.nEquations();
    }

    /**
     * initialize. Base class method does nothing.
     */
    virtual void init(doublereal t0) {}

    /**
     * Set a solver-specific input parameter.
     */
    virtual void setInputParameter(int flag, doublereal value) {
        warn("setInputParameter");
    }

    /**
     * Get the value of a solver-specific output parameter.
     */
    virtual doublereal getOutputParameter(int flag) const {
        warn("getOutputParameter");
        return 0.0;
    }

    /// the current value of solution component k.
    virtual doublereal solution(int k) const {
        warn("solution");
        return 0.0;
    }

    virtual const doublereal* solutionVector() const {
        warn("solutionVector");
        return &m_dummy;
    }


protected:
    doublereal m_dummy;
    FuncEval& m_func;

    //! Number of total equations in the system
    integer m_neq;

private:
    void warn(const std::string& msg) const {
        writelog(">>>> Warning: method "+msg+" of base class "
                 +"NonLinearSolver called. Nothing done.\n");
    }
};


//! Factor method for choosing a DAE solver
/*!
 * @param itype  String identifying the type (KIN is the only option)
 * @param f      Residual function to be solved by the Nonlinear solver
 * @returns a point to the instantiated NonLinearSolver object
 */
NonLinearSolver* newNonLinearSolver(const std::string& itype, ResidJacEval& f);

}

#endif
