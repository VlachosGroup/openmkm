/**
 *  @file KIN_Solver.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_KINSOLVER_H
#define CT_KINSOLVER_H

#include "NonLinearSolver.h"
#include "cantera/base/ctexceptions.h"

#include "sundials/sundials_nvector.h"

namespace Cantera
{

/**
 * Wrapper class for 'cvodes' integrator from LLNL.
 *
 * @see FuncEval.h. Classes that use KIN_Solver:
 * ImplicitChem, ImplicitSurfChem, Reactor
 */
class KIN_Solver : public NonLinearSolver
{
public:
    /**
     *  Constructor. Default settings: dense Jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     */
    KIN_Solver();
    virtual ~KIN_Solver();
    //virtual void setTolerances(double reltol, size_t n, double* abstol);
    //virtual void setTolerances(double reltol, double abstol);
    //virtual void setSensitivityTolerances(double reltol, double abstol);
    virtual void setProblemType(int probtype);
    virtual void initialize(FuncEval& func);
    //virtual void reinitialize(FuncEval& func);
    virtual int solve();
    virtual double& solution(size_t k);
    virtual double* solution();
    virtual int nEquations() const {
        return static_cast<int>(m_neq);
    }
    virtual int nEvals() const;
    /*
    virtual void setMaxOrder(int n) {
        m_maxord = n;
    }
    */
    //virtual void setMethod(MethodType t);
    //virtual void setIterator(IterType t);
    virtual void setMaxStepSize(double hmax);
    //virtual void setMinStepSize(double hmin);
    virtual void setMaxSteps(int nmax);
    virtual int maxSteps();
    //virtual void setMaxErrTestFails(int n);
    virtual void setBandwidth(int N_Upper, int N_Lower) {
        m_mupper = N_Upper;
        m_mlower = N_Lower;
    }
    /* No sensitivity with steady state
    virtual int nSensParams() {
        return static_cast<int>(m_np);
    }
    virtual double sensitivity(size_t k, size_t p);
    */

    //! Returns a string listing the weighted error estimates associated
    //! with each solution component.
    //! This information can be used to identify which variables are
    //! responsible for solver failures. 
    //virtual std::string getErrorInfo(int N);

    //! Error message information provide by CVodes
    std::string m_error_message;

    /* Disabled due to older version of CVodes
    virtual void setConstraint(const int k, const int flag);

    virtual void setConstraints(const int * const flags);
    */

protected:
    //! Applies user-specified options to the underlying CVODES solver. Called
    //! during integrator initialization or reinitialization.
    void applyOptions();

private:
    //void sensInit(FuncEval& func);

    size_t m_neq;
    void* m_kin_mem;
    void* m_linsol; //!< Sundials linear solver object
    void* m_linsol_matrix; //!< matrix used by Sundials
    FuncEval* m_func;
    N_Vector m_y;
    N_Vector m_scale;
    //N_Vector m_abstol;
    N_Vector m_constraints;
    int m_type;
    //int m_itol;
    //int m_method;
    //int m_iter;
    //int m_maxord;
    //double m_reltol;
    //double m_abstols;
    //double m_reltolsens, m_abstolsens;
    size_t m_nabs;
    double m_hmax, m_hmin;
    int m_maxsteps;
    int m_maxErrTestFails;
    //N_Vector* m_yS;
    //size_t m_np;
    int m_mupper, m_mlower;

    //! Indicates whether the sensitivities stored in m_yS have been updated
    //! for at the current integrator time.
    //bool m_sens_ok;
};

} // namespace

#endif
