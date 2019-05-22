/**
 *  @file KIN_Solver.h
 *  Header file for class KIN_Solver
 */

// This file is part of OpenMKM. See License.txt in the top-level directory.

#ifndef CT_KIN_SOLVER_H
#define CT_KIN_SOLVER_H

#include "NonLinearSolver.h"
//#include "NumUtil.h"

#include "sundials/sundials_nvector.h"

#define REAL_WORKSPACE_SIZE 0

namespace Cantera
{

class ResidData;

/**
 * Wrapper for Sundials KIN solver
 *
 */
class KIN_Solver : public NonLinearSolver
{
public:
    //! Constructor.
    /*!
     * Default settings: dense Jacobian, no user-supplied Jacobian function,
     * Newton iteration.
     *
     * @param f  Function that will supply the time dependent residual to be solved
     */
    KIN_Solver(ResidJacEval& f);

    virtual ~KIN_Solver();

    virtual void setTolerances(doublereal reltol,
                               doublereal* abstol);

    virtual void setTolerances(doublereal reltol, doublereal abstol);

    virtual void setLinearSolverType(int solverType);

    //! Set up the problem to use a dense linear direct solver
    virtual void setDenseLinearSolver();

    //! Set up the problem to use a band solver
    /*!
     * @param m_upper   upper band width of the matrix
     * @param m_lower   lower band width of the matrix
     */
    virtual void setBandedLinearSolver(int m_upper, int m_lower);


    //! Set the form of the Jacobian
    /*!
     * @param formJac  Form of the Jacobian
     *                 0 numerical Jacobian
     *                 1 analytical Jacobian given by the evalJacobianDP() function
     */
    virtual void setJacobianType(int formJac);


    virtual void setMaxErrTestFailures(int n);


    //! Set the maximum number of nonlinear solver convergence failures
    /*!
     * @param n  Value of nonlin failures. If value is exceeded, the calculation terminates.
     */
    virtual void setMaxNonlinConvFailures(int n);


    virtual void inclAlgebraicInErrorTest(bool yesno);

    virtual void setConstraint(const int k, const int flag);  

    virtual void setConstraints(const int * const flags);  

    /**
     * Get the value of a solver-specific output parameter.
     */
    virtual doublereal getOutputParameter(int flag) const;


    //! Solver the system 
    /*!
     * @returns the KINSolve() return flag
     *
     * The return values for KINSolve are described below. (The numerical return
     * values are defined above in this file.) All unsuccessful returns give a
     * negative return value.
     *
     * KIN_SUCCESS
     *   KINSolve succeeded and no roots were found.
     *
     * KIN_ROOT_RETURN:  KINSolve succeeded, and found one or more roots.
     *   If nrtfn > 1, call KINGetRootInfo to see which g_i were found
     *   to have a root at (*tret).
     *
     * KIN_TSTOP_RETURN:
     *   KINSolve returns computed results for the independent variable
     *   value tstop. That is, tstop was reached.
     *
     * KIN_MEM_NULL:
     *   The KIN_mem argument was NULL.
     *
     * KIN_ILL_INPUT:
     *   One of the inputs to KINSolve is illegal. This includes the
     *   situation when a component of the error weight vectors
     *   becomes < 0 during internal stepping.  It also includes the
     *   situation where a root of one of the root functions was found
     *   both at t0 and very near t0.  The ILL_INPUT flag
     *   will also be returned if the linear solver function KIN---
     *   (called by the user after calling KINCreate) failed to set one
     *   of the linear solver-related fields in ida_mem or if the linear
     *   solver's init routine failed. In any case, the user should see
     *   the printed error message for more details.
     *
     * KIN_TOO_MUCH_WORK:
     *   The solver took mxstep internal steps but could not reach tout.
     *   The default value for mxstep is MXSTEP_DEFAULT = 500.
     *
     * KIN_TOO_MUCH_ACC:
     *   The solver could not satisfy the accuracy demanded by the user
     *   for some internal step.
     *
     * KIN_ERR_FAIL:
     *   Error test failures occurred too many times (=MXETF = 10) during
     *   one internal step.
     *
     * KIN_CONV_FAIL:
     *   Convergence test failures occurred too many times (= MXNCF = 10)
     *   during one internal step.
     *
     * KIN_LSETUP_FAIL:
     *   The linear solver's setup routine failed
     *   in an unrecoverable manner.
     *
     * KIN_LSOLVE_FAIL:
     *   The linear solver's solve routine failed
     *   in an unrecoverable manner.
     *
     * KIN_CONSTR_FAIL:
     *    The inequality constraints were violated,
     *    and the solver was unable to recover.
     *
     * KIN_REP_RES_ERR:
     *    The user's residual function repeatedly returned a recoverable
     *    error flag, but the solver was unable to recover.
     *
     * KIN_RES_FAIL:
     *    The user's residual function returned a nonrecoverable error
     *    flag.
     */
    virtual int solve();

    virtual void init();

    virtual doublereal solution(int k) const;

    virtual const doublereal* solutionVector() const;

    void* KINMemory() {
        return m_kin_mem;
    }

protected:
    //! Pointer to the KIN memory for the problem
    void* m_kin_mem;
    void* m_linsol; //!< Sundials linear solver object
    void* m_linsol_matrix; //!< matrix used by Sundials

    //! Current value of the solution vector
    N_Vector m_y;

    N_Vector m_id;
    N_Vector m_constraints;
    N_Vector m_abstol;
    int m_type;

    int m_itol;
    int m_iter;
    doublereal m_reltol;
    doublereal m_abstols;
    int m_nabs;

    //! Form of the Jacobian
    /*!
     *  0 numerical Jacobian created by KIN
     *  1 analytical Jacobian. Must have populated the evalJacobianDP()
     *    function in the ResidJacEval class.
     *  2 numerical Jacobian formed by the ResidJacEval class (unimplemented)
     */
    int m_formJac;

    //! maximum number of error test failures
    int m_maxErrTestFails;

    //! Maximum number of nonlinear solver iterations at one solution
    /*!
     *  If zero, this is the default of 4.
     */
    int m_maxNonlinIters;

    //! Maximum number of nonlinear convergence failures
    int m_maxNonlinConvFails;

    //! If true, the algebraic variables don't contribute to error tolerances
    int m_setSuppressAlg;

    std::unique_ptr<ResidData> m_fdata;
    int m_mupper;
    int m_mlower;
};

}

#endif
