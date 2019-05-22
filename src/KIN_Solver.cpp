//! @file KIN_Solver.cpp

// This file is part of OpenMKM. See License.txt in the top-level directory.

#include "KIN_Solver.h"
#include "cantera/numerics/NumUtil.h"
#include "cantera/base/stringUtils.h"

#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "ida/ida.h"
#if CT_SUNDIALS_VERSION >= 30
    #if CT_SUNDIALS_USE_LAPACK
        #include "sunlinsol/sunlinsol_lapackdense.h"
        #include "sunlinsol/sunlinsol_lapackband.h"
    #else
        #include "sunlinsol/sunlinsol_dense.h"
        #include "sunlinsol/sunlinsol_band.h"
    #endif
    #include "sunlinsol/sunlinsol_spgmr.h"
    #include "kinsol/kinsol_direct.h"
    #include "kinsol/kinsol_spils.h"
#else
    #include "kinsol/kinsol_dense.h"
    #include "ida/kinsol_spgmr.h"
    #include "ida/kinsol_band.h"
#endif
#include "nvector/nvector_serial.h"

using namespace std;

#if CT_SUNDIALS_VERSION < 25
typedef int sd_size_t;
#else
typedef long int sd_size_t;
#endif

namespace Cantera
{

//! A simple class to hold an array of parameter values and a pointer to an
//! instance of a subclass of ResidEval.
class ResidData
{
public:
    ResidData(ResidJacEval* f, KIN_Solver* s, int npar = 0) {
        m_func = f;
        m_solver = s;
    }

    virtual ~ResidData() {
    }

    ResidJacEval* m_func;
    KIN_Solver* m_solver;
};
}

extern "C" {
    //! Function called by KIN to evaluate the residual, given y.
    /*!
     * KIN allows passing in a void* pointer to access external data. Instead of
     * requiring the user to provide a residual function directly to KIN (which
     * would require using the sundials data types N_Vector, etc.), we define
     * this function as the single function that KIN always calls. The real
     * evaluation of the residual is done by an instance of a subclass of
     * ResidEval, passed in to this function as a pointer in the parameters.
     *
     * FROM KIN WRITEUP -> What the KIN solver expects as a return flag from its
     * residual routines:
     *
     * A KINResFn res should return a value of 0 if successful, a positive value
     * if a recoverable error occured (e.g. yy has an illegal value), or a
     * negative value if a nonrecoverable error occured. In the latter case, the
     * program halts. If a recoverable error occured, the integrator will
     * attempt to correct and retry.
     */
    static int kin_resid(realtype t, N_Vector y,  N_Vector r, void* f_data)
    {
        Cantera::ResidData* d = (Cantera::ResidData*) f_data;
        Cantera::ResidJacEval* f = d->m_func;
        Cantera::KIN_Solver* s = d->m_solver;
        // TODO evaluate evalType. Assumed to be Base_ResidEval
        int flag = f->evalResidNJ(NV_DATA_S(y), NV_DATA_S(r));
        if (flag < 0) {
            // This signals to KIN that a nonrecoverable error has occurred.
            return flag;
        } else {
            return 0;
        }
    }

    //! Function called by by KIN to evaluate the Jacobian, given y.
    /*!
     * typedef int (*KINDlsDenseJacFn)(sd_size_t N, realtype t, realtype c_j,
     *                             N_Vector y, N_Vector yp, N_Vector r,
     *                             DlsMat Jac, void *user_data,
     *                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
     *
     * A KINDlsDenseJacFn should return
     * - 0 if successful,
     * - a positive int if a recoverable error occurred, or
     * - a negative int if a nonrecoverable error occurred.
     *
     * In the case of a recoverable error return, the integrator will attempt to
     * recover by reducing the stepsize (which changes cj).
     */
#if CT_SUNDIALS_VERSION >= 30
    static int kin_jacobian(realtype c_j, N_Vector y, N_Vector yp,
                            N_Vector r, SUNMatrix Jac, void *f_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    {
        Cantera::ResidData* d = (Cantera::ResidData*) f_data;
        Cantera::ResidJacEval* f = d->m_func;
        Cantera::KIN_Solver* s = d->m_solver;
        double** cols;
        if (SUNMatGetID(Jac) == SUNMATRIX_DENSE) {
            cols = SM_COLS_D(Jac);
        } else if (SUNMatGetID(Jac) == SUNMATRIX_BAND) {
            cols = SM_COLS_B(Jac);
        } else {
            return 1; // Unknown SUNMatrix type
        }
        f->evalJacobianDP(c_j, NV_DATA_S(y), NV_DATA_S(yp),
                          cols, NV_DATA_S(r));
        return 0;
    }
#else
    static int kin_jacobian(sd_size_t nrows, realtype c_j, N_Vector y, N_Vector ydot, N_Vector r,
                            DlsMat Jac, void* f_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    {
        Cantera::ResidData* d = (Cantera::ResidData*) f_data;
        Cantera::ResidJacEval* f = d->m_func;
        Cantera::KIN_Solver* s = d->m_solver;
        f->evalJacobianDP(c_j, NV_DATA_S(y), NV_DATA_S(ydot),
                          Jac->cols, NV_DATA_S(r));
        return 0;
    }
#endif

}

namespace Cantera
{

KIN_Solver::KIN_Solver(ResidJacEval& f) :
    NonLinearSolver(f),
    m_kin_mem(0),
    m_y(0),
    m_ydot(0),
    m_id(0),
    m_constraints(0),
    m_abstol(0),
    m_type(0),
    m_itol(KIN_SS),
    m_iter(0),
    m_reltol(1.e-9),
    m_abstols(1.e-15),
    m_nabs(0),
    m_maxsteps(20000),
    m_maxord(0),
    m_formJac(0),
    m_maxErrTestFails(-1),
    m_maxNonlinIters(0),
    m_maxNonlinConvFails(-1),
    m_setSuppressAlg(0),
    m_mupper(0),
    m_mlower(0)
{
}

KIN_Solver::~KIN_Solver()
{
    if (m_kin_mem) {
        KINFree(&m_kin_mem);
    }
    if (m_y) {
        N_VDestroy_Serial(m_y);
    }
    if (m_ydot) {
        N_VDestroy_Serial(m_ydot);
    }
    if (m_abstol) {
        N_VDestroy_Serial(m_abstol);
    }
    if (m_constraints) {
        N_VDestroy_Serial(m_constraints);
    }
}

doublereal KIN_Solver::solution(int k) const
{
    return NV_Ith_S(m_y,k);
}

const doublereal* KIN_Solver::solutionVector() const
{
    return NV_DATA_S(m_y);
}

void KIN_Solver::setTolerances(double reltol, double* abstol)
{
    m_itol = KIN_SV;
    if (!m_abstol) {
        m_abstol = N_VNew_Serial(m_neq);
    }
    for (int i = 0; i < m_neq; i++) {
        NV_Ith_S(m_abstol, i) = abstol[i];
    }
    m_reltol = reltol;
    if (m_kin_mem) {
        int flag = KINSVtolerances(m_kin_mem, m_reltol, m_abstol);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::setTolerances",
                               "Memory allocation failed.");
        }
    }
}

void KIN_Solver::setTolerances(doublereal reltol, doublereal abstol)
{
    m_itol = KIN_SS;
    m_reltol = reltol;
    m_abstols = abstol;
    if (m_kin_mem) {
        int flag = KINSStolerances(m_kin_mem, m_reltol, m_abstols);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::setTolerances",
                               "Memory allocation failed.");
        }
    }
}

void KIN_Solver::setLinearSolverType(int solverType)
{
    m_type = solverType;
}

void KIN_Solver::setDenseLinearSolver()
{
    setLinearSolverType(0);
}

void KIN_Solver::setBandedLinearSolver(int m_upper, int m_lower)
{
    m_type = 2;
    m_upper = m_mupper;
    m_mlower = m_lower;
}

void KIN_Solver::setMaxOrder(int n)
{
    m_maxord = n;
}

void KIN_Solver::setMaxNumSteps(int n)
{
    m_maxsteps = n;
}

void KIN_Solver::setInitialStepSize(doublereal h0)
{
    m_h0 = h0;
}

void KIN_Solver::setStopTime(doublereal tstop)
{
    m_tstop = tstop;
}

void KIN_Solver::setJacobianType(int formJac)
{
    m_formJac = formJac;
    if (m_kin_mem && m_formJac == 1) {
        #if CT_SUNDIALS_VERSION >= 30
            int flag = KINDlsSetJacFn(m_kin_mem, ida_jacobian);
        #else
            int flag = KINDlsSetDenseJacFn(m_kin_mem, ida_jacobian);
        #endif
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::setJacobianType",
                               "KINDlsSetDenseJacFn failed.");
        }
    }
}

void KIN_Solver::setConstraint(const int k, const int constraintFlag)
{
    if (checkFlag(constraintFlag)){
        if(!m_constraints) {
            m_constraints = N_VNew_Serial(m_neq);
        }
        NV_Ith_S(m_constraints, k) = constraintFlag;
        if (m_kin_mem){
            auto flag = KINSetConstraints(m_kin_mem, m_constraints);
            if (flag != KIN_SUCCESS) {
                throw CanteraError("KIN_Solver::setConstraint", 
                                   "KINSetConstraint failed.");
            }
        }
    } else { 
        throw CanteraError("KIN_Solver::setConstraint", 
                           "Invalid Constraint vaue");
    }
}

void KIN_Solver::setConstraints(const int * const constraintFlags)
{
    if(!m_constraints) {
        m_constraints = N_VNew_Serial(m_neq);
    }
    for(size_t i = 0; i < m_neq; i++){
        auto cflag = constraintFlags[i];
        if(!checkFlag(cflag)){
            throw CanteraError("KIN_Solver::setConstraints", 
                               "Invalid Constraint vaue detected");
        }
        NV_Ith_S(m_constraints, i) = cflag;
    }
    if (m_kin_mem){
        auto flag = KINSetConstraints(m_kin_mem, m_constraints);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::setConstraints", 
                               "KINSetConstraint failed.");
        }
    }
}

void KIN_Solver::setMaxErrTestFailures(int maxErrTestFails)
{
    m_maxErrTestFails = maxErrTestFails;
}

void KIN_Solver::setMaxNonlinIterations(int n)
{
    m_maxNonlinIters = n;
}

void KIN_Solver::setMaxNonlinConvFailures(int n)
{
    m_maxNonlinConvFails = n;
}

void KIN_Solver::inclAlgebraicInErrorTest(bool yesno)
{
    if (yesno) {
        m_setSuppressAlg = 0;
    } else {
        m_setSuppressAlg = 1;
    }
}

void KIN_Solver::init()
{
    if (m_y) {
        N_VDestroy_Serial(m_y);
    }
    if (m_ydot) {
        N_VDestroy_Serial(m_ydot);
    }
    if (m_id) {
        N_VDestroy_Serial(m_id);
    }
    if (m_constraints) {
        N_VDestroy_Serial(m_constraints);
    }

    m_y = N_VNew_Serial(m_neq);
    m_ydot = N_VNew_Serial(m_neq);
    m_constraints = N_VNew_Serial(m_neq);

    for (int i=0; i<m_neq; i++) {
        NV_Ith_S(m_y, i) = 0.0;
        NV_Ith_S(m_ydot, i) = 0.0;
        NV_Ith_S(m_constraints, i) = 0.0;
    }

    // get the initial conditions
    m_resid.getInitialConditions(NV_DATA_S(m_y), NV_DATA_S(m_ydot));

    if (m_kin_mem) {
        KINFree(&m_kin_mem);
    }

    /* Call KINCreate */
    m_kin_mem = KINCreate();

    int flag = KINInit(m_kin_mem, kin_resid, m_y, m_ydot);
    if (flag != KIN_SUCCESS) {
        if (flag == KIN_MEM_FAIL) {
            throw CanteraError("KIN_Solver::init",
                               "Memory allocation failed.");
        } else if (flag == KIN_ILL_INPUT) {
            throw CanteraError("KIN_Solver::init",
                "Illegal value for KINInit input argument.");
        } else {
            throw CanteraError("KIN_Solver::init", "KINInit failed.");
        }
    }
    if (m_itol == KIN_SV) {
        flag = KINSVtolerances(m_kin_mem, m_reltol, m_abstol);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init", "Memory allocation failed.");
        }
    } else {
        flag = KINSStolerances(m_kin_mem, m_reltol, m_abstols);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init", "Memory allocation failed.");
        }
    }

    // set the linear solver type
    if (m_type == 1 || m_type == 0) {
        long int N = m_neq;
        int flag;
        #if CT_SUNDIALS_VERSION >= 30
            //SUNLinSolFree((SUNLinearSolver) m_linsol);
            //SUNMatDestroy((SUNMatrix) m_linsol_matrix);
            m_linsol_matrix = SUNDenseMatrix(N, N);
            if (m_linsol_matrix == nullptr) {
                throw CanteraError("KIN_Solver::init",
                    "Unable to create SUNDenseMatrix of size {0} x {0}", N);
            }
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLapackDense(m_y, (SUNMatrix) m_linsol_matrix);
            #else
                m_linsol = SUNDenseLinearSolver(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
            flag = KINDlsSetLinearSolver(m_kin_mem, (SUNLinearSolver) m_linsol,
                                         (SUNMatrix) m_linsol_matrix);
        #else
            flag = KINDense(m_kin_mem, N);
        #endif
        if (flag) {
            throw CanteraError("KIN_Solver::init", "KINDense failed");
        }
    } else if (m_type == 2) {
        long int N = m_neq;
        long int nu = m_mupper;
        long int nl = m_mlower;
        #if CT_SUNDIALS_VERSION >= 30
            //SUNLinSolFree((SUNLinearSolver) m_linsol);
            //SUNMatDestroy((SUNMatrix) m_linsol_matrix);
            m_linsol_matrix = SUNBandMatrix(N, nu, nl, nu+nl);
            if (m_linsol_matrix == nullptr) {
                throw CanteraError("KIN_Solver::init",
                    "Unable to create SUNBandMatrix of size {} with bandwidths "
                    "{} and {}", N, nu, nl);
            }
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLapackBand(m_y, (SUNMatrix) m_linsol_matrix);
            #else
                m_linsol = SUNBandLinearSolver(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
            KINDlsSetLinearSolver(m_kin_mem, (SUNLinearSolver) m_linsol,
                                  (SUNMatrix) m_linsol_matrix);
        #else
            KINBand(m_kin_mem, N, nu, nl);
        #endif
    } else {
        throw CanteraError("KIN_Solver::init",
                           "unsupported linear solver type");
    }

    if (m_formJac == 1) {
        #if CT_SUNDIALS_VERSION >= 30
            flag = KINDlsSetJacFn(m_kin_mem, ida_jacobian);
        #else
            flag = KINDlsSetDenseJacFn(m_kin_mem, ida_jacobian);
        #endif
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init",
                               "KINDlsSetDenseJacFn failed.");
        }
    }

    // pass a pointer to func in m_data
    m_fdata.reset(new ResidData(&m_resid, this, m_resid.nparams()));
    flag = KINSetUserData(m_kin_mem, m_fdata.get());
    if (flag != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::init", "KINSetUserData failed.");
    }

    // set options
    if (m_maxord > 0) {
        flag = KINSetMaxOrd(m_kin_mem, m_maxord);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init", "KINSetMaxOrd failed.");
        }
    }
    if (m_maxsteps > 0) {
        flag = KINSetMaxNumSteps(m_kin_mem, m_maxsteps);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init", "KINSetMaxNumSteps failed.");
        }
    }
    if (m_h0 > 0.0) {
        flag = KINSetInitStep(m_kin_mem, m_h0);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init", "KINSetInitStep failed.");
        }
    }
    if (m_tstop > 0.0) {
        flag = KINSetStopTime(m_kin_mem, m_tstop);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init", "KINSetStopTime failed.");
        }
    }
    if (m_maxErrTestFails >= 0) {
        flag = KINSetMaxErrTestFails(m_kin_mem, m_maxErrTestFails);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init",
                               "KINSetMaxErrTestFails failed.");
        }
    }
    if (m_maxNonlinIters > 0) {
        flag = KINSetMaxNonlinIters(m_kin_mem, m_maxNonlinIters);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init",
                               "KINSetmaxNonlinIters failed.");
        }
    }
    if (m_maxNonlinConvFails >= 0) {
        flag = KINSetMaxConvFails(m_kin_mem, m_maxNonlinConvFails);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init",
                               "KINSetMaxConvFails failed.");
        }
    }
    if (m_setSuppressAlg != 0) {
        flag = KINSetSuppressAlg(m_kin_mem, m_setSuppressAlg);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init", "KINSetSuppressAlg failed.");
        }
    }
    if (m_resid.nConstraints()){ // Constraints are defined. 
        for (size_t i = 0; i < m_neq; i++) {
            NV_Ith_S(m_constraints, i) = m_resid.constraint(i);
        }
        flag = KINSetConstraints(m_kin_mem, m_constraints);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::init", "KINSetConstraints failed");
        }
    }
}

void KIN_Solver::correctInitial_Y_given_Yp(doublereal* y, doublereal* yp, doublereal tout)
{
    doublereal tout1 = tout;
    if (tout == 0.0) {
        double h0 = 1.0E-5;
        if (m_h0 > 0.0) {
            h0 = m_h0;
        }
        tout1 = m_t0 + h0;
    }

    int flag = KINCalcIC(m_kin_mem, KIN_Y_INIT, tout1);
    if (flag != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::correctInitial_Y_given_Yp",
                           "KINCalcIC failed: error = {}", flag);
    }

    flag = KINGetConsistentIC(m_kin_mem, m_y, m_ydot);
    if (flag != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::correctInitial_Y_given_Yp",
                           "KINGetSolution failed: error = {}", flag);
    }
    for (int i = 0; i < m_neq; i++) {
        y[i] = NV_Ith_S(m_y, i);
        yp[i] = NV_Ith_S(m_ydot, i);
    }
}

void KIN_Solver::correctInitial_YaYp_given_Yd(doublereal* y, doublereal* yp, doublereal tout)
{
    int icopt = KIN_YA_YDP_INIT;
    doublereal tout1 = tout;
    if (tout == 0.0) {
        double h0 = 1.0E-5;
        if (m_h0 > 0.0) {
            h0 = m_h0;
        }
        tout1 = m_t0 + h0;
    }

    int flag = KINCalcIC(m_kin_mem, icopt, tout1);
    if (flag != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::correctInitial_YaYp_given_Yd",
                           "KINCalcIC failed: error = {}", flag);
    }

    flag = KINGetConsistentIC(m_kin_mem, m_y, m_ydot);
    if (flag != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::correctInitial_YaYp_given_Yd",
                           "KINGetSolution failed: error = {}", flag);
    }
    for (int i = 0; i < m_neq; i++) {
        y[i] = NV_Ith_S(m_y, i);
        yp[i] = NV_Ith_S(m_ydot, i);
    }
}

int KIN_Solver::solve(double tout)
{
    double tretn = tout - 1000;
    int flag;
    flag = KINSetStopTime(m_kin_mem, tout);
    if (flag != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::solve", "KIN error encountered.");
    }
    while (tretn < tout) {
        if (tout <= m_tcurrent) {
            throw CanteraError("KIN_Solver::solve", "tout <= tcurrent");
        }
        m_told_old = m_told;
        m_told = m_tcurrent;
        flag = KINSolve(m_kin_mem, tout, &tretn, m_y, m_ydot, KIN_ONE_STEP);
        if (flag < 0) {
            throw CanteraError("KIN_Solver::solve", "KIN error encountered.");
        } else if (flag == KIN_TSTOP_RETURN) {
            // we've reached our goal, and have actually integrated past it
        } else if (flag == KIN_ROOT_RETURN) {
            // not sure what to do with this yet
        } else if (flag == KIN_WARNING) {
            throw CanteraError("KIN_Solver::solve", "KIN Warning encountered.");
        }
        m_tcurrent = tretn;
        m_deltat = m_tcurrent - m_told;
    };

    if (flag != KIN_SUCCESS && flag != KIN_TSTOP_RETURN) {
        throw CanteraError("KIN_Solver::solve", "KIN error encountered.");
    }
    return flag;
}

double KIN_Solver::step(double tout)
{
    double t;
    int flag;
    if (tout <= m_tcurrent) {
        throw CanteraError("KIN_Solver::step", "tout <= tcurrent");
    }
    m_told_old = m_told;
    m_told = m_tcurrent;
    flag = KINSolve(m_kin_mem, tout, &t, m_y, m_ydot, KIN_ONE_STEP);
    if (flag < 0) {
        throw CanteraError("KIN_Solver::step", "KIN error encountered.");
    } else if (flag == KIN_TSTOP_RETURN) {
        // we've reached our goal, and have actually integrated past it
    } else if (flag == KIN_ROOT_RETURN) {
        // not sure what to do with this yet
    } else if (flag == KIN_WARNING) {
        throw CanteraError("KIN_Solver::step", "KIN Warning encountered.");
    }
    m_tcurrent = t;
    m_deltat = m_tcurrent - m_told;
    return t;
}

doublereal KIN_Solver::getOutputParameter(int flag) const
{
    long int lenrw, leniw;
    switch (flag) {
    case REAL_WORKSPACE_SIZE:
        flag = KINGetWorkSpace(m_kin_mem, &lenrw, &leniw);
        return doublereal(lenrw);
        break;
    }
    return 0.0;
}

}
