//! @file KIN_Solver.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/NumUtil.h"
#include "cantera/base/stringUtils.h"
#include "KIN_Solver.h"

#include <iostream>
using namespace std;

#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"
#include "kinsol/kinsol.h"
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
    //#include "kinsol_diag.h"
    #include "kinsol/kinsol_spils.h"
#else
    #if CT_SUNDIALS_USE_LAPACK
        #include "kinsol/kinsol_lapack.h"
    #else
        #include "kinsol/kinsol_dense.h"
        #include "kinsol/kinsol_band.h"
    #endif
    #include "kinsol/kinsol_diag.h"
    #include "kinsol/kinsol_spgmr.h"
#endif

#define KIN_SS 1
#define KIN_SV 2

#if CT_SUNDIALS_VERSION < 25
typedef int sd_size_t;
#else
typedef long int sd_size_t;
#endif

namespace Cantera
{

extern "C" {
    /**
     * Function called by KINSOL to evaluate rhs given y.  The KINSOL nonlinear solver
     * allows passing in a void* pointer to access external data. This pointer
     * is cast to a pointer to a instance of class FuncEval. The equations to be
     * solved should be specified by deriving a class from FuncEval that
     * evaluates the desired equations.
     * @ingroup odeGroup
     */
    static int kin_rhs(N_Vector y, N_Vector rhs, void* f_data)
    {
        FuncEval* f = (FuncEval*) f_data;
        return f->eval_nothrow(0, NV_DATA_S(y), NV_DATA_S(rhs));
    }

    //! Function called by KINSOL when an error is encountered instead of
    //! writing to stdout. Here, save the error message provided by KINSOL so
    //! that it can be included in the subsequently raised CanteraError.
    static void kin_err(int error_code, const char* module,
                           const char* function, char* msg, void* eh_data)
    {
        KIN_Solver* nonlin_solver = (KIN_Solver*) eh_data;
        nonlin_solver->m_error_message = msg;
        nonlin_solver->m_error_message += "\n";
    }
}

KIN_Solver::KIN_Solver() :
    m_neq(0),
    m_kin_mem(0),
    m_linsol(0),
    m_linsol_matrix(0),
    m_func(0),
    m_y(0),
    m_scale(0),
    m_constraints(0),
    //m_abstol(0),
    m_type(DENSE+NOJAC),
    //m_itol(KIN_SS),
    //m_method(CV_BDF),
    //m_iter(KIN_NEWTON),
    //m_maxord(0),
    //m_reltol(1.e-9),
    //m_abstols(1.e-15),
    //m_reltolsens(1.0e-5),
    //m_abstolsens(1.0e-4),
    m_nabs(0),
    m_hmax(0.0),
    m_hmin(0.0),
    m_maxsteps(20000),
    m_maxErrTestFails(0),
    //m_yS(nullptr),
    //m_np(0),
    m_mupper(0), m_mlower(0)
    //m_sens_ok(false)
{
}

KIN_Solver::~KIN_Solver()
{
    if (m_kin_mem) {
        /*
        if (m_np > 0) {
            CVodeSensFree(m_kin_mem);
        }
        */
        KINFree(&m_kin_mem);
    }

    #if CT_SUNDIALS_VERSION >= 30
        SUNLinSolFree((SUNLinearSolver) m_linsol);
        SUNMatDestroy((SUNMatrix) m_linsol_matrix);
    #endif

    if (m_y) {
        N_VDestroy_Serial(m_y);
    }
    if (m_scale) {
        N_VDestroy_Serial(m_scale);
    }
    /*
    if (m_abstol) {
        N_VDestroy_Serial(m_abstol);
    }
    */
    /*
    if (m_yS) {
        N_VDestroyVectorArray_Serial(m_yS, static_cast<sd_size_t>(m_np));
    }
    */
    if (m_constraints) {
        N_VDestroy_Serial(m_constraints);
    }
}

double& KIN_Solver::solution(size_t k)
{
    return NV_Ith_S(m_y, k);
}

double* KIN_Solver::solution()
{
    return NV_DATA_S(m_y);
}

/*
void KIN_Solver::setTolerances(double reltol, size_t n, double* abstol)
{
    m_itol = KIN_SV;
    m_nabs = n;
    if (n != m_neq) {
        if (m_abstol) {
            N_VDestroy_Serial(m_abstol);
        }
        m_abstol = N_VNew_Serial(static_cast<sd_size_t>(n));
    }
    for (size_t i=0; i<n; i++) {
        NV_Ith_S(m_abstol, i) = abstol[i];
    }
    m_reltol = reltol;
}

void KIN_Solver::setTolerances(double reltol, double abstol)
{
    m_itol = KIN_SS;
    m_reltol = reltol;
    m_abstols = abstol;
}
*/

/*
void KIN_Solver::setSensitivityTolerances(double reltol, double abstol)
{
    m_reltolsens = reltol;
    m_abstolsens = abstol;
}
*/

void KIN_Solver::setProblemType(int probtype)
{
    m_type = probtype;
}

/*
void KIN_Solver::setMethod(MethodType t)
{
    if (t == BDF_Method) {
        m_method = CV_BDF;
    } else if (t == Adams_Method) {
        m_method = CV_ADAMS;
    } else {
        throw CanteraError("KIN_Solver::setMethod", "unknown method");
    }
}
*/

void KIN_Solver::setMaxStepSize(doublereal hmax)
{
    m_hmax = hmax;
    if (m_kin_mem) {
        KINSetMaxNewtonStep(m_kin_mem, hmax);
    }
}

/*void KIN_Solver::setMinStepSize(doublereal hmin)
{
    m_hmin = hmin;
    if (m_kin_mem) {
        KINSetMinStep(m_kin_mem, hmin);
    }
}
*/

void KIN_Solver::setMaxSteps(int nmax)
{
    m_maxsteps = nmax;
    if (m_kin_mem) {
        KINSetNumMaxIters(m_kin_mem, m_maxsteps);
    }
}

int KIN_Solver::maxSteps()
{
    return m_maxsteps;
}

/*
void KIN_Solver::setMaxErrTestFails(int n)
{
    m_maxErrTestFails = n;
    if (m_kin_mem) {
        KINSetMaxErrTestFails(m_kin_mem, n);
    }
}
*/

/*
void KIN_Solver::setIterator(IterType t)
{
    if (t == Newton_Iter) {
        m_iter = KIN_NEWTON;
    } else if (t == Functional_Iter) {
        m_iter = KIN_FUNCTIONAL;
    } else {
        throw CanteraError("KIN_Solver::setIterator", "unknown iterator");
    }
}
*/

/*
bool checkFlag(const int constraintFlag)
 {
     auto cflag = constraintFlag;
     bool valid_cflag = false;
     if (cflag == c_NONE || cflag == c_GE_ZERO || cflag == c_GT_ZERO ||
         cflag == c_LE_ZERO || cflag == c_LT_ZERO)
         valid_cflag = true;
     return valid_cflag;
 }
 */
 

/* The function KINSetConstraints is available only in latest SUNDIALS package
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
 */

/*
void KIN_Solver::sensInit(double t0, FuncEval& func)
{
    m_np = func.nparams();
    m_sens_ok = false;

    N_Vector y = N_VNew_Serial(static_cast<sd_size_t>(func.neq()));
    m_yS = N_VCloneVectorArray_Serial(static_cast<sd_size_t>(m_np), y);
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_yS[n]);
    }
    N_VDestroy_Serial(y);

    int flag = CVodeSensInit(m_kin_mem, static_cast<sd_size_t>(m_np),
                             CV_STAGGERED, CVSensRhsFn(0), m_yS);

    if (flag != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::sensInit", "Error in CVodeSensInit");
    }
    vector_fp atol(m_np);
    for (size_t n = 0; n < m_np; n++) {
        // This scaling factor is tuned so that reaction and species enthalpy
        // sensitivities can be computed simultaneously with the same abstol.
        atol[n] = m_abstolsens / func.m_paramScales[n];
    }
    flag = CVodeSensSStolerances(m_kin_mem, m_reltolsens, atol.data());
}
*/

void KIN_Solver::initialize(FuncEval& func)
{
    m_neq = func.neq();
    //m_t0 = t0;
    //m_time = t0;
    m_func = &func;
    func.clearErrors();

    if (m_y) {
        N_VDestroy_Serial(m_y); // free solution vector if already allocated
    }
    m_y = N_VNew_Serial(static_cast<sd_size_t>(m_neq)); // allocate solution vector
    N_VConst(0.0, m_y);      //TODO: Load this with initial values

    if (m_scale) {
        N_VDestroy_Serial(m_scale); // free solution vector if already allocated
    }
    m_scale = N_VNew_Serial(static_cast<sd_size_t>(m_neq)); // allocate solution vector
    N_VConst(1.0, m_scale);
    // check abs tolerance array size
    /*
    if (m_itol == KIN_SV && m_nabs < m_neq) {
        throw CanteraError("KIN_Solver::initialize",
                           "not enough absolute tolerance values specified.");
    }
    */
    
    if (m_constraints) {
        N_VDestroy_Serial(m_constraints);
    }
    m_constraints = N_VNew_Serial(static_cast<sd_size_t>(m_neq)); // allocate solution vector
    N_VConst(0.0, m_constraints);

    func.getState(NV_DATA_S(m_y));

    if (m_kin_mem) {
        KINFree(&m_kin_mem);
    }

    m_kin_mem = KINCreate();
    if (!m_kin_mem) {
        throw CanteraError("KIN_Solver::initialize",
                           "KINCreate failed.");
    }

    int flag = KINInit(m_kin_mem, kin_rhs, m_y);
    if (flag != KIN_SUCCESS) {
        if (flag == KIN_MEM_FAIL) {
            throw CanteraError("KIN_Solver::initialize",
                               "Memory allocation failed.");
        } else if (flag == KIN_ILL_INPUT) {
            throw CanteraError("KIN_Solver::initialize",
                               "Illegal value for CVodeInit input argument.");
        } else {
            throw CanteraError("KIN_Solver::initialize",
                               "KINInit failed.");
        }
    }
    KINSetErrHandlerFn(m_kin_mem, &kin_err, this);

    /*
    if (m_itol == KIN_SV) {
        flag = CVodeSVtolerances(m_kin_mem, m_reltol, m_abstol);
    } else {
        flag = CVodeSStolerances(m_kin_mem, m_reltol, m_abstols);
    }
    if (flag != KIN_SUCCESS) {
        if (flag == CV_MEM_FAIL) {
            throw CanteraError("KIN_Solver::initialize",
                               "Memory allocation failed.");
        } else if (flag == CV_ILL_INPUT) {
            throw CanteraError("KIN_Solver::initialize",
                               "Illegal value for CVodeInit input argument.");
        } else {
            throw CanteraError("KIN_Solver::initialize",
                               "CVodeInit failed.");
        }
    }
    */

    flag = KINSetUserData(m_kin_mem, &func);
    if (flag != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::initialize",
                           "KINSetUserData failed.");
    }
    /*
    if (func.nparams() > 0) {
        sensInit(t0, func);
        flag = KINSetSensParams(m_kin_mem, func.m_sens_params.data(),
                                  func.m_paramScales.data(), NULL);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::initialize",
                               "KINSetSensParams failed.");
        }
    }
    */
    applyOptions();
}

/*
void KIN_Solver::reinitialize(FuncEval& func)
{
    
    //m_t0 = t0;
    //m_time = t0;
    func.getState(NV_DATA_S(m_y));
    m_func = &func;
    func.clearErrors();

    int result = KINReInit(m_kin_mem, m_y);
    if (result != KIN_SUCCESS) {
        throw CanteraError("KIN_Solver::reinitialize",
                           "CVodeReInit failed. result = {}", result);
    }
    applyOptions();
}
*/

void KIN_Solver::applyOptions()
{
    if (m_type == DENSE + NOJAC) {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        #if CT_SUNDIALS_VERSION >= 30
            SUNLinSolFree((SUNLinearSolver) m_linsol);
            SUNMatDestroy((SUNMatrix) m_linsol_matrix);
            m_linsol_matrix = SUNDenseMatrix(N, N);
            if (m_linsol_matrix == nullptr) {
                throw CanteraError("KIN_Solver::applyOptions",
                    "Unable to create SUNDenseMatrix of size {0} x {0}", N);
            }
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLapackDense(m_y, (SUNMatrix) m_linsol_matrix);
            #else
                m_linsol = SUNDenseLinearSolver(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
            KINDlsSetLinearSolver(m_kin_mem, (SUNLinearSolver) m_linsol,
                                 (SUNMatrix) m_linsol_matrix);
        #else
            #if CT_SUNDIALS_USE_LAPACK
                KINLapackDense(m_kin_mem, N);
            #else
                KINDense(m_kin_mem, N);
            #endif
        #endif
                /*
    } else if (m_type == DIAG) {
        KINDiag(m_kin_mem);
        */
    } else if (m_type == GMRES) {
        #if CT_SUNDIALS_VERSION >= 30
            m_linsol = SUNSPGMR(m_y, PREC_NONE, 0);
            KINSpilsSetLinearSolver(m_kin_mem, (SUNLinearSolver) m_linsol);
        #else
            KINSpgmr(m_kin_mem, PREC_NONE, 0);
        #endif
    } else if (m_type == BAND + NOJAC) {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        long int nu = m_mupper;
        long int nl = m_mlower;
        #if CT_SUNDIALS_VERSION >= 30
            SUNLinSolFree((SUNLinearSolver) m_linsol);
            SUNMatDestroy((SUNMatrix) m_linsol_matrix);
            m_linsol_matrix = SUNBandMatrix(N, nu, nl, nu+nl);
            if (m_linsol_matrix == nullptr) {
                throw CanteraError("KIN_Solver::applyOptions",
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
            #if CT_SUNDIALS_USE_LAPACK
                KINLapackBand(m_kin_mem, N, nu, nl);
            #else
                KINBand(m_kin_mem, N, nu, nl);
            #endif
        #endif
    } else {
        throw CanteraError("KIN_Solver::applyOptions",
                           "unsupported option");
    }

    /*
    if (m_maxord > 0) {
        KINSetMaxOrd(m_kin_mem, m_maxord);
    }
    */
    if (m_maxsteps > 0) {
        KINSetNumMaxIters(m_kin_mem, m_maxsteps);
    }
    if (m_hmax > 0) {
        KINSetMaxNewtonStep(m_kin_mem, m_hmax);
    }
    /*
    if (m_hmin > 0) {
        KINSetMinStep(m_kin_mem, m_hmin);
    }
    */
    /*
    if (m_maxErrTestFails > 0) {
        KINSetMaxErrTestFails(m_kin_mem, m_maxErrTestFails);
    }
    */
}

void KIN_Solver::solve()
{
    int flag = KINSol(m_kin_mem, m_y, KIN_LINESEARCH, m_scale, m_scale);
    if (flag != KIN_SUCCESS) {
        string f_errs = m_func->getErrors();
        if (!f_errs.empty()) {
            f_errs = "Exceptions caught during RHS evaluation:\n" + f_errs;
        }
        throw CanteraError("KIN_Solver::solve",
            "KINSol error encountered. Error code: {}\n{}\n{}\n",
            flag, m_error_message, f_errs);
        /*
        throw CanteraError("KIN_Solver::solve",
            "KINSol error encountered. Error code: {}\n{}\n"
            "{}"
            "Components with largest weighted error estimates:\n{}",
            flag, m_error_message, f_errs, getErrorInfo(10));
        */
    }
}

/*
double KIN_Solver::step(double tout)
{
    int flag = CVode(m_kin_mem, tout, m_y, &m_time, CV_ONE_STEP);
    if (flag != KIN_SUCCESS) {
        string f_errs = m_func->getErrors();
        if (!f_errs.empty()) {
            f_errs = "Exceptions caught during RHS evaluation:\n" + f_errs;
        }
        throw CanteraError("KIN_Solver::step",
            "CVodes error encountered. Error code: {}\n{}\n"
            "{}"
            "Components with largest weighted error estimates:\n{}",
            flag, f_errs, m_error_message, getErrorInfo(10));

    }
    m_sens_ok = false;
    return m_time;
}
*/

int KIN_Solver::nEvals() const
{
    long int ne;
    KINGetNumFuncEvals(m_kin_mem, &ne);
    return ne;
}

/*
double KIN_Solver::sensitivity(size_t k, size_t p)
{
    if (m_time == m_t0) {
        // calls to CVodeGetSens are only allowed after a successful time step.
        return 0.0;
    }
    if (!m_sens_ok && m_np) {
        int flag = CVodeGetSens(m_kin_mem, &m_time, m_yS);
        if (flag != KIN_SUCCESS) {
            throw CanteraError("KIN_Solver::sensitivity",
                               "CVodeGetSens failed. Error code: {}", flag);
        }
        m_sens_ok = true;
    }

    if (k >= m_neq) {
        throw CanteraError("KIN_Solver::sensitivity",
                           "sensitivity: k out of range ({})", k);
    }
    if (p >= m_np) {
        throw CanteraError("KIN_Solver::sensitivity",
                           "sensitivity: p out of range ({})", p);
    }
    return NV_Ith_S(m_yS[p],k);
}
*/

/*
string KIN_Solver::getErrorInfo(int N)
{
    N_Vector errs = N_VNew_Serial(static_cast<sd_size_t>(m_neq));
    N_Vector errw = N_VNew_Serial(static_cast<sd_size_t>(m_neq));
    CVodeGetErrWeights(m_kin_mem, errw);
    CVodeGetEstLocalErrors(m_kin_mem, errs);

    vector<tuple<double, double, size_t> > weightedErrors;
    for (size_t i=0; i<m_neq; i++) {
        double err = NV_Ith_S(errs, i) * NV_Ith_S(errw, i);
        weightedErrors.emplace_back(-abs(err), err, i);
    }
    N_VDestroy(errs);
    N_VDestroy(errw);

    N = std::min(N, static_cast<int>(m_neq));
    sort(weightedErrors.begin(), weightedErrors.end());
    fmt::memory_buffer s;
    for (int i=0; i<N; i++) {
        format_to(s, "{}: {}\n",
                get<2>(weightedErrors[i]), get<1>(weightedErrors[i]));
    }
    return to_string(s);
}
*/

}
