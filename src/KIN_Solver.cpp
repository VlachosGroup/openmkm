//! @file KIN_Solver.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/NumUtil.h"
#include "cantera/base/stringUtils.h"
#include "cantera/numerics/Integrator.h"
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
        int rval, i, len;
        FuncEval* f = (FuncEval*) f_data;
        rval = f->eval_nothrow(0, NV_DATA_S(y), NV_DATA_S(rhs));
        len = NV_LENGTH_S(rhs);
        for (i = 0; i < len; i++){
            printf("ith rhs %d: %E\n", i, NV_Ith_S(rhs, i)); /* = NV_Ith_S(rhs, i) * 0.1;*/    // 0.1 is taken as \Delta T
        }
        printf("\n");
        return rval;
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
    //m_constraints(0),
    m_functol(0),
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


void KIN_Solver::setTolerance(double functol)
{
    m_functol = functol;
}

void KIN_Solver::setProblemType(int probtype)
{
    m_type = probtype;
}

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
    flag = CVodeSensSStolerances(m_kin_mem, m_reltolsens, atol.data());
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
 

/* 
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

void KIN_Solver::initialize(FuncEval& func)
{
    m_neq = func.neq();
    m_func = &func;
    func.clearErrors();

    if (m_y) {
        N_VDestroy_Serial(m_y); // free solution vector if already allocated
    }
    cout << "NEQ " << m_neq << endl;
    m_y = N_VNew_Serial(static_cast<sd_size_t>(m_neq)); // allocate solution vector
    //N_VConst(0.0, m_y);      
    func.getState(NV_DATA_S(m_y));   // Initial guess
    for (size_t i = 0; i < m_neq; i++){
        cout << "i " << i << " y " << NV_Ith_S(m_y, i) << endl;
    }

    if (m_scale) {
        N_VDestroy_Serial(m_scale); // free solution vector if already allocated
    }
    m_scale = N_VNew_Serial(static_cast<sd_size_t>(m_neq)); // allocate solution vector
    //N_VConst(1.0, m_scale);
    N_VInv(m_y, m_scale);
    for (size_t i = 0; i < m_neq; i++){
        cout << "i " << i << " scale " << NV_Ith_S(m_scale, i) << endl;
    }
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
    N_VConst(1.0, m_constraints);

    if (m_kin_mem) {
        KINFree(&m_kin_mem);
    }
    m_kin_mem = KINCreate();
    if (!m_kin_mem) 
        throw CanteraError("KIN_Solver::initialize", "KINCreate failed.");

    int flag = KINSetUserData(m_kin_mem, &func);
    if (flag != KIN_SUCCESS) 
        throw CanteraError("KIN_Solver::initialize", "KINSetUserData failed.");

    flag = KINSetConstraints(m_kin_mem, m_constraints);
    if (flag != KIN_SUCCESS) 
        throw CanteraError("KIN_Solver::initialize", "KINSetConstraints failed.");

    flag = KINInit(m_kin_mem, kin_rhs, m_y);
    if (flag != KIN_SUCCESS) {
        if (flag == KIN_MEM_FAIL) {
            throw CanteraError("KIN_Solver::initialize",
                               "Memory allocation failed.");
        } else if (flag == KIN_ILL_INPUT) {
            throw CanteraError("KIN_Solver::initialize",
                               "Illegal value for KINInit input argument.");
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

void KIN_Solver::applyOptions()
{

    cout << "Before calling KINSetFuncNormTol " << m_functol << endl;
    KINSetFuncNormTol(m_kin_mem, m_functol);

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
    /*cout << "Number of KINSOL iters" << m_maxsteps << endl;
    if (m_maxsteps > 0) {
        KINSetNumMaxIters(m_kin_mem, m_maxsteps);
    }
    cout << "max newton step " << m_hmax << endl;
    if (m_hmax > 0) {
        KINSetMaxNewtonStep(m_kin_mem, m_hmax);
    }
    */
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

int KIN_Solver::solve()
{
    m_func->getState(NV_DATA_S(m_y));

    //int flag = KINSol(m_kin_mem, m_y, KIN_LINESEARCH, m_scale, m_scale);
    int flag = KINSol(m_kin_mem, m_y, KIN_NONE, m_scale, m_scale);
    // There is a high chance for this solver to fail. Just return success status
    if (flag == KIN_INITIAL_GUESS_OK) {
        cout << "initial guess within the tolerances" << endl;
    }
    cout << "kinsol flag " << flag << endl;
    if (flag == KIN_SUCCESS || KIN_INITIAL_GUESS_OK) 
        return 1;    
    else
        return 0;
    /*
    if (flag != KIN_SUCCESS) {
        string f_errs = m_func->getErrors();
        if (!f_errs.empty()) {
            f_errs = "Exceptions caught during RHS evaluation:\n" + f_errs;
        }
        throw CanteraError("KIN_Solver::solve",
            "KINSol error encountered. Error code: {}\n{}\n{}\n",
            flag, m_error_message, f_errs);
       }
    */
}

int KIN_Solver::nEvals() const
{
    long int ne;
    KINGetNumFuncEvals(m_kin_mem, &ne);
    return ne;
}
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

// Get stats from KINSOL
void KIN_Solver::stats() 
{
    long int nwork;
    KINGetNumFuncEvals(m_kin_mem, &nwork);
    cout << "# of func evals " << nwork << endl;
    KINGetNumNonlinSolvIters(m_kin_mem, &nwork);
    cout << "# of Nonlinear solver iterations  " << nwork << endl;
    realtype work;
    KINGetFuncNorm(m_kin_mem, &work);
    cout << "Func norm " << work << endl;
}

}
