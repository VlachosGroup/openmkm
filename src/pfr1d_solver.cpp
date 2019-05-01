/****************************************************************************
 * Solves the DAE numerical equations w.r.t z using IDAS
***************************************************************************
*/

#include "pfr1d_solver.h"
#include "pfr1d.h"
#include "cantera/numerics/IDA_Solver.h"

using namespace Cantera;
using namespace std;

namespace HeteroCt 
{

PFR1dSolver::PFR1dSolver(PFR1d* pfr)
{
    m_neq = pfr->nEquations();
    m_vec.resize(m_neq);
    m_var = pfr->variablesNames();

    try
    {
        cout << "Before solver initialized" << endl;
        m_solver = new IDA_Solver {*pfr};
        m_solver->setJacobianType(0);
        m_solver->setDenseLinearSolver();
        m_solver->init(0.0);
        cout << "Solver initialized" << endl;
    }
    catch (Cantera::CanteraError& err)
    {
        std::cerr << err.what() << std::endl;
    }
}

PFR1dSolver::~PFR1dSolver()
{
    if (m_solver != nullptr) delete m_solver;
}

void PFR1dSolver::setTolerances(double rtol, double atol)
{
    m_solver->setTolerances(rtol, atol);
}

void PFR1dSolver::setMaxNumSteps(unsigned maxsteps)
{
    m_solver->setMaxNumSteps(maxsteps);
}

void PFR1dSolver::setInitialStepSize(double h0)
{
    m_solver->setInitialStepSize(h0);
}

void PFR1dSolver::setStopPosition(double tstop)
{
    m_solver->setStopTime(tstop);
}

int PFR1dSolver::solve(double xout)
{ 
    int retcode = 0;

    if (!m_ss_started)
    {
        m_ss << "z (m),";
        for (auto var : m_var) { 
            m_ss << var << ","; 
        }
        m_ss << endl;

        m_ss << 0.0 << ","; 
        for (unsigned i = 0; i != m_neq; ++i) {
            m_ss << m_solver->solution(i) << ",";
        }
        m_ss << endl;

        m_ss_started = true;
    }

    // TODO Manage return codes from IDA_Solver.solve.
    try
    {
        retcode = m_solver->solve(xout);
    }
    catch (CanteraError& err)
    {
        std::cerr << err.what() << std::endl;
        retcode = -99;
    }

    // TODO get pointer to sol instead.
    m_ss << xout << ",";
    cout << xout << endl;
    for (unsigned i = 0; i != m_neq; ++i) {
        m_ss << m_solver->solution(i) << ",";
    }
    m_ss << std::endl;

    return retcode;
}

double PFR1dSolver::solution(unsigned num) const
{
    return m_solver->solution(num);
}

vector<double> PFR1dSolver::solutionVector()
{
    // TODO make this with STL algorithm.
    const double* sol = m_solver->solutionVector();
    for (unsigned i = 0; i != m_vec.size(); ++i) { 
        m_vec[i] = sol[i]; 
    }
    return m_vec;
}

vector<string> PFR1dSolver::variablesNames() const
{
    return m_var;
}

void PFR1dSolver::writeResults(const string & saveas)
{
    std::ofstream ofs(saveas, std::ios::out | std::ios::binary);
    if (!ofs)
    {
        throw std::runtime_error("Cannot output to file");
    }
    ofs << m_ss.str();
    ofs.close();
}

} 

