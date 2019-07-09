/****************************************************************************
 * Solves the DAE numerical equations w.r.t z using IDAS
***************************************************************************
*/

#include "pfr1d_solver.h"
#include "pfr1d.h"
#include "cantera/numerics/IDA_Solver.h"

using namespace Cantera;
using namespace std;

namespace OpenMKM 
{

PFR1dSolver::PFR1dSolver(PFR1d* pfr)
{
    m_neq = pfr->nEquations();
    m_vec.resize(m_neq);
    m_var = pfr->variablesNames();
    m_var_gas = pfr->gasVariableNames();
    m_var_state = pfr->stateVariableNames();
    m_var_surf = pfr->surfaceVariableNames();

    m_ss_state.precision(6);
    m_ss_gas.precision(6);
    m_ss_surf.precision(6);
    try
    {
        m_solver = new IDA_Solver {*pfr};
        m_solver->setJacobianType(0);
        m_solver->setDenseLinearSolver();
        m_solver->init(0.0);
    }
    catch (Cantera::CanteraError& err)
    {
        std::cerr << err.what() << std::endl;
    }

    /*
    vector<int> constraints(pfr->nEquations());
    for (size_t i = 0; i < pfr->nEquations(); i++){
        constraints[i] = pfr->constraint(i);
        cout << i << " PFR1dSolver::init " << constraints[i] << endl;
    }
    setConstraints(constraints);
    */

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

/*
void PFR1dSolver::setConstraints(const vector<int> constraints){
    m_solver->setConstraints(constraints.data());
}*/

int PFR1dSolver::solve(double xout)
{ 
    int retcode = 0;

    size_t nstate = m_var_state.size();
    size_t ngas = m_var_gas.size();
    size_t nsurf = m_var_surf.size();

    auto print_state = [&](double x, string sep){

        m_ss_state << scientific << x;;
        m_ss_gas << scientific << x;
        m_ss_surf << scientific << x;
        cout << x << endl;
        for (unsigned i = 0; i != nstate; ++i) {
            m_ss_state << sep << m_solver->solution(i);
        }
        m_ss_state << std::endl;
        for (unsigned i = 0; i != ngas; ++i) {
            m_ss_gas << sep << m_solver->solution(i + nstate);
        }
        m_ss_gas << std::endl;
        for (unsigned i = 0; i != nsurf; ++i) {
            m_ss_surf << sep << m_solver->solution(i + nstate + ngas);
        }
        m_ss_surf << std::endl;

    };

    if (!m_ss_started)
    {
        auto print_header = [&](string sep) -> void {
            m_ss_gas  << "z(m)";
            m_ss_surf << "z(m)";
            m_ss_state  << "z(m)";
            for (auto var : m_var_gas) { 
                m_ss_gas << sep << var; 
            }
            for (auto var : m_var_surf) { 
                m_ss_surf << sep << var; 
            }
            for (auto var : m_var_state) { 
                m_ss_state << sep << var; 
            }
            m_ss_gas << endl;
            m_ss_surf << endl;
            m_ss_state << endl;
        };
        print_header(",");
        print_state(0.0, ",");
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
    print_state(xout, ",");
    
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

double PFR1dSolver::derivative(unsigned num) const
{
    return m_solver->derivative(num);
}

vector<double> PFR1dSolver::derivativeVector()
{
    // TODO make this with STL algorithm.
    const double* der = m_solver->derivativeVector();
    for (unsigned i = 0; i != m_vec.size(); ++i) {
        m_vec[i] = der[i];
    }
    return m_vec;
}

vector<string> PFR1dSolver::variablesNames() const
{
    return m_var;
}

void PFR1dSolver::writeStateData(const string & saveas)
{
    std::ofstream ofs(saveas, std::ios::out | std::ios::binary);
    if (!ofs)
    {
        throw std::runtime_error("Cannot output to file");
    }
    ofs << m_ss_state.str();
    ofs.close();
}

void PFR1dSolver::writeGasData(const string & saveas)
{
    std::ofstream ofs(saveas, std::ios::out | std::ios::binary);
    if (!ofs)
    {
        throw std::runtime_error("Cannot output to file");
    }
    ofs << m_ss_gas.str();
    ofs.close();
}

void PFR1dSolver::writeSurfaceData(const string & saveas)
{
    std::ofstream ofs(saveas, std::ios::out | std::ios::binary);
    if (!ofs)
    {
        throw std::runtime_error("Cannot output to file");
    }
    ofs << m_ss_surf.str();
    ofs.close();
}

} 

