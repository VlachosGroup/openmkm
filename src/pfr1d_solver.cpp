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

PFR1dSolver::PFR1dSolver(shared_ptr<PFR1d> pfr) : 
    m_solver(nullptr), m_ss_started(false),
    m_z(0.0), m_h0(0.0), m_rtol(0.0), m_atol(0.0),
    m_maxSteps(0), m_tStop(0.0), m_nquad(0),
    m_atolsens(0), m_rtolsens(0)
{
    m_neq = pfr->nEquations();
    m_vec.resize(m_neq);
    m_var = pfr->variablesNames();
    m_var_gas = pfr->gasVariableNames();
    m_var_state = pfr->stateVariableNames();
    m_var_surf = pfr->surfaceVariableNames();

    m_pfr = pfr;
    /*
    vector<int> constraints(pfr->nEquations());
    for (size_t i = 0; i < pfr->nEquations(); i++){
        constraints[i] = pfr->constraint(i);
        cout << i << " PFR1dSolver::init " << constraints[i] << endl;
    }
    setConstraints(constraints);
    */

}

void PFR1dSolver::init()
{
    try
    {
        m_solver = new IDA_Solver {*m_pfr};
        m_solver->setJacobianType(0);
        m_solver->setDenseLinearSolver();
        applyOptions();
        m_solver->setQuadratureVarSize(m_nquad);
        m_solver->init(0.0);
    }
    catch (Cantera::CanteraError& err)
    {
        std::cerr << err.what() << std::endl;
    }
    m_ss_gas.str(std::string());
    m_ss_surf.str(std::string());
    m_ss_state.str(std::string());
    m_ss_state.precision(6);
    m_ss_gas.precision(6);
    m_ss_surf.precision(6);
    m_ss_started = false;
}

void PFR1dSolver::reinit()
{
    if (m_solver != nullptr){
        delete m_solver;
    }
    init();
}

void PFR1dSolver::applyOptions()
{
    if (m_rtol) // atol and rtol are supplied in pair
        m_solver->setTolerances(m_rtol, m_atol);
    if (m_maxSteps)
        m_solver->setMaxNumSteps(m_maxSteps);
    if (m_h0)
        m_solver->setInitialStepSize(m_h0);
    if (m_tStop)
        m_solver->setStopTime(m_tStop);
    if (m_rtolsens)
        m_solver->setSensitivityTolerances(m_rtolsens, m_atolsens);
}

PFR1dSolver::~PFR1dSolver()
{
    if (m_solver != nullptr) delete m_solver;
}

void PFR1dSolver::setTolerances(double rtol, double atol)
{
    m_rtol = rtol;
    m_atol = atol;
}

void PFR1dSolver::setSensitivityTolerances(double rtol, double atol)
{
    m_rtolsens = rtol;
    m_atolsens = atol;
}

void PFR1dSolver::setMaxNumSteps(unsigned maxSteps)
{
    m_maxSteps = maxSteps;
}

void PFR1dSolver::setInitialStepSize(double h0)
{
    m_h0 = h0;
}

void PFR1dSolver::setStopPosition(double tStop)
{
    //m_solver->setStopTime(tstop);
    m_tStop = tStop;
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
        //cout << x << endl;
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
    } // TODO get pointer to sol instead.
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

void PFR1dSolver::writeSensitivityData(const string & saveas, 
                                       const vector<string>& rxnids, 
                                       string sep)
{
    std::ofstream ofs(saveas, std::ios::out | std::ios::binary);
    if (!ofs)
    {
        throw std::runtime_error("Cannot output to file");
    }

    stringstream sensi_strm;
    sensi_strm.str(std::string());
    sensi_strm.precision(6);
    sensi_strm  << "Rxn-Species-ids";
    for (auto var : m_var_state) { 
	sensi_strm << sep << var; 
    }
    for (auto var : m_var_gas) { 
	sensi_strm << sep << var; 
    }
    for (auto var : m_var_surf) { 
	sensi_strm << sep << var; 
    }
    sensi_strm << endl;

    for (size_t j = 0; j < rxnids.size(); j++){
        sensi_strm << rxnids[j] << scientific ;
        //cout << x << endl;
        for (unsigned i = 0; i < neq(); ++i) {
            // The sensitivity is multiplied by 2 to conform with the manual evaluation of the 
            // sensitivity coefficients. The difference in 2 is attributable to the use of 
            // central difference scheme. Central difference scheme was found to be not resulting
            // in accurate derivatives in some cases, where the -ve perturbation is essentially ignored
            // because some variables & parameters can't be < 0.  
            sensi_strm << sep << m_solver->sensitivity(i, j)/m_solver->solution(i);
        }
        sensi_strm << std::endl;
    }

    ofs << sensi_strm.str();
}


void PFR1dSolver::writeFisherInformationMatrixDiag(const string & saveas, 
                                                   const vector<std::string>& rxnids, 
                                                   std::string sep)
{
    std::ofstream ofs(saveas, std::ios::out | std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Cannot output to file");
    }

    stringstream sensi_strm;
    sensi_strm.str(std::string());
    sensi_strm.precision(6);
    sensi_strm  << "Rxnid";
    /*
    for (auto var : m_var_state) { 
	sensi_strm << sep << var; 
    }
    for (auto var : m_var_gas) { 
	sensi_strm << sep << var; 
    }
    for (auto var : m_var_surf) { 
	sensi_strm << sep << var; 
    }
    */
    sensi_strm << sep << "FIM_Diag";
    sensi_strm << endl;

    auto fim = m_solver->quadratureVector();
    auto tm = m_solver->getCurrentTimeFromIDA();
    if (tm == 0.0) {
        throw std::runtime_error("Time should be greater than 0 to output Fisher Information Matrix");
    }

    for (size_t j = 0; j < rxnids.size(); j++){
        sensi_strm << rxnids[j] << sep << scientific << fim[j]/tm << endl;
    }
    /*
    for (size_t j = 0; j < rxnids.size(); j++){
        sensi_strm << rxnids[j] << scientific ;
        //cout << x << endl;
        for (unsigned i = 0; i < neq(); ++i) {
            // The sensitivity is multiplied by 2 to conform with the manual evaluation of the 
            // sensitivity coefficients. The difference in 2 is attributable to the use of 
            // central difference scheme. Central difference scheme was found to be not resulting
            // in accurate derivatives in some cases, where the -ve perturbation is essentially ignored
            // because some variables & parameters can't be < 0.  
            sensi_strm << sep << m_solver->sensitivity(i, j)/m_solver->solution(i);
        }
        sensi_strm << std::endl;
    }
    */
    ofs << sensi_strm.str();
}

} 

