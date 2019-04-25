/****************************************************************************
 * Solves the DAE numerical equations w.r.t z using IDAS
***************************************************************************
*/

#ifndef HTRCT_SOLVER_1D
#define HTRCT_SOLVER_1D

#include "pfr1d.h"
#include "cantera/numerics/IDA_Solver.h"

using namespace Cantera;
using namespace std;

namespace HeteroCt {

class PFR1dSolver
{
public:
    PFR1dSolver(PFR1d* pfr)
    {
        m_neq = pfr->nEquations();
        m_vec.resize(m_neq);
        m_var = pfr->variablesNames();

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
    }

    ~PFR1dSolver()
    {
        if (m_solver != nullptr) delete m_solver;
    }

    void setTolerances(double rtol, double atol)
    {
        m_solver->setTolerances(rtol, atol);
    }

    void setMaxNumSteps(unsigned maxsteps)
    {
        m_solver->setMaxNumSteps(maxsteps);
    }

    void setInitialStepSize(double h0)
    {
        m_solver->setInitialStepSize(h0);
    }

    void setStopPosition(double tstop)
    {
        m_solver->setStopTime(tstop);
    }

    int solve(double xout)
    { 
        int retcode = 0;

        if (!m_ss_started)
        {
            for (auto var : m_var) { m_ss << var << ","; }
            m_ss << "x" << std::endl;

            for (unsigned i = 0; i != m_neq; ++i) {
                m_ss << m_solver->solution(i) << ",";
            }
            m_ss << 0.0 << std::endl;

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
        for (unsigned i = 0; i != m_neq; ++i) {
            m_ss << m_solver->solution(i) << ",";
        }
        m_ss << xout << std::endl;

        return retcode;
    }

    double solution(unsigned num) const
    {
        return m_solver->solution(num);
    }

    std::vector<double> solutionVector()
    {
        // TODO make this with STL algorithm.
        const double* sol = m_solver->solutionVector();
        for (unsigned i = 0; i != m_vec.size(); ++i) { m_vec[i] = sol[i]; }
        return m_vec;
    }

    std::vector<std::string> variablesNames() const
    {
        return m_var;
    }

    void writeResults(const std::string & saveas)
    {
        std::ofstream ofs(saveas, std::ios::out | std::ios::binary);
        if (!ofs)
        {
            throw std::runtime_error("Cannot output to file");
        }
        ofs << m_ss.str();
        ofs.close();
    }
protected:
    //! To provide access to results.
    std::vector<double> m_vec;

    //! Provides access to variables names.
    std::vector<std::string> m_var;

    //! Pointer to IDA solver.
    IDA_Solver* m_solver;

    //! Buffer for results file.
    std::stringstream m_ss;

    //! Check if writer was started.
    bool m_ss_started = false;

    //! Number of equations in reactor model.
    unsigned m_neq;

}; 

} 

#endif 
