/****************************************************************************
 * Solves the DAE numerical equations w.r.t z using IDAS
***************************************************************************
*/

#ifndef HTRCT_SOLVER_1D_H
#define HTRCT_SOLVER_1D_H

#include "pfr1d.h"
#include "cantera/numerics/IDA_Solver.h"

namespace HeteroCt {

class PFR1dSolver
{
public:
    PFR1dSolver(PFR1d* pfr);

    ~PFR1dSolver();

    void setTolerances(double rtol, double atol);

    void setMaxNumSteps(unsigned maxsteps);

    void setMaxTimeStep(double maxtimestep);

    void setInitialStepSize(double h0);

    void setStopPosition(double tstop);

    int solve(double xout);

    double z() {
        return m_z;
    }

    double rtol(){
        return m_rtol;
    }

    double atol(){
        return m_atol;
    }

    double solution(unsigned num) const;

    std::vector<double> solutionVector();

    double derivative(unsigned num) const;

    std::vector<double> derivativeVector();

    std::vector<std::string> variablesNames() const;

    void writeResults(const std::string & saveas);

    Cantera::DAE_Solver& solver() {
        return *m_solver;
    }

    size_t neq(){
        return m_neq;
    }

protected:
    //! To provide access to results.
    std::vector<double> m_vec;

    //! Provides access to variables names.
    std::vector<std::string> m_var;

    //! Pointer to IDA solver.
    Cantera::IDA_Solver* m_solver;

    //! Buffer for results file.
    std::stringstream m_ss;

    //! Check if writer was started.
    bool m_ss_started = false;

    //! Current z position during advancing
    double m_z;

    //! Number of equations in reactor model.
    unsigned m_neq;

    //! Absolute DAE solver tolerance
    double m_atol;

    //! Relative DAE solver tolerance
    double m_rtol;


}; 

} 

#endif 
