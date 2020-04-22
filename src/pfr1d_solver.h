//! @file pfr1d_solver.h

// Solves the DAE numerical equations w.r.t z using IDAS

// This file is part of OpenMKM. See License.txt in the top-level directory 
// for license and copyright information.

#ifndef OMKM_SOLVER_1D_H
#define OMKM_SOLVER_1D_H

#include <vector>
#include <memory>
#include "cantera/numerics/IDA_Solver.h"
#include "pfr1d.h"

namespace OpenMKM {

class PFR1dSolver
{
public:
    PFR1dSolver(std::shared_ptr<PFR1d> pfr);

    ~PFR1dSolver();

    void init();

    void reinit();
    
    void setTolerances(double rtol, double atol);

    void setSensitivityTolerances(double rtol, double atol);

    void setQuadratureSize(size_t nquad) {
        m_nquad = nquad;
    }

    void setMaxNumSteps(unsigned maxsteps);

    void setMaxTimeStep(double maxtimestep);

    void setInitialStepSize(double h0);

    void setStopPosition(double tstop);

    //void setConstraints(const std::vector<int> constraints);

    int solve(double xout);

    double z() {
        return m_z;
    }

    //! Relative tolerance
    double rtol(){
        return m_rtol;
    }

    //! Absolute tolerance
    double atol(){
        return m_atol;
    }

    //! Relative sensitivity tolerance
    doublereal rtolSensitivity() const {
        return m_rtolsens;
    }

    //! Absolute sensitivity tolerance
    doublereal atolSensitivity() const {
        return m_atolsens;
    }

    double solution(unsigned num) const;

    std::vector<double> solutionVector();

    double derivative(unsigned num) const;

    std::vector<double> derivativeVector();

    std::vector<std::string> variablesNames() const;

    void writeStateData(const std::string & saveas);

    void writeGasData(const std::string & saveas);

    void writeSurfaceData(const std::string & saveas);

    void writeSensitivityData(const std::string & saveas, 
                              const std::vector<std::string> & rxnids, 
                              std::string sep);

    void writeFisherInformationMatrixDiag(const std::string & saveas, 
                                          const std::vector<std::string> & rxnids, 
                                          std::string sep);

    Cantera::DAE_Solver& solver() {
        return *m_solver;
    }

    size_t neq(){
        return m_neq;
    }

protected:
    void applyOptions();

    //! To provide access to results.
    std::vector<double> m_vec;

    //! Provides access to variables names.
    std::vector<std::string> m_var;
    std::vector<std::string> m_var_gas;
    std::vector<std::string> m_var_surf;
    std::vector<std::string> m_var_state;

    //! Pointer to IDA solver.
    Cantera::IDA_Solver* m_solver;

    //! Buffers for results.
    std::stringstream m_ss_gas;
    std::stringstream m_ss_surf;
    std::stringstream m_ss_state;

    //! Check if writer was started.
    bool m_ss_started = false;

    //! Current z position during advancing
    double m_z;

    //! Number of equations in reactor model.
    unsigned m_neq;

    //! Number of quadrature equations to be computed.
    unsigned m_nquad;

    //! Absolute DAE solver tolerance
    double m_atol;

    //! Relative DAE solver tolerance
    double m_rtol;

    //! Absolute DAE solver tolerance for sensitivity
    double m_atolsens;

    //! Relative DAE solver tolerance for sensitivity
    double m_rtolsens;

    //! Initial Step Size
    double m_h0;

    //! Max solver internal steps 
    size_t m_maxSteps;

    //! Simulation end time
    double m_tStop;

    //! PFR object
    std::shared_ptr<PFR1d> m_pfr;
}; 

} 

#endif 
