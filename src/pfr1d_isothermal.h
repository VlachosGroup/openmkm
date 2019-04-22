// ***************************************************************************
// Provides an isothermal plug-flow reactor model.
// ***************************************************************************

#ifndef HTRCT_PFR1D_ISOTHERMAL
#define HTRCT_PFR1D_ISOTHERMAL

#include "pfr1d_constA.h"

namespace HeteroCt
{

class PFR1dIsothermal : public PFR1dConstArea
{
public:
    PFR1dIsothermal(std::string const& mech,
                    std::string phase,
                    double Di,
                    double T0,
                    double p0,
                    std::string X0,
                    double Q0)
        : PFR1dConstArea{mech, phase, Di, T0, p0, X0, Q0, neqs_extra_},
          idx0{nspec_gas_+0},
          idx1{nspec_gas_+1},
          idx2{nspec_gas_+2},
          m_T0{T0}
    {
        m_var.push_back("u");
        m_var.push_back("rho");
        m_var.push_back("p");

        std::cout << "\nStarting solver : " << "IsothermalPFR"
                  << "\nInitial temperature (K) . " << m_gas->temperature()
                  << "\nInitial pressure (Pa) ... " << m_gas->pressure()
                  << "\nInitial velocity (m/s) .. " << m_u0
                  << "\nNumber of equations ..... " << neq_
                  << std::endl;
    }

    int getInitialConditions(const double t0,
                             double *const y,
                             double *const ydot);

    int evalResidNJ(const double t,
                    const double delta_t,
                    const double* const y,
                    const double* const ydot,
                    double* const resid,
                    const Cantera::ResidEval_Type_Enum evalType = Cantera::Base_ResidEval,
                    const int id_x = -1,
                    const double delta_x = 0.0);

private:
    //! Number of extra equations.
    static constexpr unsigned neqs_extra_ = 3;

    //! Index of extra equations.
    const unsigned idx0 = 0, idx1 = 0, idx2 = 0;

    //! Reactor temperature.
    const doublereal m_T0 = 0.0;

}; 

} 

#endif 
