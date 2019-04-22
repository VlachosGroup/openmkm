// ***************************************************************************
// Provides a generic constant area plug-flow reactor model.
// ***************************************************************************

#ifndef HTRCT_PFR1D_CONSTAREA
#define HTRCT_PFR1D_CONSTAREA

#include "pfr1d_base.h"

namespace HeteroCt
{

class PFR1dConstArea : public PFR1dBase
{
public:
    PFR1dConstArea(std::string const& mech, std::string phase, double Di,
                   double T0, double p0, std::string X0, double Q0,
                   unsigned neqs_extra_)
        : PFR1dBase{mech, phase, T0, p0, X0, neqs_extra_}, m_Ac{circleArea(Di)},
          m_u0{setVelocity(Q0)}
    {

    }

    virtual int getInitialConditions(const double t0,
                                     double *const y,
                                     double *const ydot) = 0;

    virtual int evalResidNJ(const double t,
                            const double delta_t,
                            const double* const y,
                            const double* const ydot,
                            double* const resid,
                            const Cantera::ResidEval_Type_Enum evalType = Cantera::Base_ResidEval,
                            const int id_x = -1,
                            const double delta_x = 0.0) = 0;

    //! Compute inlet velocity.
    virtual const double setVelocity(double Q0) const
    {
        return (m_rho_ref / m_gas->density()) * sccmTocmps(Q0) / m_Ac;
    }


protected:
    //! Reactor cross-sectional area.
    const double m_Ac = 0.0;

    //! Inlet velocity.
    const double m_u0 = 0.0;

}; 

} 

#endif 
