// ***************************************************************************
// Provides a generic plug-flow reactor model.
// ***************************************************************************

#ifndef PFR_1D_BASE
#define PFR_1D_BASE

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include "cantera/IdealGasMix.h"
#include "cantera/numerics/IDA_Solver.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/transport.h"

#ifndef SUPPRESS_WARNINGS
#define SUPPRESS_WARNINGS true
#endif

namespace HeteroCt
{

static inline double circleArea(doublereal Di)
{
    return Cantera::Pi * Di * Di / 4;
}

static inline double sccmTocmps(doublereal sccm)
{
    return sccm / 60000000;
}

class PFR1dBase : public Cantera::ResidJacEval
{
public:
    PFR1dBase(std::string const& mech, std::string phase, double T0,
              double p0, std::string X0, unsigned neqs_extra_)
        : Cantera::ResidJacEval{} {

        Cantera::suppress_thermo_warnings(SUPPRESS_WARNINGS);

        std::cout << std::boolalpha
                  << "\nIntegrating PFR"
                  << "\nUsing Sundials : " << CT_SUNDIALS_VERSION
                  << "\nUsing LAPACK   : " << bool(CT_SUNDIALS_USE_LAPACK)
                  << std::endl;

        m_gas = new Cantera::IdealGasMix {mech, phase};
        m_gas->setState_TPX(273.15, Cantera::OneAtm, X0);
        m_rho_ref = m_gas->density();

        m_gas->setState_TPX(T0, p0, X0);

        nspec_gas_ = m_gas->nSpecies();
        neq_ = nspec_gas_ + neqs_extra_;

        m_W.resize(nspec_gas_);
        m_wdot.resize(nspec_gas_);
        m_gas->getMolecularWeights(m_W);

        m_var = m_gas->speciesNames();

        
    }

    ~PFR1dBase()
    {
        if (m_gas != nullptr) delete m_gas;
        Cantera::appdelete();
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

    unsigned getSpeciesIndex(std::string name) const
    {
        return m_gas->kineticsSpeciesIndex(name);
    }

    double getIntEnergyMass() const
    {
        return m_gas->intEnergy_mass();
    }


    std::vector<std::string> variablesNames() const
    {
        return m_var;
    }

protected:
    //! Pointer to the gas phase object.
    Cantera::IdealGasMix *m_gas = nullptr;

    //! Species molar weights.
    std::vector<doublereal> m_W;

    //! Species net production rates.
    std::vector<double> m_wdot;

    //! Species names and variables.
    std::vector<std::string> m_var;

    //! Number of gas phase species.
    unsigned nspec_gas_;

    //! Reference state inlet density.
    double m_rho_ref;

}; 

} 
#endif 
