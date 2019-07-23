/**
 * @file pfr1d.h
 */

// Provides a 1d plug-flow reactor model.

// This file is part of OpenMKM. See License.txt in the top-level directory 
// for license and copyright information.

#ifndef OMKM_PFR_1D_H
#define OMKM_PFR_1D_H

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include <boost/math/interpolators/barycentric_rational.hpp>

#include "cantera/IdealGasMix.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/SurfLatIntPhase.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/numerics/ResidJacEval.h"
#include "cantera/transport.h"

#ifndef SUPPRESS_WARNINGS
#define SUPPRESS_WARNINGS true
#endif

namespace OpenMKM
{

// Helper functions
static inline double circleArea(doublereal Di)
{
    return Cantera::Pi * Di * Di / 4;
}

static inline double sccmTocmps(doublereal sccm)
{
    return sccm / 60000000;
}

//! A plug flow reactor (PFR) model implemented in 1d. 
//!
/*!The model calculates the 
 * steady state conditions of the PFR as a function of z (axial direction).
 * To evaluate the steady state conditions, first a pseudo steady state is solved
 * for the surfaces at the inlet for given T and P conditions. The resulting 
 * state of PFR at the inlet is used to propagate the state as a function of z by
 * solving differential algebraic governing equations of the PFR.
 */
class PFR1d : public Cantera::ResidJacEval
{
public:
    //! Constructor
    /*!
     * @param gas Gas phase (IdealGasMix object) containing both thermo 
     *            properties and kinetics manager to initialize the class.
     * @param surf_kins Kinetics managers of the surfaces. These have to 
     *                  instances of either Interface or InterfaceInteractions 
     *                  classes. Supply them as vector<InterfaceKinetics*>
     * @param surf_phases Surface phases corresponding to the interface kinetics
     *                    managers. These have to instances of either Interface 
     *                    or InterfaceInteractions classes. Supply the same 
     *                    objects used for surf_kins. Supply them as 
     *                    vector<SurfPhase*>.
     */
    PFR1d(Cantera::IdealGasMix *gas, 
          std::vector<Cantera::InterfaceKinetics*> surf_kins,       
          std::vector<Cantera::SurfPhase*> surf_phases,
          double pfr_crosssection_area,
          double cat_abyv, 
          double inlet_gas_flowrate);

    ~PFR1d()
    {
        Cantera::appdelete();
        m_T_interp = nullptr;
    }

    void reinit();

    virtual int getInitialConditions(const double t0,
                                     double *const y,
                                     double *const ydot);

    //! Evalueate the residual functional F(t, y, y') = 0 of differential 
    //! algebraic equations corresponding to 1d PFR.
    /*!
     * @param t z value from the PFR entrance. In DAE parlance, this is the
     *          individual variable and is typically denoted as t, because 
     *          often time is the independent variable. 
     * @param delta_t Step in z to evaluate the jacobian
     * @param y State of the PFR. The state consists of gas velocity, density,
     *          pressure, gas mass fractions and coverages of the surfaces.
     * @param ydot First order derivates of the state variables w.r.t. z.
     * @param resid Residual function F(t, y, y') corresponding to the DAEs 
     *              of PFR
     * @param evalType
     */
    virtual int evalResidNJ(const double t,
                    const double delta_t,
                    const double* const y,
                    const double* const ydot,
                    double* const resid,
                    const Cantera::ResidEval_Type_Enum evalType = Cantera::Base_ResidEval,
                    const int id_x = -1,
                    const double delta_x = 0.0);

    //! Evaluate the production rates of the species at the surfaces
    double evalSurfaces();

    unsigned getSpeciesIndex(std::string name) const
    {
        return m_gas->kineticsSpeciesIndex(name);
    }

    void setConstraints();

    double getIntEnergyMass() const
    {
        return m_gas->intEnergy_mass();
    }

    std::vector<std::string> variablesNames() const
    {
        return m_var;
    }

    std::vector<std::string> stateVariableNames() const
    {
        std::vector<std::string> state_var;
        for (size_t i = 0; i < m_neqs_extra; i++)
            state_var.push_back(m_var[i]);
        return state_var;
    }

    std::vector<std::string> gasVariableNames() const
    {
        std::vector<std::string> gas_var;
        for (size_t i = 0; i < m_nsp; i++){
            auto k = i + m_neqs_extra;
            gas_var.push_back(m_var[k]);
        }
        return gas_var;
    }

    std::vector<std::string> surfaceVariableNames()
    {
        size_t nsurf = neq_ - m_neqs_extra - m_nsp;
        std::vector<std::string> surf_var;
        for (size_t i = 0; i < nsurf; i++){
            auto k = i + m_neqs_extra + m_nsp;
            surf_var.push_back(m_var[k]);
        }
        return surf_var;
    }

    void setFlowRate(double flow_rate) 
    {
        m_u0 = flow_rate;
    }

    void getSurfaceInitialConditions(double* y);

    void getSurfaceProductionRates(double* y);

    void setEnergy(int eflag) 
    {
        if (eflag > 0) {
            m_energy = true;
            m_neqs_extra = 4;
        } else {
            m_energy = false;
            m_neqs_extra = 3;
        }
    }

    bool energyEnabled() const 
    {
        return m_energy;
    }

    void setHeatTransfer(double htc, double Text, double wall_abyv);

    double getHeat(double Tint) const;

    //! User supplied T profile in the PFR.
    //! The profile is supplied as map of distance and temperature values.
    //! If the first z is 0.0, then T has to be equal to m_T0 at inlet.
    //! The code doesn't do the check and if the condition is not hold,
    //! the behavior is undefined.
    void setTProfile(const std::map<double, double>& T_profile);

    //! Applicable only for cases where energy equation is not solved.
    double getT(double z);

    Cantera::thermo_t& contents() 
    {
        if (!m_gas) {
            throw Cantera::CanteraError("PFR1d::contents",
                               "Reactor contents not defined.");
        }
        return *m_gas;
    }

    const Cantera::thermo_t& contents() const 
    {
        if (!m_gas) {
            throw Cantera::CanteraError("PFR1d::contents",
                               "Reactor contents not defined.");
        }
        return *m_gas;
    }

    //! Return a reference to the *n*-th SurfPhase connected to the PFR
    Cantera::SurfPhase* surface(size_t n) 
    {
        return m_surf_phases[n];
    }

protected:
    //! Pointer to the gas phase object.
    Cantera::IdealGasMix *m_gas = nullptr;

    //! Surface kinetics objects
    std::vector<Cantera::InterfaceKinetics*> m_surf_kins;

    //! Surface phases objects.
    //! Both surface kinetics and phases have to refer
    //! to same objects, which are instances of ether
    //! Interface or InterfaceInteractions classes
    std::vector<Cantera::SurfPhase*> m_surf_phases;

    //! Species molar weights.
    std::vector<doublereal> m_W;

    //! Species net production rates in  gas phase reactions.
    std::vector<double> m_wdot;

    //! Species net production rates in surface reactions.
    std::vector<double> m_sdot;

    //! Species names and variables.
    std::vector<std::string> m_var;

    //! Number of equations in residual.
    unsigned m_neq;

    //! Number of gas species.
    unsigned m_nsp;

    //! Number of extra equations to solve beyond gas and surface species number
    unsigned m_neqs_extra; 

    //! Catalyst area by volume
    double m_cat_abyv;

    //! Reference state inlet density.
    double m_rho_ref;

    //! Reactor cross-sectional area
    const double m_Ac = 0.0;

    //! Solve Energy Equation
    bool m_energy;

    //! Inlet flow rate
    double m_u0 = 0.0;

    //! Inlet temperature
    double m_T0 = 0.0;

    //! Tprofile
    // The idea is to use m_T_profile_iind to keep track of current index of 
    // T_profile. For distance between m_T_profile_ind[m_T_profile_iind] and 
    // m_T_profile_ind[m_T_profile_iind+1], T is given as linear interpolation 
    // of m_T_profile[m_T_profile_iind] and m_T_profile[m_T_profile_iind+1].
    // If m_T_profile_ind[0] > 0.0, then m_T_profile_iind=-1 and T is given as
    // linear interpolation of m_T0 and m_T_profile[0]
    std::vector<double> m_T_profile;

    //! Index of Tprofile in terms of distance
    std::vector<double> m_T_profile_ind;
    
    //! Current index of Tprofile_ind
    int m_T_profile_iind;

    //! External Heat supplied 
    bool m_heat;

    //! External temperature
    double m_Text = 0;

    //! Heat transfer coefficent
    double m_htc = 0;

    //! Reactor Wall area where heat is transferred
    double m_surf_ext_abyv = 0;

    //! Inlet pressure
    double m_P0 = 0;

    std::shared_ptr<boost::math::barycentric_rational<double>> m_T_interp;

}; 

} 
#endif 
