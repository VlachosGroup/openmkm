/**
 * @file pfr1d.h
 */
// Provides a 1d plug-flow reactor model.

// This file is part of hetero_ct

#ifndef PFR_1D_BASE_H
#define PFR_1D_BASE_H

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
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/SurfLatIntPhase.h"
#include "cantera/numerics/IDA_Solver.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/transport.h"

#ifndef SUPPRESS_WARNINGS
#define SUPPRESS_WARNINGS true
#endif

namespace HeteroCt
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

//! A plug flow reactor (PFR) model implemented in 1d. The model calculates the 
//! steady state conditions of the PFR as a function of z (axial direction).
/*!
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
          double inlet_gas_velocity);

    ~PFR1d()
    {
        Cantera::appdelete();
    }

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
    double evalSurfaces(const double* const ydot);

    double evalSurfaceInit();

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

    void setVelocity(double velocity) 
    {
        m_u0 = velocity;
    }

    void getSurfaceInitialConditions(double* y);

    void getSurfaceProductionRates(double* y);

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

    //! Number of extra equations for IDAs solver
    unsigned m_neqs_extra; 

    //! Catalyst area by volume
    double m_surf_sp_area;

    //! Reference state inlet density.
    double m_rho_ref;

    //! Reactor cross-sectional area
    const double m_Ac = 0.0;

    //! Inlet velocity
    double m_u0 = 0.0;

    //! Inlet temperature
    double m_T0 = 0.0;

    //! Inlet pressure
    double m_P0 = 0.0;

}; 

} 
#endif 
