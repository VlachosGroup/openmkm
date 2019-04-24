//! @file pfr1d.cpp 
//! Provides an plug-flow reactor model in 1d.

// This file is part of Hetero_ct. 
// Source is adapted from CanteraPFR codebase.

#include "pfr1d.h"

using namespace std;
using namespace Cantera;

namespace HeteroCt 
{

PFR1d::PFR1d(IdealGasMix *gas, vector<InterfaceKinetics*> surf_kins,
             vector<SurfPhase*> surf_phases, double area, 
             double inlet_gas_velocity) 
   : ResidJacEval{}, m_gas(gas), m_surf_kins(surf_kins), 
     m_surf_phases(surf_phases), m_Ac(area), m_u0(inlet_gas_velocity)
{
   
    suppress_thermo_warnings(SUPPRESS_WARNINGS);

    cout << boolalpha
             << "\nIntegrating PFR"
             << "\nUsing Sundials : " << CT_SUNDIALS_VERSION
             << "\nUsing LAPACK   : " << bool(CT_SUNDIALS_USE_LAPACK)
             << endl;

    m_rho_ref = m_gas->density();

    m_nsp = m_gas->nSpecies();
    neq_ = m_nsp + 3;
    for (const auto s_ph : m_surf_phases) {
        neq_ += s_ph->nSpecies();
    }

    m_W.resize(m_nsp);
    m_gas->getMolecularWeights(m_W);
    m_wdot.resize(m_nsp);
    fill(m_wdot.begin(), m_wdot.end(), 0.0);
    m_sdot.resize(m_nsp);
    fill(m_sdot.begin(), m_sdot.end(), 0.0);

    m_var = m_gas->speciesNames();

    m_T0 = gas->temperature();
    m_P0 = gas->pressure();
}   


int PFR1d::getInitialConditions(const double t0,
                                double *const y, 
                                double *const ydot)
{
    const double P0 = m_gas->pressure();
    const double rho0 = m_gas->density();
    const double Wavg = m_gas->meanMolecularWeight();
    const double RT = m_gas->temperature() * GasConstant;


    Eigen::MatrixXd A(m_nsp+3, m_nsp+3);
    Eigen::VectorXd b(m_nsp+3);

    y[0] = m_u0;
    y[1] = rho0;
    y[2] = P0;

    m_gas->getMassFractions(y+3);
    auto loc = m_nsp + 3;
    for (const auto s_ph : m_surf_phases) {
        s_ph->getCoverages(ydot+loc);
        loc += s_ph->nSpecies();
    }

    m_gas->getNetProductionRates(&m_wdot[0]);


    
    // Gas phase species
    /*
    for (unsigned i = 0; i != idx0; ++i)
    {
        // For species equations.
        A(i, i) = rho0 * m_u0;
        b(i) = m_wdot[i] * m_W[i];

        // Yk' for other equations, exceptionally here!
        A(idx2, i) = P0 * Wavg * Wavg / m_W[i];
    }
    */
    

    // Continuity equation elements.
    A(0, 0) = rho0;           // u'
    A(0, 1) = m_u0;           // rho'
    A(0, 2) = 0;              // p'

    // Momentum equation elements.
    A(1, 0) = 2* rho0 * m_u0; // u'
    A(1, 1) = m_u0 * m_u0;    // rho'
    A(1, 2) = 1;              // p'

    // State equation elements.
    A(2, 0) = 0;              // u'
    A(2, 1) = RT;             // rho'
    A(2, 2) = -Wavg;          // p'

    // TODO: change b(0) 
    b(0) = 0;                    // RHS continuity
    b(1) = 0;                    // RHS momentum
    b(2) = 0;                    // RHS state

    //Eigen::VectorXd x = A.fullPivLu().solve(b);
    //Eigen::VectorXd::Map(ydot, x.rows()) = x;

    return 0;
}

int PFR1d::evalResidNJ(const double t, const double delta_t, 
                       const double* const y, const double* const ydot, 
                       double* const resid,
                       const Cantera::ResidEval_Type_Enum evalType, 
                       const int id_x, const double delta_x)
{
    const double u = y[0];       // Velocity
    const double r = y[1];       // Density
    const double p = y[2];       // Pressure

    const double dudz = ydot[0];
    const double drdz = ydot[1];
    const double dpdz = ydot[2];

    m_gas->setMassFractions_NoNorm(y+3);
    m_gas->setState_TP(m_T0, y[2]);

    auto loc = 3 + m_nsp;
    for (auto s_ph : m_surf_phases) {
        s_ph->setState_TP(m_T0, y[2]);
        s_ph->setCoveragesNoNorm(y+loc);
        loc += s_ph->nSpecies();
    }

    // Get species production rates
    m_gas->getNetProductionRates(&m_wdot[0]);

    double mdot_surf = evalSurfaces(ydot + m_nsp + 3);

    resid[0] = u * drdz + r * dudz - m_surf_sp_area * mdot_surf;
    resid[1] = 2 * r * u * dudz + u *u * drdz + dpdz;
    resid[2] = m_gas->density() - r;
    for (unsigned k = 0; k < m_gas->nSpecies(); ++k)
    {
        resid[3+k] = u * r * ydot[3+k] + y[3+k] * mdot_surf * m_surf_sp_area 
                     - (m_wdot[k] + m_sdot[k] * m_surf_sp_area) * m_W[k];
    }

    loc = m_nsp; 
    for (auto s_ph : m_surf_phases) {
        double cov_sum = 0.0;
        for (size_t k = 0; k < s_ph->nSpecies()-1; k++) {
            resid[loc + 3 + k] = m_sdot[loc + k];
            cov_sum += y[loc + 3 + k];
        }
        loc += s_ph->nSpecies();
        resid[loc + 2] = 1.0 - cov_sum;

    }
    //resid[2+m_nsp] = 1 - 
    //resid[idx0] = r * dudz + u * drdz;
    //resid[idx1] = u * r * dudz;
    //resid[idx2] = m_gas->density() - r;

    return 0;
}

double PFR1d::evalSurfaceInit() 
{
    const vector_fp& mw = m_gas->molecularWeights();
    vector_fp work(m_nsp);
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    double mdot_surf = 0.0; // net mass flux from surface


    for (auto i = 0; i < m_surf_phases.size(); i++) {
        InterfaceKinetics* kin = m_surf_kins[i];
        SurfPhase* surf = m_surf_phases[i];

        double rs0 = 1.0/surf->siteDensity();
        //surf->setTemperature(m_state[0]);
        kin->getNetProductionRates(work.data());

        for (size_t k = 0; k < m_nsp; k++) {
            m_sdot[k] += work[k];
            mdot_surf += m_sdot[k] * mw[k];
        }
    }

    return mdot_surf;
}

double PFR1d::evalSurfaces(const double* const ydot)
{
    const vector_fp& mw = m_gas->molecularWeights();
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    size_t loc = 0; // offset into ydot
    double mdot_surf = 0.0; // net mass flux from surface

    /*
    for (auto i = 0; i < m_surf_phases.size(); i++) {
        InterfaceKinetics* kin = m_surf_kins[i];
        SurfPhase* surf = m_surf_phases[i];

        double rs0 = 1.0/surf->siteDensity();
        size_t nk = surf->nSpecies();
        double sum = 0.0;
        surf->setTemperature(m_state[0]);
        // S->syncCoverages(); TODO: Add replacement statement
        kin->getNetProductionRates(&m_work[0]);
        size_t ns = kin->surfacePhaseIndex();
        size_t surfloc = kin->kineticsSpeciesIndex(0, ns);
        for (size_t k = 1; k < nk; k++) {
            ydot[loc + k] = m_work[surfloc+k]*rs0*surf->size(k);
            sum -= ydot[loc + k];
        }
        ydot[loc] = sum;
        loc += nk;

        double wallarea = S->area();
        for (size_t k = 0; k < m_nsp; k++) {
            m_sdot[k] += m_work[k]*wallarea;
            mdot_surf += m_sdot[k] * mw[k];
        }
    }
    */
    return mdot_surf;
}

void PFR1d::getSurfaceProductionRates(double* y)
{
    for (size_t i=0; i < m_nsp; i++) {
        y[i] = m_sdot[i];
    }
}

void PFR1d::getSurfaceInitialConditions(double* y)
{
    size_t loc = 0;
    for (const auto s_ph : m_surf_phases) {
        s_ph->getCoverages(y+loc);
        loc += s_ph->nSpecies();
    }
}


}
