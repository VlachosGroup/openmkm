//! @file pfr1d.cpp 
//! Provides an plug-flow reactor model in 1d.

// This file is part of Hetero_ct. 
// Source is adapted from CanteraPFR codebase.
//#include <boost/math/interpolators/barycentric_rational_interpolation.hpp>
#include <boost/range/adaptors.hpp>

#include "pfr1d.h"

using namespace std;
using namespace Cantera;


namespace OpenMKM 
{

PFR1d::PFR1d(IdealGasMix *gas, vector<InterfaceKinetics*> surf_kins,
             vector<SurfPhase*> surf_phases, double area, double cat_abyv,
             double inlet_gas_velocity) 
   : ResidJacEval{}, m_gas(gas), m_surf_kins(surf_kins), 
     m_surf_phases(surf_phases), m_Ac(area), m_cat_abyv(cat_abyv), 
     m_u0(inlet_gas_velocity), m_heat(false), m_energy(0), 
     m_T_interp(nullptr)
{
   
    suppress_thermo_warnings(SUPPRESS_WARNINGS);

    cout << boolalpha
             << "\nIntegrating PFR"
             << "\nUsing Sundials : " << CT_SUNDIALS_VERSION
             << "\nUsing LAPACK   : " << bool(CT_SUNDIALS_USE_LAPACK)
             << endl;

    m_rho_ref = m_gas->density();

    cout << "energy enabled? " << energyEnabled() << endl;
    m_nsp = m_gas->nSpecies();
    if (energyEnabled()) 
        m_neqs_extra = 4;
    else
        m_neqs_extra = 3;
    cout << "m_neqs_extra: " << m_neqs_extra << endl;

    neq_ = m_nsp + m_neqs_extra;
    cout << "neq " << neq_ << endl;
    for (const auto s_ph : m_surf_phases) {
        neq_ += s_ph->nSpecies();
    }
    cout << "neq " << neq_ << endl;

    m_W.resize(m_nsp);
    m_gas->getMolecularWeights(m_W.data());
    m_wdot.resize(m_nsp);
    fill(m_wdot.begin(), m_wdot.end(), 0.0);
    m_sdot.resize(m_nsp);
    fill(m_sdot.begin(), m_sdot.end(), 0.0);

    //m_var.resize(neq_ - m_neqs_extra);
    //m_var = m_gas->speciesNames();
    //m_var.resize(neq_);
    m_var.clear();
    m_var.push_back("Velocity(m/s)");
    m_var.push_back("Density");
    m_var.push_back("Pressure(Pa)");
    if (energyEnabled())
        m_var.push_back("Temperature(K)");
    auto sp_names = m_gas->speciesNames();
    for (auto i = 0; i < sp_names.size(); i++){
        m_var.push_back(sp_names[i]);
    }
    for (const auto s_ph : m_surf_phases) {
        auto sp_nm = s_ph->speciesNames();
        for (auto i = 0; i < sp_nm.size(); i++){
            m_var.push_back(sp_nm[i]);
        }
    }

    m_T0 = gas->temperature();
    m_P0 = gas->pressure();
    cout << "get heat: " << getHeat(m_T0) << endl;
}   

void PFR1d::reinit()
{

    m_rho_ref = m_gas->density();

    cout << "energy enabled? " << energyEnabled() << endl;
    m_nsp = m_gas->nSpecies();
    if (energyEnabled())
        m_neqs_extra = 4;
    else
        m_neqs_extra = 3;
    cout << "m_neqs_extra: " << m_neqs_extra << endl;

    neq_ = m_nsp + m_neqs_extra;
    cout << "neq " << neq_ << endl;
    for (const auto s_ph : m_surf_phases) {
        neq_ += s_ph->nSpecies();
    }
    cout << "neq " << neq_ << endl;

    m_W.resize(m_nsp);
    m_gas->getMolecularWeights(m_W.data());
    m_wdot.resize(m_nsp);
    fill(m_wdot.begin(), m_wdot.end(), 0.0);
    m_sdot.resize(m_nsp);
    fill(m_sdot.begin(), m_sdot.end(), 0.0);

    //m_var.resize(neq_);
    m_var.clear();
    m_var.push_back("Velocity(m/s)");
    m_var.push_back("Density");
    m_var.push_back("Pressure(Pa)");
    if (energyEnabled())
        m_var.push_back("Temperature(K)");
    auto sp_names = m_gas->speciesNames();
    //m_var.insert(m_var.end(), sp_names.begin(), sp_names.end());
    for (auto i = 0; i < sp_names.size(); i++){
        m_var.push_back(sp_names[i]);
    }
    for (const auto s_ph : m_surf_phases) {
        auto sp_nm = s_ph->speciesNames();
        for (auto i = 0; i < sp_nm.size(); i++){
            m_var.push_back(sp_nm[i]);
        }
    }

    m_T0 = m_gas->temperature();
    m_P0 = m_gas->pressure();
}

void PFR1d::setConstraints()
{
    constrain(0, c_GT_ZERO);
    constrain(1, c_GT_ZERO);
    constrain(2, c_GT_ZERO);
    if (energyEnabled()) {
        constrain(3, c_GT_ZERO);
    }
    for (size_t i = 0; i < m_nsp; i++){
        auto k = i + m_neqs_extra;
        constrain(k, c_GE_ZERO);
    }
    auto loc = m_nsp + m_neqs_extra;
    for (const auto s_ph : m_surf_phases) {
        auto nSpecies = s_ph->nSpecies();
        for (size_t i = 0; i < nSpecies; i++){
            auto k = i + loc;
            constrain(k, c_GE_ZERO);
        }
        loc += nSpecies;
    }
}

int PFR1d::getInitialConditions(const double t0,
                                double *const y, 
                                double *const ydot)
{
    //const double P0 = m_gas->pressure();
    const double rho0 = m_gas->density();
    const double Wavg = m_gas->meanMolecularWeight();
    const double RT = m_T0 * GasConstant;
    const double Rrho = rho0 * GasConstant;


    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m_nsp + m_neqs_extra, 
                                              m_nsp + m_neqs_extra);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(m_nsp + m_neqs_extra);

    y[0] = m_u0;
    y[1] = rho0;
    y[2] = m_P0;
    cout << 0 << " y " << y[0] << endl
         << 1 << " y " << y[1] << endl
         << 2 << " y " << y[2] << endl;
    int gas_start_loc;
    if (energyEnabled()) {
        y[3] = m_T0;
        cout << 3 << " y " << y[3] << endl;
    }

    m_gas->getMassFractions(y + m_neqs_extra);
    for (size_t i = 0; i < m_nsp; i++){
        cout << i << " y " << y[i + m_neqs_extra] << endl;
    }

    auto loc = m_nsp + m_neqs_extra;
    for (const auto s_ph : m_surf_phases) {
        s_ph->getCoverages(y + loc);
        loc += s_ph->nSpecies();
    }

    // Continuity equation elements.
    A(0, 0) = rho0;           // u'
    A(0, 1) = m_u0;           // rho'
    A(0, 2) = 0;              // p'
    if (energyEnabled())
        A(0, 3) = 0;          // T'
    for (unsigned k = m_neqs_extra; k < m_nsp + m_neqs_extra; ++k)
        A(0, k) = 0.0;

    // Momentum equation elements.
    A(1, 0) = 2* rho0 * m_u0; // u'
    A(1, 1) = m_u0 * m_u0;    // rho'
    A(1, 2) = 1;              // p'
    if (energyEnabled())
        A(1, 3) = 0;          // T'
    for (unsigned k = m_neqs_extra; k < m_nsp + m_neqs_extra; ++k)
        A(1, k) = 0.0;

    // State equation elements.
    A(2, 0) = 0;              // u'
    A(2, 1) = RT;             // rho'
    A(2, 2) = -Wavg;          // p'
    if (energyEnabled())
        //A(2, 3) = Rrho;       // T'   // This is the source of the difference
        A(2, 3) = 0;       // T'   // This is the source of the difference
    // Yk' for other equations, exceptionally here!
    for (unsigned k = m_neqs_extra; k < m_nsp + m_neqs_extra; ++k)
        A(2, k) = m_P0 * Wavg * Wavg / m_W[k - m_neqs_extra];

    double mdot_surf = evalSurfaces();
    m_gas->getNetProductionRates(&m_wdot[0]);

    if (energyEnabled()) {
        //vector<double> cpr_k(m_nsp);
        //m_gas->getCp_R(cpr_k.data());
        //double cpr = 0.0;
        //for (size_t i = 0; i < m_nsp; i++){
            //cpr += cpr_k[i] * y[i + m_neqs_extra];
        //}
        //cout << "cp" << cpr*GasConstant << endl;
        auto cp = m_gas->cp_mass();

        A(3,0) = 0;            // u'
        A(3,1) = 0;            // rho'
        A(3,2) = 0;            // p'
        A(3,3) = rho0 * m_u0 * cp;// * GasConstant;            // T'
        cout << "A(3,3) " << A(3, 3) << endl;

        vector<double> H_rt(m_nsp);
        m_gas->getEnthalpy_RT(H_rt.data());
        b(3) = 0;                          // RHS energy
        for (size_t i = 0; i < m_nsp; i++){
            //b(3) -= (m_wdot[i] + m_sdot[i] * m_cat_abyv) * m_W[i] * cpr_k[i];
            b(3) -= (m_wdot[i] + m_sdot[i] * m_cat_abyv) *  H_rt[i];
            cout << i << " msdot " << m_sdot[i] << endl
                 //<< i << " m_W " << m_W[i] << endl
                 << i << " H_rt " <<  H_rt[i] << endl;
            cout << "Running b(3) " << b(3) << endl;
            cout << m_cat_abyv << endl;
        }

        cout << "b(3) before multiplication : " << b(3) << endl;
        //b(3) *= m_T0;
        b(3) *= RT;
        cout << "b(3) after multiplication : " << b(3) << endl;
        b(3) += getHeat(m_T0);///GasConstant;
        cout << "get heat: " << getHeat(m_T0) << endl;
        cout << "b(3): " << b(3) << endl;
    }


    b(0) = m_cat_abyv * mdot_surf;          // RHS continuity
    cout << "b(0): " << b(0) << endl;
    b(1) = 0;                               // RHS momentum
    b(2) = 0;                               // RHS state
    //if (energyEnabled()){
    //}

    // Gas phase species
    for (unsigned k = m_neqs_extra; k < m_nsp + m_neqs_extra; ++k)
    {
        cout << " y[k] " << y[k] << endl;
        auto i = k - m_neqs_extra;
        // For species equations.
        A(k, k) = rho0 * m_u0;
        b(k) = (m_wdot[i] + m_sdot[i] * m_cat_abyv) * m_W[i] -
               y[k] * mdot_surf * m_cat_abyv;
        cout << k << " b(k) " << b(k) << endl; 

    }
    
    Eigen::VectorXd x = A.fullPivLu().solve(b);
    Eigen::VectorXd::Map(ydot, x.rows()) = x;
    for (size_t i = 0; i < neq_; i++)
        cout << "i " << i << " ydot[i] " << ydot[i] << endl;

    return 0;
}

int PFR1d::evalResidNJ(const double t, const double delta_t, 
                       const double* const y, const double* const ydot, 
                       double* const resid,
                       const Cantera::ResidEval_Type_Enum evalType, 
                       const int id_x, const double delta_x)
{
    const double u = y[0];       // Flow rate
    const double r = y[1];       // Density
    const double p = y[2];       // Pressure
    const double temp = energyEnabled() ? y[3] : getT(t);   // Temperature
    double RT = GasConstant * temp;

    const double dudz = ydot[0];
    const double drdz = ydot[1];
    const double dpdz = ydot[2];
    double dtempdz = energyEnabled() ? ydot[3] : 0;

    m_gas->setMassFractions_NoNorm(y + m_neqs_extra);
    m_gas->setState_TP(temp, p);

    auto loc = m_neqs_extra + m_nsp;
    for (auto s_ph : m_surf_phases) {
        s_ph->setState_TP(temp, p);
        s_ph->setCoveragesNoNorm(y + loc);
        loc += s_ph->nSpecies();
    }

    // Get species production rates
    m_gas->getNetProductionRates(&m_wdot[0]);
    //double mdot_surf = s();
    vector_fp work(m_nsp);
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    double mdot_surf = 0.0; // net mass flux from surface

    loc = m_nsp;
    for (auto i = 0; i < m_surf_phases.size(); i++) {
        double cov_sum = 0.0;
        InterfaceKinetics* kin = m_surf_kins[i];
        SurfPhase* surf = m_surf_phases[i];
        work.resize(kin->nTotalSpecies());

        kin->getNetProductionRates(work.data());
        for (size_t k = m_nsp + 1; k < kin->nTotalSpecies() - 1; k++){ // Gas + bulk + surface species
            resid[m_neqs_extra + loc + k - m_nsp - 1] = work[k];
            cov_sum += y[loc + m_neqs_extra + k - m_nsp -1];
        }
        loc += surf->nSpecies();
        cov_sum += y[loc + m_neqs_extra - 1];
        resid[loc + m_neqs_extra - 1] = 1.0 - cov_sum; 

        for (size_t k = 0; k < m_nsp; k++) {
            m_sdot[k] += work[k];
            mdot_surf += m_sdot[k] * m_W[k];
        }
    }


    resid[0] = u * drdz + r * dudz - m_cat_abyv * mdot_surf;
    resid[1] = 2 * r * u * dudz + u *u * drdz + dpdz;
    resid[2] = m_gas->density() - r;  // Density is set as pW/RT

    if (energyEnabled()){
        /*
        vector<double> cpr(m_nsp);
        m_gas->getCp_R(cpr.data());
        double cp;
        for (size_t i = 0; i < m_nsp; i++){
            cp += cpr[i] * y[i + m_neqs_extra];
        }
        */
        auto cp = m_gas->cp_mass();
        resid[3] = cp *  r * u * dtempdz;

        double h_term  = 0;                          
        vector<double> H_rt(m_nsp);
        m_gas->getEnthalpy_RT(H_rt.data());
        for (size_t i = 0; i < m_nsp; i++){
            //h_term += (m_wdot[i] + m_sdot[i] * m_cat_abyv) * m_W[i] * cpr[i];
            h_term += (m_wdot[i] + m_sdot[i] * m_cat_abyv)  * H_rt[i];
        }
        //resid[3] += h_term * temp;
        resid[3] += h_term * RT;
        resid[3] -= getHeat(temp);///GasConstant;
    }

    for (unsigned k = 0; k < m_nsp; ++k)
    {
        auto k1 = k + m_neqs_extra;
        resid[k1] = u * r * ydot[k1] + y[k1] * mdot_surf * m_cat_abyv - 
                    (m_wdot[k] + m_sdot[k] * m_cat_abyv) * m_W[k];
    }

    return 0;
}

/*
void PFR1d::setTProfile(map<double, double> T_profile)
{
    if (T_profile.begin()->first > 0)
        m_T_profile_iind = -1;
    else
        m_T_profile_iind = 0;

    m_T_profile.resize(T_profile.size());
    m_T_profile_ind.resize(T_profile.size());
    auto i = 0;
    for (const auto& dist_Tprof : T_profile) {
        m_T_profile_ind[i] = dist_Tprof.first;
        m_T_profile[i++] = dist_Tprof.second;
    }
}

double PFR1d::getT(double z)
{
    if (!m_T_profile.size()){
        return m_T0;
    }

    if (m_T_profile_iind >= 0 && z == m_T_profile_ind[m_T_profile_iind])
        return  m_T_profile[m_T_profile_iind];
    if (m_T_profile_iind < m_T_profile_ind.size() - 1  &&
        z == m_T_profile_ind[m_T_profile_iind + 1])
        return  m_T_profile[m_T_profile_iind + 1];  
    if (m_T_profile_iind == m_T_profile_ind.size() - 1 && 
        z > m_T_profile_ind[m_T_profile_iind])
        return m_T_profile[m_T_profile_iind];

    if (m_T_profile_iind < m_T_profile_ind.size() - 1 && 
        z > m_T_profile_ind[m_T_profile_iind + 1]){
        m_T_profile_iind += 1;
    }
 
    double lz, lT, hz, hT;
    if (m_T_profile_iind < 0) {
        lz = 0;
        lT = m_T0;
    } else {
        lz = m_T_profile_ind[m_T_profile_iind];
        lT = m_T_profile[m_T_profile_iind];
    }
    hz = m_T_profile_ind[m_T_profile_iind + 1];
    hT = m_T_profile[m_T_profile_iind + 1];
    
    return lT + (hT - lT) / (hz - lz) * (z - lz); 
}
*/

/* Even this complexity is not needed 
void PFR1d::setTProfile(const map<double, double>& T_profile)
{
    // The profile given as map in the input is stored two 
    // vectors for boost library interpolators 
    // Clear the existing profile
    if (m_T_profile_ind.size())
        m_T_profile_ind.clear();
    if (m_T_profile.size())
        m_T_profile.clear();

    // If T @ z=0 is not given, use inlet T
    if (T_profile.begin()->first > 0){
        m_T_profile_ind.push_back(0);
        m_T_profile.push_back(m_T0);
    }
    for (const auto& dist_Tprof : T_profile) {
        m_T_profile_ind.push_back(dist_Tprof.first);
        m_T_profile.push_back(dist_Tprof.second);
    }

    m_T_interp = make_shared<boost::math::barycentric_rational<double>>(m_T_profile_ind, m_T_profile, m_T_profile.size());
}
*/

void PFR1d::setTProfile(const map<double, double>& T_profile)
{
    /*
    // The profile given as map in the input is temporarily stored into two 
    // vectors for boost library interpolators 
    // Clear the existing profile
    vector<double> Ts, Zs;

    // If T @ z=0 is not given, use inlet T
    if (T_profile.begin()->first > 0){
        Zs.push_back(0);
        Ts.push_back(m_T0);
    }
    for (const auto& dist_Tprof : T_profile) {
        Zs.push_back(dist_Tprof.first);
        Ts.push_back(dist_Tprof.second);
    }
    */
    auto Zs = boost::adaptors::keys(T_profile);
    auto Ts = boost::adaptors::values(T_profile);

    //if (m_T_interp != nullptr){
    //    m_T_interp = nullptr;
        m_T_interp = make_shared<boost::math::barycentric_rational<double>>(
                Zs.begin(), Zs.end(), Ts.begin());
    //}
}

double PFR1d::getT(double z)
{
    if (m_T_interp == nullptr){
        return m_T0;
    } else{
        return m_T_interp->operator()(z);
    }
}

double PFR1d::evalSurfaces() 
{
    vector_fp work(m_nsp);
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    double mdot_surf = 0.0; // net mass flux from surface


    for (auto i = 0; i < m_surf_phases.size(); i++) {
        InterfaceKinetics* kin = m_surf_kins[i];
        SurfPhase* surf = m_surf_phases[i];
        work.resize(kin->nTotalSpecies());
        fill(work.begin(), work.end(), 0.0);

        kin->getNetProductionRates(work.data());

        for (size_t k = 0; k < m_nsp; k++) {
            m_sdot[k] += work[k];
            mdot_surf += m_sdot[k] * m_W[k];
        }
    }

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

void PFR1d::setHeatTransfer(double htc, double Text, double wall_abyv)
{
    m_heat = true;
    m_htc = htc;
    m_Text = Text;
    m_surf_ext_abyv = wall_abyv;
}

double PFR1d::getHeat(double Tint) const
{
    if (m_heat)
        return m_htc * (m_Text - Tint) * m_surf_ext_abyv;
    else 
        return 0;
   
}

}
