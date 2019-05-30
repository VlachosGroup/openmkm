#include <string>
#include <vector>
#include <memory>
#include <map>
#include <iostream>
#include <fstream>

#include <yaml-cpp/yaml.h>
#include "cantera/base/stringUtils.h"
#include "cantera/base/ct_defs.h"
#include "cantera/IdealGasMix.h"
#include "cantera/InterfaceLatInt.h"
#include "cantera/Interface.h"
#include "cantera/thermo/SurfLatIntPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/LateralInteraction.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/flowControllers.h"

#include "util.h"
#include "reactor.h"
#include "pfr1d.h"
#include "pfr1d_solver.h"


using namespace std;
using namespace Cantera;

namespace HeteroCt 
{

   
void print_rctr_state(double z, PFR1d* rctr, vector<SurfPhase*> surfaces,
                      ofstream& gas_mole_out, ofstream& gas_mass_out,
                      ofstream& gas_msdot_out, ofstream& surf_cov_out,
                      ofstream& state_var_out)
{
 
    vector<double> work(rctr->contents().nSpecies());
 
    
    state_var_out << setw(16) << left  << z
                  << setw(16) << left  << rctr->contents().temperature()
                  << setw(16) << left  << rctr->contents().pressure()
                  << setw(16) << left  << rctr->contents().density()
                  //<< setw(16) << left  << rctr->intEnergy_mass()
                  << endl;
 
    auto gas_var_print = [&work, &z](ostream& out) -> void
    {
        out << setw(16) << left << z;
        for (auto k = 0; k < work.size(); k++) {
            out << setw(16) << left << work[k];
        }
        out << endl;
    };
 
    rctr->contents().getMoleFractions(work.data());
    gas_var_print(gas_mole_out);
 
    rctr->contents().getMassFractions(work.data());
    gas_var_print(gas_mass_out);
 
    rctr->getSurfaceProductionRates(work.data());
    gas_var_print(gas_msdot_out);
    
    
    surf_cov_out << setw(16) << left << z;
    for (auto j = 0;  j <  surfaces.size(); j++) {
        work.resize(surfaces[j]->nSpecies());
        rctr->surface(j)->getCoverages(work.data());
        for (auto k = 0; k < work.size(); k++) {
            surf_cov_out << setw(16) << left << work[k];
        }
    }
    surf_cov_out << endl;
    
}


void run_1d_reactor(YAML::Node& tube_node,
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<Interface>> surfaces,
                    ofstream& gen_info)
{
    //Define the reactor based on the input file
    auto rctr_node = tube_node["reactor"];

    // Read the reactor dimensions
    auto area_node = rctr_node["area"];
    auto len_node = rctr_node["length"];
    //double rctr_vol = strSItoDbl(rctr_node["volume"].as<string>());
    if (!area_node || area_node.IsNull()){
        ;//TODO: Raise Error
    }
    if (!len_node || len_node.IsNull()){
        ;//TODO: Raise Error
    }
    double rctr_xc_area = strSItoDbl(rctr_node["area"].as<string>());
    double rctr_len = strSItoDbl(rctr_node["length"].as<string>());

    auto inlet_node = tube_node["inlet_gas"];
    auto flowrate_node = inlet_node["flow_rate"];     
    auto mfr_node = inlet_node["mass_flow_rate"];
    double velocity{0}, flow_rate{0},  mfr{0};
    if (flowrate_node && !flowrate_node.IsNull()) {
        flow_rate = strSItoDbl(inlet_node["flow_rate"].as<string>());
        velocity = flow_rate / rctr_xc_area;
        mfr = gas->density() * flow_rate;
    }
    else if (mfr_node && !mfr_node.IsNull()) {
        mfr = strSItoDbl(inlet_node["mass_flow_rate"].as<string>());
        flow_rate = mfr / gas->density();
        velocity = flow_rate /rctr_xc_area;
    }

    auto cat_node = rctr_node["cat_abyv"];
    if (surfaces.size() > 0 && (!cat_node || cat_node.IsNull())){
         ; // Throw error
    }
    double cat_abyv = 0.0;
    if (cat_node && !cat_node.IsNull()) {
        cat_abyv = strSItoDbl(rctr_node["cat_abyv"].as<string>());
    }

    vector<InterfaceKinetics*> ikin;
    vector<SurfPhase*> surf_ph;
    for (const auto surf: surfaces) {
        ikin.push_back(surf.get());
        surf_ph.push_back(surf.get());
    }

    // Start the simulation
    cout << "Solving for equilibirum surface coverages at PFR inlet" << endl;
    for (const auto surf: surfaces) {
        surf->solvePseudoSteadyStateProblem();
        vector<double> cov(surf->nSpecies());
        surf->getCoverages(cov.data());

        cout << "Equilibrium surface coverages on Surface: " <<  surf->name() << endl;
        for (auto i = 0; i < surf->nSpecies(); i++)
            gen_info << surf->speciesSPName(i) << " coverage: " << cov[i] << endl;
    }

    cout << "Reached here" << endl;
    cout << flush;
    auto pfr = PFR1d(gas.get(), ikin, surf_ph, rctr_xc_area, cat_abyv, velocity);
    cout << "Reached after PFR initialization" << endl;
    cout << flush;
    string mode = rctr_node["mode"].as<string>();
    cout << "mode " << mode << endl;
    if (mode == "isothermal") {
        pfr.setEnergy(0);
    }
    else if (mode == "tprofile") {
        pfr.setEnergy(0);
        map<double, double> T_profile;
        auto tprofile_nd = rctr_node["TProfile"];
        if (!tprofile_nd || tprofile_nd.IsNull() || !tprofile_nd.IsSequence()){
            ;//TODO: Raise Error
        }
        for(YAML::const_iterator it = tprofile_nd.begin(); it != tprofile_nd.end(); ++it) {
            T_profile.insert(pair<double, double>(it->first.as<double>(),
                                                 it->second.as<double>()));
        }
        pfr.setTProfile(T_profile);
    } 
    else {
        pfr.setEnergy(1);
        if (mode == "heat") {
            double htc = strSItoDbl(rctr_node["htc"].as<string>());
            double wall_abyv = strSItoDbl(rctr_node["wall_abyv"].as<string>());
            double ext_temp = strSItoDbl(rctr_node["Text"].as<string>());
            pfr.setHeatTransfer(htc, ext_temp, wall_abyv);
        }
        pfr.reinit();
    }
    gen_info << "Reactor operating mode: "  << mode << endl;
    gen_info << "Energy enabled? "  << pfr.energyEnabled() << endl;
    
    /*
    vector<double> ydot(25);
    vector<double> y(25);
    pfr.getInitialConditions(0, y.data(), ydot.data());
    for (size_t i = 0; i < 25; i++){
        cout << "i:   " << i << "   y: " << y[i] << "   ydot:   " << ydot[i] << endl;
    }
    */

    PFR1dSolver pfr_solver {&pfr};

    auto simul_node = tube_node["simulation"];
    //auto end_time = strSItoDbl(simul_node["end_time"].as<string>());
    auto solver_node = simul_node["solver"];
    auto abs_tol = solver_node["atol"].as<double>();
    auto rel_tol = solver_node["rtol"].as<double>();
    pfr_solver.setTolerances(rel_tol, abs_tol);

    auto solver_init_step_size_nd = solver_node["init_step_size"];
    if (solver_init_step_size_nd && !solver_init_step_size_nd.IsNull()) {
        double solver_init_step_size = solver_init_step_size_nd.as<double>();
        pfr_solver.setInitialStepSize(solver_init_step_size);
    }

    auto solver_max_steps_nd = solver_node["max_steps"];
    if (solver_max_steps_nd && !solver_max_steps_nd.IsNull()) {
        int solver_max_steps = solver_max_steps_nd.as<int>();
        pfr_solver.setMaxNumSteps(solver_max_steps);
    }
    //pfr_solver.setMaxNumSteps(5000);  //TODO: Make this optional input option
    //
    double simul_init_step_size = 1e-6;
    auto simul_init_step_size_nd = simul_node["init_step"];
    if (simul_init_step_size_nd && !simul_init_step_size_nd.IsNull()) {
        simul_init_step_size = simul_init_step_size_nd.as<double>();
    }
    //vector<double> zvals = get_log10_intervals(rctr_len, 1e-7); //Use the same function to get z steps
    vector<double> zvals = get_log10_intervals(rctr_len, simul_init_step_size); //Use the same function to get z steps

    ofstream gas_mole_out("gas_mole_ss.out");
    ofstream gas_mass_out("gas_mass_ss.out");
    ofstream gas_msdot_out("gas_msdot_ss.out");
    ofstream surf_cov_out("surf_cov_ss.out");
    ofstream state_var_out("rctr_state_ss.out");

    gas_mole_out.precision(6);
    gas_mass_out.precision(6);
    gas_msdot_out.precision(6);
    surf_cov_out.precision(6);
    state_var_out.precision(6);

    for (const auto& z : zvals) {
        pfr_solver.solve(z);
        print_rctr_state(z, &pfr, surf_ph, gas_mole_out, gas_mass_out, 
                         gas_msdot_out, surf_cov_out, state_var_out);
    }
    
    pfr_solver.writeStateData("1d_pfr_state.out");
    pfr_solver.writeGasData("1d_pfr_gas.out");
    pfr_solver.writeSurfaceData("1d_pfr_surface.out");

    //vector<double> gas_X(gas->nSpecies());

}

void run_1d_reactor(YAML::Node& tube_node,
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<InterfaceInteractions>> surfaces,
                    ofstream& gen_info)
{
    //Define the reactor based on the input file
    auto rctr_node = tube_node["reactor"];

    // Read the reactor dimensions
    auto area_node = rctr_node["area"];
    auto len_node = rctr_node["length"];
    //double rctr_vol = strSItoDbl(rctr_node["volume"].as<string>());
    if (!area_node || area_node.IsNull()){
        ;//TODO: Raise Error
    }
    if (!len_node || len_node.IsNull()){
        ;//TODO: Raise Error
    }
    double rctr_xc_area = strSItoDbl(rctr_node["area"].as<string>());
    double rctr_len = strSItoDbl(rctr_node["length"].as<string>());

    auto inlet_node = tube_node["inlet_gas"];
    auto flowrate_node = inlet_node["flow_rate"];     
    auto mfr_node = inlet_node["mass_flow_rate"];
    double velocity{0}, flow_rate{0},  mfr{0};
    if (flowrate_node && !flowrate_node.IsNull()) {
        flow_rate = strSItoDbl(inlet_node["flow_rate"].as<string>());
        velocity = flow_rate / rctr_xc_area;
        mfr = gas->density() * flow_rate;
    }
    else if (mfr_node && !mfr_node.IsNull()) {
        mfr = strSItoDbl(inlet_node["mass_flow_rate"].as<string>());
        flow_rate = mfr / gas->density();
        velocity = flow_rate /rctr_xc_area;
    }

    auto cat_node = rctr_node["cat_abyv"];
    if (surfaces.size() > 0 && (!cat_node || cat_node.IsNull())){
         ; // Throw error
    }
    double cat_abyv = 0.0;
    if (cat_node && !cat_node.IsNull()) {
        cat_abyv = strSItoDbl(rctr_node["cat_abyv"].as<string>());
    }

    vector<InterfaceKinetics*> ikin;
    vector<SurfPhase*> surf_ph;
    for (const auto surf: surfaces) {
        ikin.push_back(surf.get());
        surf_ph.push_back(surf.get());
    }

    // Start the simulation
    cout << "Solving for equilibirum surface coverages at PFR inlet" << endl;
    for (const auto surf: surfaces) {
        surf->solvePseudoSteadyStateProblem();
        vector<double> cov(surf->nSpecies());
        surf->getCoverages(cov.data());

        cout << "Equilibrium surface coverages on Surface: " <<  surf->name() << endl;
        for (auto i = 0; i < surf->nSpecies(); i++)
            gen_info << surf->speciesSPName(i) << " coverage: " << cov[i] << endl;
    }

    cout << "Reached here" << endl;
    cout << flush;
    auto pfr = PFR1d(gas.get(), ikin, surf_ph, rctr_xc_area, cat_abyv, velocity);
    cout << "Reached after PFR initialization" << endl;
    cout << flush;
    string mode = rctr_node["mode"].as<string>();
    cout << "mode " << mode << endl;
    if (mode == "isothermal") {
        pfr.setEnergy(0);
    }
    else if (mode == "tprofile") {
        pfr.setEnergy(0);
        map<double, double> T_profile;
        auto tprofile_nd = rctr_node["TProfile"];
        if (!tprofile_nd || tprofile_nd.IsNull() || !tprofile_nd.IsSequence()){
            ;//TODO: Raise Error
        }
        for(YAML::const_iterator it = tprofile_nd.begin(); it != tprofile_nd.end(); ++it) {
            T_profile.insert(pair<double, double>(it->first.as<double>(),
                                                 it->second.as<double>()));
        }
        pfr.setTProfile(T_profile);
    } 
    else {
        pfr.setEnergy(1);
        if (mode == "heat") {
            double htc = strSItoDbl(rctr_node["htc"].as<string>());
            double wall_abyv = strSItoDbl(rctr_node["wall_abyv"].as<string>());
            double ext_temp = strSItoDbl(rctr_node["Text"].as<string>());
            pfr.setHeatTransfer(htc, ext_temp, wall_abyv);
        }
        pfr.reinit();
    }
    gen_info << "Reactor operating mode: "  << mode << endl;
    gen_info << "Energy enabled? "  << pfr.energyEnabled() << endl;
    
    /*
    vector<double> ydot(25);
    vector<double> y(25);
    pfr.getInitialConditions(0, y.data(), ydot.data());
    for (size_t i = 0; i < 25; i++){
        cout << "i:   " << i << "   y: " << y[i] << "   ydot:   " << ydot[i] << endl;
    }
    */

    PFR1dSolver pfr_solver {&pfr};

    auto simul_node = tube_node["simulation"];
    //auto end_time = strSItoDbl(simul_node["end_time"].as<string>());
    auto solver_node = simul_node["solver"];
    auto abs_tol = solver_node["atol"].as<double>();
    auto rel_tol = solver_node["rtol"].as<double>();
    pfr_solver.setTolerances(rel_tol, abs_tol);

    auto solver_init_step_size_nd = solver_node["init_step_size"];
    if (solver_init_step_size_nd && !solver_init_step_size_nd.IsNull()) {
        double solver_init_step_size = solver_init_step_size_nd.as<double>();
        pfr_solver.setInitialStepSize(solver_init_step_size);
    }

    auto solver_max_steps_nd = solver_node["max_steps"];
    if (solver_max_steps_nd && !solver_max_steps_nd.IsNull()) {
        int solver_max_steps = solver_max_steps_nd.as<int>();
        pfr_solver.setMaxNumSteps(solver_max_steps);
    }
    //pfr_solver.setMaxNumSteps(5000);  //TODO: Make this optional input option
    //
    double simul_init_step_size = 1e-6;
    auto simul_init_step_size_nd = simul_node["init_step"];
    if (simul_init_step_size_nd && !simul_init_step_size_nd.IsNull()) {
        simul_init_step_size = simul_init_step_size_nd.as<double>();
    }
    //vector<double> zvals = get_log10_intervals(rctr_len, 1e-7); //Use the same function to get z steps
    vector<double> zvals = get_log10_intervals(rctr_len, simul_init_step_size); //Use the same function to get z steps

    ofstream gas_mole_out("gas_mole_ss.out");
    ofstream gas_mass_out("gas_mass_ss.out");
    ofstream gas_msdot_out("gas_msdot_ss.out");
    ofstream surf_cov_out("surf_cov_ss.out");
    ofstream state_var_out("rctr_state_ss.out");

    gas_mole_out.precision(6);
    gas_mass_out.precision(6);
    gas_msdot_out.precision(6);
    surf_cov_out.precision(6);
    state_var_out.precision(6);

    for (const auto& z : zvals) {
        pfr_solver.solve(z);
        print_rctr_state(z, &pfr, surf_ph, gas_mole_out, gas_mass_out, 
                         gas_msdot_out, surf_cov_out, state_var_out);
    }
    
    pfr_solver.writeStateData("1d_pfr_state.out");
    pfr_solver.writeGasData("1d_pfr_gas.out");
    pfr_solver.writeSurfaceData("1d_pfr_surface.out");

    //vector<double> gas_X(gas->nSpecies());

}

}
