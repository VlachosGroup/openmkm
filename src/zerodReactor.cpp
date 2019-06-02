#include <string>
#include <vector>
#include <memory>
#include <map>
#include <iostream>
#include <iomanip>
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
#include "cantera/zeroD/Wall.h"
#include "IdealGasTRampReactor.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/flowControllers.h"

#include "util.h"
#include "reactor.h"
#include "hctexceptions.h"


namespace HeteroCt 
{

using namespace std;
using namespace Cantera;

//TODO: Add nSurfaces() function to Reactor/ReactorBase to eliminate the need
// to pass surfaces argument
void print_rctr_state(double z, Reactor* rctr, vector<SurfPhase*> surfaces, 
                      ofstream& gas_mole_out, ofstream& gas_mass_out, 
                      ofstream& gas_msdot_out, ofstream& surf_cov_out,
                      ofstream& state_var_out)
{

    vector<double> work(rctr->contents().nSpecies());

    state_var_out << setw(16) << left  << z
                  << setw(16) << left  << rctr->temperature()
                  << setw(16) << left  << rctr->pressure() 
                  << setw(16) << left  << rctr->density() 
                  << setw(16) << left  << rctr->intEnergy_mass() 
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


void run_0d_reactor(RctrType rctr_type, 
                    YAML::Node& tube_node,
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<InterfaceInteractions>> surfaces,
                    ofstream& gen_info_out) 
{
    //Define the reactor based on the input file
    auto rctr_node = tube_node["reactor"];

    auto rctr = make_shared<IdealGasTRampReactor>();
    rctr->insert(*gas);

    // Read the reactor dimensions

    double rctr_vol = strSItoDbl(rctr_node["volume"].as<string>());
    if (!rctr_vol) {
        // TODO: Raise Error
        ;
    }
    size_t rctr_nos = 1;
    if (rctr_type == PFR_0D) { // Modify rctr_nos only for PFR_0D
        auto nd_node = rctr_node["nodes"];
        if (nd_node)
            if (!nd_node.IsNull())
                rctr_nos = nd_node.as<int>();

        if (rctr_nos == 1) { // Raise warning to increase cstr number
            cout << "Number of nodes in 0d PFR simulation is 1. \n " 
                 << "Suggestion: Input 'nodes' parameter if not given already or " 
                 << "increase its value to greater than 1.";
        }
        rctr_vol /= rctr_nos;
    }

    rctr->setInitialVolume(rctr_vol);
    vector<shared_ptr<ReactorSurface>> cat_surfs;
    vector<SurfPhase*> surf_phases;
    auto cat_node = rctr_node["cat_abyv"];
    if (surfaces.size() > 0 && (!cat_node || cat_node.IsNull())){
        throw YAMLParserError("zerodReactor.cpp::runZerodReactor", "cat_abyv", 
                              "Not found or null");
    }

    if (cat_node && !cat_node.IsNull()) {
        double cat_area = strSItoDbl(cat_node.as<string>());
        cat_area *= rctr_vol;
        size_t surf_spec_no = 0;
        for (const auto surf : surfaces) {
            surf_spec_no += surf->nSpecies();
            auto cat_surf = make_shared<ReactorSurface>(); 
            cat_surf->setKinetics(surf.get());
            surf_phases.push_back(dynamic_cast<SurfPhase*> (surf.get()));
            vector<double> coverages(surf->nSpecies());
            surf->getCoverages(coverages.data());
            cat_surf->setCoverages(coverages.data());
            cat_surf->setArea(cat_area);
            cat_surfs.push_back(cat_surf);
            rctr->addSurface(cat_surf.get());
        }
    }

    

    // Only CSTR and PFR_0D require valves, controllers and 
    // reservoirs to facilitate mass transfer. Attach them to reactor

    // The following are nominally defined in the case of BATCH reactor also 
    // to prevent crashing the program
    auto in_rsrv = make_shared<Reservoir>();
    auto exhst = make_shared<Reservoir>();
    auto inlet_mfc = make_shared<MassFlowController>();
    auto outlet = make_shared<PressureController>();
    
    if (rctr_type != BATCH) {
        auto inlet_node = tube_node["inlet_gas"];
        auto flowrate_node = inlet_node["flow_rate"];
        auto restime_node = inlet_node["residence_time"];
        auto mfr_node = inlet_node["mass_flow_rate"];
        double flowrate{0}, residence_time{0}, mfr{0};
        if (flowrate_node && !flowrate_node.IsNull()) {
            flowrate = strSItoDbl(inlet_node["flow_rate"].as<string>());
            //residence_time = rctr_vol/flowrate;
            mfr = rctr->mass() * flowrate / rctr_vol;
        }
        else if (restime_node && !restime_node.IsNull()) {
            residence_time = strSItoDbl(inlet_node["residence_time"].as<string>());
            mfr = rctr->mass()/residence_time;
        }
        else {
            mfr = strSItoDbl(inlet_node["mass_flow_rate"].as<string>());
        }
        inlet_mfc->setMassFlowRate(mfr);
        outlet->setMaster(inlet_mfc.get());
        outlet->setPressureCoeff(0.00001);
        in_rsrv->insert(*gas);
        exhst->insert(*gas);
        inlet_mfc->install(*in_rsrv, *rctr);
        outlet->install(*rctr, *exhst);
    }

    gen_info_out << "reactor density " << rctr->density() << endl;

    rctr->setChemistry();

    // Read the reactor mode and set corresponding parameters
    string mode = rctr_node["mode"].as<string>();
    // heat_rsrv and wall are nominally defined.
    // They are used if heat transfer is required.
    auto heat_rsrv = make_shared<Reservoir>(); 
    auto wall = make_shared<Wall>();
    if (mode == "isothermal" || mode == "tpd")  
        rctr->setEnergy(0);
    else {
        rctr->setEnergy(1);
        if (mode == "heat"){
            double htc = strSItoDbl(rctr_node["htc"].as<string>());
            double wall_abyv = strSItoDbl(rctr_node["wall_abyv"].as<string>());
            double ext_temp = strSItoDbl(rctr_node["Text"].as<string>());
            auto press = gas->pressure();
            gas->setState_TP(ext_temp, press);
            heat_rsrv->insert(*gas);
            wall->setHeatTransferCoeff(htc);
            wall->setArea(wall_abyv * rctr->volume());
            wall->install(*heat_rsrv, *rctr);
        }
    }

    // Read simulation parameters
    auto simul_node = tube_node["simulation"];
    bool transient_log = simul_node["transient"].as<bool>();
    double end_time = 0;
    string stepping;
    double step_size = 0;
    if (mode == "tpd"){
        auto beg_temp = strSItoDbl(rctr_node["temperature"].as<string>());
        auto end_temp = strSItoDbl(rctr_node["Tend"].as<string>());
        auto temp_ramp = strSItoDbl(rctr_node["Tramp"].as<string>());
        rctr->setBeta(temp_ramp);
        end_time = (end_temp - beg_temp)/temp_ramp;
        cout << "Tramp " << temp_ramp << endl;
        cout << "end time  " << end_time << endl;

        stepping = "regular";
        transient_log = true;
    } else {
        end_time = strSItoDbl(simul_node["end_time"].as<string>());
    }

    vector<double> times;
    if (transient_log && rctr_nos == 1){ // Only for singler CSTR or a batch reactor
        if (! (mode == "tpd")) { // for tpd, stepping set already
            stepping = simul_node["stepping"].as<string>(); // regular or logarithmic
        }
        auto step_size = strSItoDbl(simul_node["step_size"].as<string>());
        cout << "Step size " << step_size << endl;

        if (stepping == "logarithmic"){
            times = get_log10_intervals(end_time, step_size);
        }
        else if (stepping == "regular"){
            times = get_reg_intervals(0, end_time, step_size);
        } 
    } else{
        times.push_back(end_time);
    }
    cout << "Size of time vector " << times.size() << endl;
    
    // Setup simulation 
    ReactorNet rnet; 
    rnet.addReactor(*rctr);

    // Pass any user defined numerical options to the solver
    auto solver_node = simul_node["solver"];
    if (solver_node && !solver_node.IsNull()){
        auto abs_tol_node = solver_node["atol"];
        if (abs_tol_node && !abs_tol_node.IsNull()){
            auto abs_tol = abs_tol_node.as<double>();
            auto rel_tol_node = solver_node["rtol"]; 
            if (rel_tol_node && !rel_tol_node.IsNull()){
                auto rel_tol = rel_tol_node.as<double>();
                rnet.setTolerances(rel_tol, abs_tol);
            }
            else {
                cout << "WARNING: Both atol and rtol are required" << endl;
            }
        }
        // Implement the other options here 
    }

    auto gas_print_specie_header = [&gas](string ind_var, ostream& out) -> void
    {
        out << setw(16) << left << ind_var;
        for (const auto & sp_name : gas->speciesNames()) {
            out << setw(16) << left << sp_name;
        }
        out << endl;
    };


    // Steady state condition makes sense For CSTR and PFR_0D 
    ofstream gas_ss_mole_out ("gas_mole_ss.out", ios::out);
    ofstream gas_ss_mass_out ("gas_mass_ss.out", ios::out);
    ofstream gas_ss_msdot_out ("gas_msdot_ss.out", ios::out);
    ofstream surf_ss_out ("surf_cov_ss.out", ios::out);
    ofstream state_var_out ("rctr_state.out", ios::out);

    if (rctr_type != BATCH) {
        gas_ss_mole_out 
            << "Gas Mole fractions at Steady State "  << endl;
        gas_print_specie_header("z(m)", gas_ss_mole_out);

        gas_ss_mass_out 
            << "Gas Mass fractions at Steady State "  << endl;
        gas_print_specie_header("z(m)", gas_ss_mass_out);

        gas_ss_msdot_out 
            << "Surface Production Rates of  Gas Species at Steady State "  
            << endl;
        gas_print_specie_header("z(m)", gas_ss_msdot_out);

        surf_ss_out 
            << "Surace Coverages at Steady State: "  << endl 
            << setw(16) << left << "z(m)";
        for (const auto surf : surfaces) {
            for (const auto & sp_name : surf->speciesNames()) {
                surf_ss_out << setw(16) << left << sp_name;
            }
        }
        surf_ss_out << endl;

        state_var_out 
            << "Reactor State Variables at Steady State" << endl
            << setw(16) << left << "z(m)"
            << setw(16) << left << "Temperature(K)" 
            << setw(16) << left << "Pressure(Pa)" 
            << setw(16) << left << "Density(kg/m3)" 
            << setw(16) << left << "U(J/kg)" 
            << endl;
        
        // Print the inlet state
        rnet.reinitialize();
        print_rctr_state(0, rctr.get(), surf_phases, gas_ss_mole_out, 
                         gas_ss_mass_out, gas_ss_msdot_out, surf_ss_out, 
                         state_var_out);
    }

    // Transient state makes sense For BATCH and CSTR 
    ofstream gas_tr_mole_out ("gas_mole_tr.out", ios::out);
    ofstream gas_tr_mass_out ("gas_mass_tr.out", ios::out);
    ofstream gas_tr_msdot_out ("gas_msdot_tr.out", ios::out);
    ofstream surf_tr_out ("surf_cov_tr.out", ios::out);
    ofstream state_var_tr_out ("rctr_state_tr.out", ios::out);
    if (rctr_type != PFR_0D && times.size() > 1) {
        gas_tr_mole_out 
            << "Transient Gas Mole fractions"  << endl;
        gas_print_specie_header("t(s)", gas_tr_mole_out);

        gas_tr_mass_out 
            << "Transient Gas Mass fractions"  << endl;
        gas_print_specie_header("t(s)", gas_tr_mass_out);

        gas_tr_msdot_out 
            << "Transient Surface Production Rates of  Gas Species"  << endl;
        gas_print_specie_header("t(s)", gas_tr_msdot_out);

        surf_tr_out 
            << "Transient Surace Coverages"  << endl 
            << setw(16) << left << "t(s)";
        for (const auto surf : surfaces) {
            for (const auto & sp_name : surf->speciesNames()) {
                surf_tr_out << setw(16) << left << sp_name;
            }
        }
        surf_tr_out << endl;

        state_var_tr_out 
            << "Transient Reactor State"  << endl
            << setw(16) << left << "t(s)" 
            << setw(16) << left << "Temperature(K)" 
            << setw(16) << left << "Pressure(Pa)" 
            << setw(16) << left << "Density(kg/m3)" 
            << setw(16) << left << "U(J/kg)" 
            << endl;
     
    }

    vector<double> gas_X(gas->nSpecies());

    for (size_t i = 0; i < rctr_nos; i++) {
        if (rctr_nos > 1) {
            gen_info_out << "CSTR #: " << i << endl;
        }
        // The next 9 lines redefine the TPX of gas for each reactor
        double temp = rctr->contents().temperature();
        double pressure = rctr->contents().pressure();
        rctr->contents().getMoleFractions(gas_X.data());
        gas->setState_TPX(temp, pressure, gas_X.data());
        if (rctr_type != BATCH){
            in_rsrv->syncState();
        }
        rnet.setInitialTime(0);
        rnet.reinitialize();
        rnet.setMaxTimeStep(1e-1);

        double T_init = temp;
        double beta = 0;
        if (mode == "tpd") {
            beta = strSItoDbl(rctr_node["Tramp"].as<string>());
            cout << "beta " << beta << endl;
        }

        for (const auto & tm : times) {
            rnet.advance(tm);
            if (transient_log) {
                print_rctr_state(tm, rctr.get(), surf_phases, 
                                 gas_tr_mole_out, gas_tr_mass_out, 
                                 gas_tr_msdot_out, surf_tr_out, 
                                 state_var_tr_out);

            }
        }
        rctr->restoreState();

        print_rctr_state((i+0.5)*rctr_vol, rctr.get(), surf_phases, 
                            gas_ss_mole_out, gas_ss_mass_out, 
                            gas_ss_msdot_out, surf_ss_out, state_var_out);
    }
}

void run_0d_reactor(RctrType rctr_type, 
                    YAML::Node& tube_node,
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<Interface>> surfaces,
                    ofstream& gen_info_out) 
{
    //Define the reactor based on the input file
    auto rctr_node = tube_node["reactor"];

    auto rctr = make_shared<IdealGasTRampReactor>();
    rctr->insert(*gas);

    // Read the reactor dimensions

    double rctr_vol = strSItoDbl(rctr_node["volume"].as<string>());
    if (!rctr_vol) {
        // TODO: Raise Error
        ;
    }
    size_t rctr_nos = 1;
    if (rctr_type == PFR_0D) { // Modify rctr_nos only for PFR_0D
        auto nd_node = rctr_node["nodes"];
        if (nd_node)
            if (!nd_node.IsNull())
                rctr_nos = nd_node.as<int>();

        if (rctr_nos == 1) { // Raise warning to increase cstr number
            cout << "Number of nodes in 0d PFR simulation is 1. \n " 
                 << "Suggestion: Input 'nodes' parameter if not given already or " 
                 << "increase its value to greater than 1.";
        }
        rctr_vol /= rctr_nos;
    }

    rctr->setInitialVolume(rctr_vol);
    vector<shared_ptr<ReactorSurface>> cat_surfs;
    vector<SurfPhase*> surf_phases;
    auto cat_node = rctr_node["cat_abyv"];
    if (surfaces.size() > 0 && (!cat_node || cat_node.IsNull())){
        throw YAMLParserError("zerodReactor.cpp::runZerodReactor", "cat_abyv", 
                              "Not found or null");
    }

    if (cat_node && !cat_node.IsNull()) {
        double cat_area = strSItoDbl(cat_node.as<string>());
        cat_area *= rctr_vol;
        size_t surf_spec_no = 0;
        for (const auto surf : surfaces) {
            surf_spec_no += surf->nSpecies();
            auto cat_surf = make_shared<ReactorSurface>(); 
            cat_surf->setKinetics(surf.get());
            surf_phases.push_back(dynamic_cast<SurfPhase*> (surf.get()));
            vector<double> coverages(surf->nSpecies());
            surf->getCoverages(coverages.data());
            cat_surf->setCoverages(coverages.data());
            cat_surf->setArea(cat_area);
            cat_surfs.push_back(cat_surf);
            rctr->addSurface(cat_surf.get());
        }
    }

    rctr->setChemistry();
    string mode = rctr_node["mode"].as<string>();
    if (mode == "isothermal" || mode == "tpd")  
        rctr->setEnergy(0);
    else 
        rctr->setEnergy(1);

    // Only CSTR and PFR_0D require valves, controllers and 
    // reservoirs to facilitate mass transfer. Attach them to reactor

    // The following are nominally defined in the case of BATCH reactor also 
    // to prevent crashing the program
    auto in_rsrv = make_shared<Reservoir>();
    auto exhst = make_shared<Reservoir>();
    auto inlet_mfc = make_shared<MassFlowController>();
    auto outlet = make_shared<PressureController>();
    
    if (rctr_type != BATCH) {
        auto inlet_node = tube_node["inlet_gas"];
        auto flowrate_node = inlet_node["flow_rate"];
        auto restime_node = inlet_node["residence_time"];
        auto mfr_node = inlet_node["mass_flow_rate"];
        double flowrate{0}, residence_time{0}, mfr{0};
        if (flowrate_node && !flowrate_node.IsNull()) {
            flowrate = strSItoDbl(inlet_node["flow_rate"].as<string>());
            //residence_time = rctr_vol/flowrate;
            mfr = rctr->mass() * flowrate / rctr_vol;
        }
        else if (restime_node && !restime_node.IsNull()) {
            residence_time = strSItoDbl(inlet_node["residence_time"].as<string>());
            mfr = rctr->mass()/residence_time;
        }
        else {
            mfr = strSItoDbl(inlet_node["mass_flow_rate"].as<string>());
        }
        inlet_mfc->setMassFlowRate(mfr);
        outlet->setMaster(inlet_mfc.get());
        outlet->setPressureCoeff(0.00001);
        in_rsrv->insert(*gas);
        exhst->insert(*gas);
        inlet_mfc->install(*in_rsrv, *rctr);
        outlet->install(*rctr, *exhst);
    }

    gen_info_out << "reactor density " << rctr->density() << endl;

    // Read simulation parameters
    auto simul_node = tube_node["simulation"];
    bool transient_log = simul_node["transient"].as<bool>();
    double end_time = 0;
    string stepping;
    double step_size = 0;
    if (mode == "tpd"){
        auto beg_temp = strSItoDbl(rctr_node["temperature"].as<string>());
        auto end_temp = strSItoDbl(rctr_node["Tend"].as<string>());
        auto temp_ramp = strSItoDbl(rctr_node["Tramp"].as<string>());
        rctr->setBeta(temp_ramp);
        end_time = (end_temp - beg_temp)/temp_ramp;
        cout << "Tramp " << temp_ramp << endl;
        cout << "end time  " << end_time << endl;

        stepping = "regular";
        transient_log = true;
    } else {
        end_time = strSItoDbl(simul_node["end_time"].as<string>());
    }

    vector<double> times;
    if (transient_log && rctr_nos == 1){ // Only for singler CSTR or a batch reactor
        if (! (mode == "tpd")) { // for tpd, stepping set already
            stepping = simul_node["stepping"].as<string>(); // regular or logarithmic
        }
        auto step_size = strSItoDbl(simul_node["step_size"].as<string>());
        cout << "Step size " << step_size << endl;

        if (stepping == "logarithmic"){
            times = get_log10_intervals(end_time, step_size);
        }
        else if (stepping == "regular"){
            times = get_reg_intervals(0, end_time, step_size);
        } 
    } else{
        times.push_back(end_time);
    }
    cout << "Size of time vector " << times.size() << endl;
    
    // Setup simulation 
    ReactorNet rnet; 
    rnet.addReactor(*rctr);

    // Pass any user defined numerical options to the solver
    auto solver_node = simul_node["solver"];
    if (solver_node && !solver_node.IsNull()){
        auto abs_tol_node = solver_node["atol"];
        if (abs_tol_node && !abs_tol_node.IsNull()){
            auto abs_tol = abs_tol_node.as<double>();
            auto rel_tol_node = solver_node["rtol"]; 
            if (rel_tol_node && !rel_tol_node.IsNull()){
                auto rel_tol = rel_tol_node.as<double>();
                rnet.setTolerances(rel_tol, abs_tol);
            }
            else {
                cout << "WARNING: Both atol and rtol are required" << endl;
            }
        }
        // Implement the other options here 
    }

    auto gas_print_specie_header = [&gas](string ind_var, ostream& out) -> void
    {
        out << setw(16) << left << ind_var;
        for (const auto & sp_name : gas->speciesNames()) {
            out << setw(16) << left << sp_name;
        }
        out << endl;
    };


    // Steady state condition makes sense For CSTR and PFR_0D 
    ofstream gas_ss_mole_out ("gas_mole_ss.out", ios::out);
    ofstream gas_ss_mass_out ("gas_mass_ss.out", ios::out);
    ofstream gas_ss_msdot_out ("gas_msdot_ss.out", ios::out);
    ofstream surf_ss_out ("surf_cov_ss.out", ios::out);
    ofstream state_var_out ("rctr_state.out", ios::out);

    if (rctr_type != BATCH) {
        gas_ss_mole_out 
            << "Gas Mole fractions at Steady State "  << endl;
        gas_print_specie_header("z(m)", gas_ss_mole_out);

        gas_ss_mass_out 
            << "Gas Mass fractions at Steady State "  << endl;
        gas_print_specie_header("z(m)", gas_ss_mass_out);

        gas_ss_msdot_out 
            << "Surface Production Rates of  Gas Species at Steady State "  
            << endl;
        gas_print_specie_header("z(m)", gas_ss_msdot_out);

        surf_ss_out 
            << "Surace Coverages at Steady State: "  << endl 
            << setw(16) << left << "z(m)";
        for (const auto surf : surfaces) {
            for (const auto & sp_name : surf->speciesNames()) {
                surf_ss_out << setw(16) << left << sp_name;
            }
        }
        surf_ss_out << endl;

        state_var_out 
            << "Reactor State Variables at Steady State" << endl
            << setw(16) << left << "z(m)"
            << setw(16) << left << "Temperature(K)" 
            << setw(16) << left << "Pressure(Pa)" 
            << setw(16) << left << "Density(kg/m3)" 
            << setw(16) << left << "U(J/kg)" 
            << endl;
        
        // Print the inlet state
        rnet.reinitialize();
        print_rctr_state(0, rctr.get(), surf_phases, gas_ss_mole_out, 
                         gas_ss_mass_out, gas_ss_msdot_out, surf_ss_out, 
                         state_var_out);
    }

    // Transient state makes sense For BATCH and CSTR 
    ofstream gas_tr_mole_out ("gas_mole_tr.out", ios::out);
    ofstream gas_tr_mass_out ("gas_mass_tr.out", ios::out);
    ofstream gas_tr_msdot_out ("gas_msdot_tr.out", ios::out);
    ofstream surf_tr_out ("surf_cov_tr.out", ios::out);
    ofstream state_var_tr_out ("rctr_state_tr.out", ios::out);
    if (rctr_type != PFR_0D && times.size() > 1) {
        gas_tr_mole_out 
            << "Transient Gas Mole fractions"  << endl;
        gas_print_specie_header("t(s)", gas_tr_mole_out);

        gas_tr_mass_out 
            << "Transient Gas Mass fractions"  << endl;
        gas_print_specie_header("t(s)", gas_tr_mass_out);

        gas_tr_msdot_out 
            << "Transient Surface Production Rates of  Gas Species"  << endl;
        gas_print_specie_header("t(s)", gas_tr_msdot_out);

        surf_tr_out 
            << "Transient Surace Coverages"  << endl 
            << setw(16) << left << "t(s)";
        for (const auto surf : surfaces) {
            for (const auto & sp_name : surf->speciesNames()) {
                surf_tr_out << setw(16) << left << sp_name;
            }
        }
        surf_tr_out << endl;

        state_var_tr_out 
            << "Transient Reactor State"  << endl
            << setw(16) << left << "t(s)" 
            << setw(16) << left << "Temperature(K)" 
            << setw(16) << left << "Pressure(Pa)" 
            << setw(16) << left << "Density(kg/m3)" 
            << setw(16) << left << "U(J/kg)" 
            << endl;
     
    }

    vector<double> gas_X(gas->nSpecies());

    for (size_t i = 0; i < rctr_nos; i++) {
        if (rctr_nos > 1) {
            gen_info_out << "CSTR #: " << i << endl;
        }
        // The next 9 lines redefine the TPX of gas for each reactor
        double temp = rctr->contents().temperature();
        double pressure = rctr->contents().pressure();
        rctr->contents().getMoleFractions(gas_X.data());
        gas->setState_TPX(temp, pressure, gas_X.data());
        if (rctr_type != BATCH){
            in_rsrv->syncState();
        }
        rnet.setInitialTime(0);
        rnet.reinitialize();
        rnet.setMaxTimeStep(1e-1);

        double T_init = temp;
        double beta = 0;
        if (mode == "tpd") {
            beta = strSItoDbl(rctr_node["Tramp"].as<string>());
            cout << "beta " << beta << endl;
        }

        for (const auto & tm : times) {
            rnet.advance(tm);
            if (transient_log) {
                print_rctr_state(tm, rctr.get(), surf_phases, 
                                 gas_tr_mole_out, gas_tr_mass_out, 
                                 gas_tr_msdot_out, surf_tr_out, 
                                 state_var_tr_out);

            }
        }
        rctr->restoreState();

        print_rctr_state((i+0.5)*rctr_vol, rctr.get(), surf_phases, 
                            gas_ss_mole_out, gas_ss_mass_out, 
                            gas_ss_msdot_out, surf_ss_out, state_var_out);
    }
}

}
