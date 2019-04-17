#include <string>
#include <memory>
#include <map>
#include <iostream>
#include <fstream>

#include <yaml-cpp/yaml.h>
#include "cantera/base/stringUtils.h"
#include "cantera/base/ct_defs.h"
#include "cantera/IdealGasMix.h"
#include "cantera/InterfaceLatInt.h"
#include "cantera/thermo/SurfLatIntPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/LateralInteraction.h"
#include "cantera/IncompressibleSolid.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/flowControllers.h"


using namespace std;
using namespace Cantera;

enum RctrType {
    BATCH,
    CSTR,
    PFR_0D,
    PFR
};

map<string, RctrType> rt = {{"batch", BATCH}, 
                            {"cstr", CSTR}, 
                            {"pfr_0d", PFR_0D}, 
                            {"pfr", PFR}}; 

shared_ptr<Reactor> get_0d_reactor(YAML::Node& tube_node,
                                      IdealGasMix& gas, 
                                      vector<shared_ptr<InterfaceInteractions>> surfaces)
{
    auto rctr_node = tube_node["reactor"];
    //auto rctr_type_node = rctr_node["type"];
    //auto rctr_type = rt[rctr_type_node.as<string>()];

    auto rctr = make_shared<Reactor>();
    shared_ptr<Reservoir> in_rsrv (new Reservoir), exhst (new Reservoir);
    rctr->insert(gas);
    in_rsrv->insert(gas);
    exhst->insert(gas);

    // Read the reactor dimensions
    int nodes = 1;
    auto nd_node = rctr_node["nodes"];
    if (nd_node)
        if (!nd_node.IsNull())
            nodes = nd_node.as<int>();

    double rctr_vol = strSItoDbl(rctr_node["volume"].as<string>());
    if (!rctr_vol) {
        // Raise Error
        ;
    }
    rctr_vol /= nodes;

    rctr->setInitialVolume(rctr_vol);
    vector<shared_ptr<ReactorSurface>> cat_surfs;

    double cat_area = strSItoDbl(rctr_node["cat_abyv"].as<string>());
    cat_area *= rctr_vol;
    for (const auto surf: surfaces) {
        auto cat_surf = make_shared<ReactorSurface>();
        cat_surf->setArea(cat_area);
        cat_surf->setKinetics(surf.get());
        cat_surfs.push_back(cat_surf);
        rctr->addSurface(cat_surf.get());
    }

    rctr->setChemistry();
    string mode = rctr_node["mode"].as<string>();
    if (mode == "isothermal") 
        rctr->setEnergy(0);
    else
        rctr->setEnergy(1);

    auto inlet_mfc = make_shared<MassFlowController>();
    auto outlet = make_shared<PressureController>();
    
    auto inlet_node = tube_node["inlet_gas"];
    auto vel_node = inlet_node["velocity"];
    auto restime_node = inlet_node["residence_time"];
    auto mfr_node = inlet_node["mass_flow_rate"];
    double velocity{0}, residence_time{0}, mfr{0};
    if (vel_node && !vel_node.IsNull()) {
        velocity = strSItoDbl(inlet_node["velocity"].as<string>());
        //residence_time = rctr_vol/velocity;
        mfr = rctr->mass() * velocity /rctr_vol;
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
    outlet->setPressureCoeff(0.0001);
    inlet_mfc->install(*in_rsrv, *rctr);
    outlet->install(*rctr, *exhst);

    return rctr;
}

vector<double> get_times(double end_time)
{
    vector<double> times;
    auto scale = 1e-6;
    auto r_time=0;
    while (r_time < end_time) {
        for (size_t i=1; i < 10; i++){
            r_time = i*scale;
            if (r_time > end_time)
                break;
            times.push_back(r_time);
        }
        if (r_time > end_time)
            break;
        scale *= 10;
    }
    return times;
}

void print_formation_enthalpy(vector<ThermoPhase*> phases, string output_file) 
{
    vector<doublereal> hform; 
    ofstream out;
    out.open(output_file);
    out << "#Dimensionless formation enthalpies of species (H/RT)\n" << endl;
    for  (const auto phase: phases){
        hform.resize(phase->nSpecies());
        phase->getEnthalpy_RT(hform.data());
        for (size_t k = 0; k < phase->nSpecies(); k++) {
            out.width(12); 
            out << std::left << phase->speciesName(k); 
            out.width(12); 
            out << std::right << hform[k] << endl;
        }   
    }   
    out.close();
}

void print_formation_entropy(vector<ThermoPhase*> phases, string output_file) 
{
    vector<doublereal> sform; 
    ofstream out;
    out.open(output_file);
    out << "#Dimensionless formation entropies of species (S/R)\n" << endl;
    for  (const auto phase: phases){
        sform.resize(phase->nSpecies());
        phase->getEntropy_R(sform.data());
        for (size_t k = 0; k < phase->nSpecies(); k++) {
            out.width(12); 
            out << std::left << phase->speciesName(k); 
            out.width(12); 
            out << std::right << sform[k] << endl;
        }   
    }   
    out.close();
}

void print_rxn_enthalpy(vector<Kinetics*> kinetic_mgrs, doublereal T, string output_file)
{
    vector<doublereal> hrxn;
    ofstream out;
    out.open(output_file);
    out << "#Dimensionless enthalpies of reactions (H/RT)\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            hrxn.resize(size);
            mgr->getDeltaEnthalpy(hrxn.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12);
                out << std::left << hrxn[k]/(GasConstant*T) ;
                out.width(12);
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
    out.close();
}

void print_rxn_eq_consts(vector<Kinetics*> kinetic_mgrs, string output_file)
{
    vector<doublereal> kc;
    ofstream out;
    out.open(output_file);
    out << "#Equilibrium constants of reactions\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            kc.resize(size);
            mgr->getEquilibriumConstants(kc.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12);
                out << std::left << kc[k];
                out.width(12);
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
    out.close();
}


int main(int argc, char* argv[]) 
{
    if (argc < 3) {
        // TODO: Throw error
        ;
    };

    string tube_file_name {argv[1]};       // Tube drive file in YAML format
    string phase_file_name {argv[2]};      // Thermodata in either CTI/XML formats
    auto tube_node = YAML::LoadFile(tube_file_name);
    
    // Read the phase definitions
    auto phase_node = tube_node["phases"];
    string gas_phase_name = phase_node["gas"]["name"].as<string>();
    string gas_phase_X = phase_node["gas"]["initial_state"].as<string>();
    string blk_phase_name = phase_node["bulk"]["name"].as<string>();
    auto surface_nodes = phase_node["surfaces"];
    vector<string> surface_phase_names; 
    vector<string> surface_states;
    for (size_t i=0;  i < surface_nodes.size(); i++) {
        surface_phase_names.push_back(surface_nodes[i]["name"].as<string>());
        surface_states.push_back(surface_nodes[i]["initial_state"].as<string>());
    }
    for (const auto& surf_name: surface_phase_names) 
        cout << surf_name << endl;


    IdealGasMix gas(phase_file_name, gas_phase_name);
    IncompressibleSolid bulk(phase_file_name, blk_phase_name);
    vector<ThermoPhase*> gb_phases {&gas, &bulk};
    vector<shared_ptr<InterfaceInteractions>> surf_phases;
    for (size_t i=0; i < surface_phase_names.size(); i++) {
        auto surf_ph_name = surface_phase_names[i];
        surf_phases.push_back(make_shared<InterfaceInteractions>(
                    phase_file_name, 
                    surf_ph_name, 
                    gb_phases));
    }

    /* Read the reactor definition */
    auto rctr_node = tube_node["reactor"];
    if (!rctr_node) {
        // TODO: Throw error
        ;
    };

    // Set the temp and press for all phases
    auto temp = strSItoDbl(rctr_node["temperature"].as<string>());
    auto press = strSItoDbl(rctr_node["pressure"].as<string>());
    gas.setState_TPX(temp, press, gas_phase_X);
    bulk.setState_TP(temp, press);
    for (size_t i=0; i < surface_phase_names.size(); i++) {
        surf_phases[i]->setState_TPX(temp, press, surface_states[i]);
    }

    auto rctr_type_node = rctr_node["type"];
    auto rctr_type = rt[rctr_type_node.as<string>()];


    if (rctr_type == BATCH || rctr_type == CSTR || rctr_type == PFR_0D) { // 0d reactors
        auto rctr = get_0d_reactor(tube_node, gas, surf_phases);

        int nodes = 1;
        auto nd_node = rctr_node["nodes"];
        if (nd_node)
            if (!nd_node.IsNull())
                nodes = nd_node.as<int>();

        auto simul_node = tube_node["simulation"];
        auto end_time = strSItoDbl(simul_node["end_time"].as<string>());
        auto solver_node = simul_node["solver"];
        auto abs_tol = solver_node["atol"].as<double>();
        auto rel_tol = solver_node["rtol"].as<double>();

        ReactorNet rnet;
        rnet.addReactor(*rctr);
        rnet.setTolerances(rel_tol, abs_tol);

        //TODO: Implement the advancing function
        //advance_0d_reactor(rnet, nodes, end_time);
    }
    else if (rctr_type == PFR) { // 1d reactor
        ; //TODO: Add 1d PFR implementation
    }
    
    /*
    switch (rctr_type) {
        case BATCH:
        case CSTR:
        case PFR_0D:
            rctr_vol = strSItoDbl(rctr_node["volume"].as<string>());
            cout << rctr_vol << endl;
            if (!rctr_vol) {
                // Raise Error
                ;
            }
            break;
        case PFR:
            rctr_len = strSItoDbl(rctr_node["length"].as<string>());
            rctr_area = strSItoDbl(rctr_node["area"].as<string>());
            if (!rctr_len || !rctr_area) {
                // Raise Error
                ;
            }
            break;
    }
    */

    // Start setting the reactor
}
