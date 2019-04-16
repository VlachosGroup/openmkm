#include <string>
#include <memory>
#include <map>
#include <iostream>
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

shared_ptr<ReactorNet> get_0d_reactor(RctrType rctr_type, YAML::Node& tube_node,
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
    if (rctr_type == PFR_0D) // Update nodes
        if (rctr_node["nodes"])
            nodes = rctr_node["nodes"].as<int>();

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
    if (vel_node){
        velocity = strSItoDbl(inlet_node["velocity"].as<string>());
        //residence_time = rctr_vol/velocity;
        mfr = rctr->mass() * velocity /rctr_vol;
    }
    else if (restime_node) {
        residence_time = strSItoDbl(inlet_node["residence_time"].as<string>());
        mfr = rctr->mass()/residence_time;
    }
    else {
        mfr = strSItoDbl(inlet_node["mass_flow_rate"].as<string>());
    }
    inlet_mfc->setMassFlowRate(mfr);
    outlet->setMaster(inlet_mfc.get());
    outlet->setPressureCoeff(0.0001);

    auto rnet = make_shared<ReactorNet>();
    rnet->addReactor(*rctr);
    return rnet;
            
}

int main(int argc, char* argv[]) 
{
    if (argc < 3) {
        // Throw error
        ;
    };

    string tube_file_name {argv[1]};       // Tube drive file in YAML format
    string phase_file_name {argv[2]};      // Thermodata in either CTI/XML formats
    auto tube_node = YAML::LoadFile(tube_file_name);
    
    // Read the phase definitions
    auto phase_node = tube_node["phases"];
    cout << "Reached till here1" << endl;
    string gas_phase_name = phase_node["gas"]["name"].as<string>();
    string gas_phase_X = phase_node["gas"]["initial_state"].as<string>();
    string blk_phase_name = phase_node["bulk"]["name"].as<string>();
    cout << "Reached till here2" << endl;
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
        surf_phases.push_back(
                make_shared<InterfaceInteractions>(phase_file_name, 
                                                   surf_ph_name, 
                                                   gb_phases));
    }

    /* Read the reactor definition */
    auto rctr_node = tube_node["reactor"];
    if (!rctr_node) {
        // Throw error
        ;
    };

    cout << "Reached till here" << endl;
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

    cout << "Reached till here" << endl;
    if (rctr_type == BATCH || rctr_type == CSTR || rctr_type == PFR_0D) { // 0d reactors
        shared_ptr<ReactorNet> rnet = get_0d_reactor(rctr_type, tube_node, 
                                                     gas, surf_phases);
    }
    else if (rctr_type == PFR) { // 1d reactor
        ;
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
