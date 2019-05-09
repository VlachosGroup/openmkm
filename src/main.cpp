#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
/*
#include <memory>
#include <map>

include <yaml-cpp/yaml.h>
#include "cantera/base/ct_defs.h"
#include "cantera/thermo/SurfLatIntPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/LateralInteraction.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/Reservoir.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/zeroD/flowControllers.h"
*/

#include "cantera/base/stringUtils.h"
#include "cantera/IdealGasMix.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/InterfaceLatInt.h"
#include "cantera/Interface.h"

#include "io.h"
#include "reactor.h"
#include "hctexceptions.h"

using namespace std;
using namespace std::chrono;
using namespace Cantera;
using namespace HeteroCt;



std::map<std::string, RctrType> RctrTypeMap = {{"batch", BATCH}, 
                                               {"cstr", CSTR}, 
                                               {"pfr_0d", PFR_0D}, 
                                               {"pfr", PFR}}; 


int main(int argc, char* argv[]) 
{
    ofstream gen_info ("general_info.out", ios::out);
    print_htrct_header(gen_info);
    print_htrct_header(cout);
    if (argc < 3) {
        // TODO: Throw error
        ;
    };

    auto start_t = high_resolution_clock::now();

    string tube_file_name {argv[1]};       // Tube drive file in YAML format
    string phase_file_name {argv[2]};      // Thermodata in either CTI/XML formats
    auto tube_node = YAML::LoadFile(tube_file_name);
    
    // Read the gas phase definition
    auto phase_node = tube_node["phases"];
    string gas_phase_name = phase_node["gas"]["name"].as<string>();
    auto gas = make_shared<IdealGasMix>(phase_file_name, gas_phase_name);
    vector<shared_ptr<ThermoPhase>> all_phases {gas};
    vector<Kinetics*> all_km {gas.get()};
    string gas_phase_X = phase_node["gas"]["initial_state"].as<string>();


    /* Read the state variables */
    auto rctr_node = tube_node["reactor"];
    if (!rctr_node) {
        throw YAMLParserError("main.cpp::main", "reactor", "Node not found");
    };

    // Set the temp and press for all phases
    auto temp = strSItoDbl(rctr_node["temperature"].as<string>());
    auto press = strSItoDbl(rctr_node["pressure"].as<string>());
    gas->setState_TPX(temp, press, gas_phase_X);

    // Try to read the bulk node and if present read surface definitons as well
    //vector<shared_ptr<Interface>> surf_phases;
    vector<shared_ptr<InterfaceInteractions>> surf_phases;
    auto bulk_node = phase_node["bulk"];
    string blk_phase_name;
    if (bulk_node && !bulk_node.IsNull()) { 

        blk_phase_name = phase_node["bulk"]["name"].as<string>();
        auto bulk = make_shared<StoichSubstance>(phase_file_name, blk_phase_name);
        bulk->setState_TP(temp, press);  // Set bulk state
        vector<ThermoPhase*> gb_phases {gas.get(), bulk.get()};
        all_phases.push_back(bulk);

        vector<string> surface_phase_names; 
        vector<string> surface_states;
        auto surface_nodes = phase_node["surfaces"];
        if (surface_nodes && !surface_nodes.IsNull()) {
            for (size_t i=0;  i < surface_nodes.size(); i++) {
                surface_phase_names.push_back(surface_nodes[i]["name"].as<string>());
                surface_states.push_back(surface_nodes[i]["initial_state"].as<string>());
            }
            for (const auto& surf_name: surface_phase_names) 
                cout << surf_name << endl;
        }

        for (size_t i=0; i < surface_phase_names.size(); i++) {
            auto surf_ph_name = surface_phase_names[i];
            auto surf = make_shared<InterfaceInteractions>(phase_file_name, 
            //auto surf = make_shared<Interface>(phase_file_name, 
                                               surf_ph_name, 
                                               gb_phases);
            surf_phases.push_back(surf);
            all_km.push_back(surf.get());
            all_phases.push_back(surf);
        }

        for (size_t i=0; i < surface_phase_names.size(); i++) {
            surf_phases[i]->setState_TP(temp, press);
            surf_phases[i]->setCoveragesByName(surface_states[i]);
        }
    }

    /* Print the species thermodynamic info */
    print_formation_enthalpy(all_phases, "Hform.out");
    print_formation_entropy(all_phases, "Sform.out");

    /* Print the reaction thermodynamic info */
    size_t n_rxns = 0;
    for (const auto km : all_km) {
        n_rxns += km->nReactions();
    }
    cout << "Total # of reactions: " << n_rxns << endl;

    print_rxn_enthalpy(all_km, gas->temperature(), "Hrxn.out");
    print_rxn_entropy(all_km, "Srxn.out");
    print_rxn_gibbs(all_km, gas->temperature(), "Grxn.out");
    print_rxn_kf(all_km,  "kf.out");
    print_rxn_kc(all_km,  "kc.out");
    print_rxn_kr(all_km,  "kr.out");


    auto rctr_type_node = rctr_node["type"];
    auto rctr_type = RctrTypeMap[rctr_type_node.as<string>()];

    if (rctr_type == BATCH || rctr_type == CSTR || rctr_type == PFR_0D) { // 0d reactors
        run_0d_reactor(rctr_type, tube_node, gas, surf_phases, gen_info);

    }
    else if (rctr_type == PFR) { // 1d reactor
        run_1d_reactor(tube_node, gas, surf_phases, gen_info);
    }
 
    auto end_t = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_t - start_t);
    cout << "Program ran for " << duration.count() <<  " milliseconds" << endl;
    gen_info << "Program ran for " << duration.count() << " milliseconds" << endl;

   
}
