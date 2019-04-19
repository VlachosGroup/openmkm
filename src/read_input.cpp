#include <string>
#include <vector>
/*
#include <memory>
#include <map>
#include <iostream>
#include <fstream>

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

#include "reactor.h"

using namespace std;
using namespace Cantera;
using namespace HeteroCt;



/*


void advance_0d_reactor(shared_ptr<ReactorNet> rnet, int nodes, double end_time)
{
    cout << "Reached advance function" << endl;
    auto times = get_times(end_time);
    cout << "time 1 "  << times[0] << endl;
    cout << "rctr density " << rnet->reactor(0).density() << endl;
    //rnet->step();
    for (const auto & tm : times) {
        if (tm < 2e-6) {
            rnet->advance(tm);//(tm);
        }
    }

}
*/

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


    auto gas = make_shared<IdealGasMix>(phase_file_name, gas_phase_name);
    auto bulk = make_shared<StoichSubstance>(phase_file_name, blk_phase_name);
    vector<ThermoPhase*> gb_phases {gas.get(), bulk.get()};
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
    gas->setState_TPX(temp, press, gas_phase_X);
    bulk->setState_TP(temp, press);
    for (size_t i=0; i < surface_phase_names.size(); i++) {
        surf_phases[i]->setState_TP(temp, press);
        surf_phases[i]->setCoveragesByName(surface_states[i]);
    }

    auto rctr_type_node = rctr_node["type"];
    auto rctr_type = HeteroCt::RctrTypeMap[rctr_type_node.as<string>()];

    if (rctr_type == BATCH || rctr_type == CSTR || rctr_type == PFR_0D) { // 0d reactors
        run_0d_reactor(tube_node, gas, surf_phases);

        /*int nodes = 1;
        auto nd_node = rctr_node["nodes"];
        if (nd_node)
            if (!nd_node.IsNull())
                nodes = nd_node.as<int>();

                
        auto rnet = make_shared<ReactorNet>();
        rnet->addReactor(*rctr);
        cout << "rel tol: " << rel_tol << "  abs tol: " << abs_tol << endl;
        rnet->setTolerances(rel_tol, abs_tol);

        //TODO: Implement the advancing function
        cout << "Reached here" << endl;
        auto times = get_times(end_time);
        cout << "time 1 "  << times[0] << endl;
        cout << "rctr density " << rnet->reactor(0).density() << endl;
        //rnet->step();
        for (const auto & tm : times) {
            if (tm < 2e-6) {
                rnet->advance(tm);//(tm);
            }
        }

        advance_0d_reactor(rnet, nodes, end_time);

        */
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
