#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>

#include "cantera/base/stringUtils.h"
#include "cantera/IdealGasMix.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/InterfaceLatInt.h"
#include "cantera/Interface.h"

#include "io.h"
#include "run_reactor.h"
#include "omkmexceptions.h"
#include "reactor_parser.h"

using namespace std;
using namespace std::chrono;
using namespace Cantera;
using namespace OpenMKM;

map<RctrType, std::string> RctrTypeString  = { 
    {BATCH, "Batch Reactor"}, 
    {CSTR, "Continous Stirred Tank Reactor (CSTR)"}, 
    {PFR_0D, "PFR as series of CSTRs"}, 
    {PFR, "Plug Flow Reactor (PFR)"}
};


int main(int argc, char* argv[]) 
{
    ofstream gen_info ("general_info.out", ios::out);
    print_omkm_header(gen_info);
    print_omkm_header(cout);
    if (argc < 3) {
        string err_str("Insufficient number of arguments.\n");
        err_str += "OpenMKM requires one YAML file specifying simulation parameters\n";
        err_str += "as first argument and one XML file containing thermodynamic \n";
        err_str += "definitions in Cantera format as second argument.\n";
        throw Cantera::CanteraError("main", err_str);
    };

    auto start_t = high_resolution_clock::now();

    string tube_file_name {argv[1]};       // Tube drive file in YAML format
    string phase_filename {argv[2]};      // Thermodata in either CTI/XML formats
    
    // Read the gas phase definition
    auto rctr_parser = ReactorParser(tube_file_name);
    auto gas = rctr_parser.getGasPhase(phase_filename);
    vector<shared_ptr<ThermoPhase>> all_phases {gas};
    vector<Kinetics*> all_km {gas.get()};

    // Try to read the bulk node and if present read surface definitons as well
    bool blk_phase_defined = rctr_parser.bulkPhaseDefined(phase_filename);
    vector<ThermoPhase*> gb_phases;
    if (blk_phase_defined) {
        auto bulk = rctr_parser.getBulkPhase(phase_filename);
        all_phases.push_back(bulk);
        gb_phases.push_back(gas.get()); 
        gb_phases.push_back(bulk.get());
    }
    bool surf_phases_defined = rctr_parser.surfacePhasesDefined(phase_filename);
    vector<shared_ptr<InterfaceInteractions>> surf_phases;
    if (surf_phases_defined && blk_phase_defined) {
        surf_phases = rctr_parser.getSurfPhases(phase_filename, gb_phases);
        for (auto& surf_phase : surf_phases) {
            all_phases.push_back(surf_phase);
            all_km.push_back(surf_phase.get());
        }
    }

    /* Print the species thermodynamic info */
    print_species_number(all_phases);

    print_formation_enthalpy(all_phases, "Hform.out");
    print_formation_entropy(all_phases, "Sform.out");

    /* Print the reaction thermodynamic info */
    size_t n_rxns = 0;
    for (const auto km : all_km) {
        n_rxns += km->nReactions();
    }
    cout << "Total # of reactions: " << n_rxns << endl;

    print_rxns(all_km, "reactions.out");
    print_rxn_enthalpy(all_km, gas->temperature(), "Hrxn.out");
    print_rxn_entropy(all_km, "Srxn.out");
    print_rxn_gibbs(all_km, gas->temperature(), "Grxn.out");
    print_rxn_kf(all_km,  "kf.out");
    print_rxn_kc(all_km,  "kc.out");
    print_rxn_kr(all_km,  "kr.out");


    auto rctr_type = rctr_parser.getReactorType();
    cout << "Type of reactor: " << RctrTypeString[rctr_type] << endl;
    if (rctr_type == BATCH || rctr_type == CSTR || rctr_type == PFR_0D) { // 0d reactors
        run_0d_reactor(rctr_parser, gas, surf_phases, gen_info);

    }
    else if (rctr_type == PFR) { // 1d reactor
        run_1d_reactor(rctr_parser, gas, surf_phases, gen_info);
    }
    print_species(all_phases, "species.out");
 
    auto end_t = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_t - start_t);
    cout << "Program ran for " << duration.count() <<  " milliseconds" << endl;
    gen_info << "Program ran for " << duration.count() << " milliseconds" << endl;

}
