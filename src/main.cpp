#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
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
        string err_str("ERROR!!!!\n---------\nInsufficient number of arguments.\n");
        err_str += "OpenMKM requires one YAML file specifying simulation parameters\n";
        err_str += "as first argument and one XML file containing thermodynamic \n";
        err_str += "definitions in Cantera format as second argument.\n---------\n";
        //throw Cantera::CanteraError("main", err_str);
        cerr << err_str;
        cout << "Exiting the program due to insufficient number of arguments" << endl;
        return 1;
    };

    auto start_t = high_resolution_clock::now();

    string tube_file_name {argv[1]};       // Tube drive file in YAML format
    string phase_filename {argv[2]};      // Thermodata in either CTI/XML formats
    
    try{
        // Read the gas phase definition
        auto rctr_parser = ReactorParser(tube_file_name);
        shared_ptr<Solution> gas = rctr_parser.getGasSolution(phase_filename);
        cout << "Kinetics type " << gas->kinetics()->kineticsType() << endl;
        vector<shared_ptr<ThermoPhase>> all_phases {gas->thermo()};
        vector<Kinetics*> all_km {gas->kinetics().get()};
        vector<shared_ptr<Solution>> gb_solns {gas};

        // Try to read the bulk node 
        bool blk_phase_defined = rctr_parser.bulkPhaseDefined(phase_filename);
        if (blk_phase_defined) {
            //auto bulk = rctr_parser.getBulkPhase(phase_filename);
            auto bulk = rctr_parser.getBulkSolution(phase_filename);
            all_phases.push_back(bulk->thermo());
            gb_solns.push_back(bulk);
        }
        bool surf_phases_defined = rctr_parser.surfacePhasesDefined(phase_filename);
        //vector<shared_ptr<InterfaceInteractions>> surf_phases;
        //if (surf_phases_defined && blk_phase_defined) {
         //   surf_phases = rctr_parser.getSurfPhases(phase_filename, gb_phases);
          //  for (auto& surf_phase : surf_phases) {
           //     all_phases.push_back(surf_phase);
            //    all_km.push_back(surf_phase.get());
            //}
        //}
        vector<shared_ptr<Solution>> surf_solns;
        if (surf_phases_defined) {
            surf_solns = rctr_parser.getSurfaceSolutions(phase_filename, gb_solns);
            for (auto& surf_soln : surf_solns) {
                cout << "Kinetics type " << surf_soln->kinetics()->kineticsType() << endl;
                all_phases.push_back(surf_soln->thermo());
                all_km.push_back(surf_soln->kinetics().get());
            }
        }
        cout << "Total # of phases: " << all_phases.size() << endl;
        cout << "Surface phase defined? " << boolalpha << surf_phases_defined << endl;

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
        print_rxn_enthalpy(all_km, gas->thermo()->temperature(), "Hrxn.out");
        print_rxn_entropy(all_km, "Srxn.out");
        print_rxn_gibbs(all_km, gas->thermo()->temperature(), "Grxn.out");
        print_rxn_kf(all_km,  "kf.out");
        print_rxn_kc(all_km,  "kc.out");
        print_rxn_kr(all_km,  "kr.out");


        auto rctr_type = rctr_parser.getReactorType();
        cout << "Reactor Model: " << RctrTypeString[rctr_type] << endl;
        if (rctr_type == BATCH || rctr_type == CSTR || rctr_type == PFR_0D) { // 0d reactors
            //run_0d_reactor(rctr_parser, gas, surf_phases, gen_info);
            run_0d_reactor(rctr_parser, gas, surf_solns, gen_info);

        }
        else if (rctr_type == PFR) { // 1d reactor
            //run_1d_reactor(rctr_parser, gas, surf_phases, gen_info);
            run_1d_reactor(rctr_parser, gas, surf_solns, gen_info);
        }
        print_species(all_phases, "species.out");
    }
    catch (YAMLParserError& err) {
        cout << "Error in parsing YAML file. Exiting..." << endl;
        gen_info << "Error in parsing YAML file. Exiting..." << endl;
        cerr << err.what() << endl;
        return -3;
    }
    catch (CanteraError& err) {
        cout << "Error raised in Cantera. Exiting..." << endl;
        gen_info << "Error raised in Cantera. Exiting... " << endl;
        cerr << err.what() << endl;
        return -2;
    }
    catch (exception& e){
        cout << "Unknown Error. Exiting..." << endl;
        gen_info << "Unknown Error" << endl;
        cerr << e.what() << endl;
        return -1;
    }
 
    auto end_t = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end_t - start_t);
    cout << "Program ran for " << duration.count() <<  " milliseconds" << endl;
    gen_info << "Program ran for " << duration.count() << " milliseconds" << endl;
    return 0;
}
