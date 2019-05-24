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
#include "run_reactor.h"
#include "pfr1d.h"
#include "pfr1d_solver.h"
#include "io.h"


using namespace std;
using namespace Cantera;

<<<<<<< HEAD
namespace OpenMKM 
{

void run_1d_reactor(ReactorParser& rctr_parser,
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<InterfaceInteractions>> surfaces,
                    ofstream& gen_info)
{
    //Define the reactor based on the input file
    //auto rctr_node = tube_node["reactor"];

    // Read the reactor dimensions
    double rctr_xc_area = rctr_parser.getXCArea();
    double rctr_len = rctr_parser.getLength();

    bool fr_defined = rctr_parser.FlowRateDefined();
    bool mfr_defined = rctr_parser.MassFlowRateDefined();
    bool rt_defined = rctr_parser.ResidenceTimeDefined();
    int no_var_defined = fr_defined + mfr_defined + rt_defined;
    if (no_var_defined > 1) {
        cout << "Define only one of 'flow_rate', 'mass_flow_rate', 'residence_time'"
             << "Only one of the variables is arbitrarily used" << endl;
    }
    double velocity{0};
    if (fr_defined) {
        auto flow_rate = rctr_parser.getFlowRate();
        velocity = flow_rate / rctr_xc_area;
    }
    else if (mfr_defined) {
        auto mfr = rctr_parser.getMassFlowRate();
        auto flow_rate = mfr / gas->density();
        velocity = flow_rate / rctr_xc_area;
    } else if (rt_defined) {
        auto rt = rctr_parser.getResidenceTime();
        velocity = rctr_len / rt;
    }

    double cat_abyv = 0;
    if(rctr_parser.catalystAreaDefined()) {
        cat_abyv = rctr_parser.getCatalystAbyV();
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

    auto pfr = PFR1d(gas.get(), ikin, surf_ph, rctr_xc_area, cat_abyv, velocity);
    
    //string mode = rctr_node["mode"].as<string>();
    string mode = rctr_parser.getMode();
    cout << "mode " << mode << endl;
    if (mode == "isothermal") {
        pfr.setEnergy(0);
    }
    else if (mode == "tprofile") {
        pfr.setEnergy(0);
        /*
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
        */
        pfr.setTProfile(rctr_parser.getTProfile());
    } 
    else {
        pfr.setEnergy(1);
        if (mode == "heat") {
            double htc = rctr_parser.getWallHeatTransferCoeff();  // htc 
            double wall_abyv = rctr_parser.getWallSpecificArea(); // wall_abyv
            double ext_temp = rctr_parser.getExternalTemp();      // Text
            pfr.setHeatTransfer(htc, ext_temp, wall_abyv);  
        }
        pfr.reinit();
    }
    pfr.setConstraints();
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

    //auto simul_node = tube_node["simulation"];
    if (rctr_parser.tolerancesDefined()){
        auto abs_tol = rctr_parser.get_atol();
        auto rel_tol = rctr_parser.get_rtol();
        pfr_solver.setTolerances(rel_tol, abs_tol);
    }
    if (rctr_parser.solverInitStepSizeDefined()){
        pfr_solver.setInitialStepSize(rctr_parser.getSolverInitStepSize());
    }
    if (rctr_parser.solverMaxStepsDefined()){
        pfr_solver.setMaxNumSteps(rctr_parser.getSolverMaxSteps());
    }

    double simul_init_step = 1e-6;
    if (rctr_parser.initStepDefined()){
        simul_init_step = rctr_parser.getInitStep();
    }
    vector<double> zvals = get_log10_intervals(rctr_len, simul_init_step); //Use the same function to get z steps

    ofstream gas_mole_out("gas_mole_ss.out");
    ofstream gas_mass_out("gas_mass_ss.out");
    ofstream gas_msdot_out("gas_msdot_ss.out");
    ofstream surf_cov_out("surf_cov_ss.out");
    ofstream state_var_out("rctr_state_ss.out");
    ofstream rates_out("rates_ss.out", ios::out);

    gas_mole_out.precision(6);
    gas_mass_out.precision(6);
    gas_msdot_out.precision(6);
    surf_cov_out.precision(6);
    state_var_out.precision(6);
    rates_out.precision(6);

    auto rpa_flag = rctr_parser.RPA();

    for (const auto& z : zvals) {
        pfr_solver.solve(z);
        print_pfr_rctr_state(z, &pfr, surf_ph, gas_mole_out, gas_mass_out, 
                             gas_msdot_out, surf_cov_out, state_var_out);
        if (rpa_flag) {
            string rpa_file_name = "rates_z-";
            rpa_file_name += to_string(z);
            rpa_file_name += ".out";
            ofstream rates_out (rpa_file_name, ios::out); // Masks the name
            print_rxn_rates_hdr("Rates (mol/s) and Partial Equilibrium Analysis:",
                                rates_out);
            rates_out.precision(6);

            auto rxn_index = 0;
            print_rxn_rates(gas.get(), rxn_index, rates_out);
            rxn_index += gas->nReactions();
            for (auto surf : surfaces) {
                print_rxn_rates(surf.get(), rxn_index, rates_out);
                rxn_index += surf->nReactions();
            }
        }
    }

    
    pfr_solver.writeStateData("1d_pfr_state.out");
    pfr_solver.writeGasData("1d_pfr_gas.out");
    pfr_solver.writeSurfaceData("1d_pfr_surface.out");

    // Print final rpa data
    rates_out.precision(6);
    auto rxn_index = 0;
    print_rxn_rates(gas.get(), rxn_index, rates_out);
    rxn_index += gas->nReactions();
    for (auto surf : surfaces) {
        print_rxn_rates(surf.get(), rxn_index, rates_out);
        rxn_index += surf->nReactions();
    }


}

}
