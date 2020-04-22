#include <string>
#include <vector>
#include <memory>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>
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
#include "run_reactor.h"
#include "omkmexceptions.h"
#include "io.h"
#include "ReactorNetHybrid.h"


namespace OpenMKM 
{

using namespace std;
using namespace Cantera;
namespace fs = boost::filesystem;

//TODO: Add nSurfaces() function to Reactor/ReactorBase to eliminate the need
// to pass surfaces argument
/*
void print_rctr_state(double z, Reactor* rctr, vector<SurfPhase*> surfaces, 
                      ofstream& gas_mole_out, ofstream& gas_mass_out, 
                      ofstream& gas_msdot_out, ofstream& surf_cov_out,
                      ofstream& state_var_out)
*/


void run_0d_reactor(ReactorParser& rctr_parser,
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<InterfaceInteractions>> surfaces,
                    ofstream& gen_info_out) 
{
    //Define the reactor based on the input file
    auto rctr = make_shared<IdealGasTRampReactor>();
    rctr->insert(*gas);
    gen_info_out << "Reactor density " << rctr->density() << endl;

    // Read the reactor dimensions
    double rctr_vol = rctr_parser.getVolume();
    auto rctr_type = rctr_parser.getReactorType();
    size_t rctr_nos = 1;
    if (rctr_type == PFR_0D) { // Modify rctr_nos only for PFR_0D
        rctr_nos = rctr_parser.getNodes();
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

    double cat_abyv = rctr_parser.getCatalystAbyV();
    cout << "Catalyst loading (Area/Volume): " << cat_abyv << endl;

    if (cat_abyv == 0.0 && surfaces.size() > 0) {
        cout << "WARNING!!!\nCatalyst loading is zero.\n"
             << "Ignoring the surface phases given in the YAML file\n"
             << "--------------\n";
    }

    if (cat_abyv != 0.0) {
        auto cat_area = cat_abyv * rctr_vol;
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
    cout << "# of surface phases considered: " << cat_surfs.size() << endl;

    // Only CSTR and PFR_0D require valves, controllers and 
    // reservoirs to facilitate mass transfer. Attach them to reactor

    // The following are nominally defined in the case of BATCH reactor also 
    // to prevent crashing the program
    auto in_rsrv = make_shared<Reservoir>();
    auto exhst = make_shared<Reservoir>();
    auto inlet_mfc = make_shared<MassFlowController>();
    auto outlet = make_shared<PressureController>();
    
    // Set the output type
    OutputFormat data_format = rctr_parser.printFormat();
    setOutputFormat(data_format);

    if (rctr_type != BATCH) {
        bool fr_defined = rctr_parser.FlowRateDefined();
        bool mfr_defined = rctr_parser.MassFlowRateDefined();
        bool rt_defined = rctr_parser.ResidenceTimeDefined();
        int no_var_defined = fr_defined + mfr_defined + rt_defined;
        if (no_var_defined > 1) {
            cout << "WARNING: Define only one of 'flow_rate', 'mass_flow_rate', "
                 << "'residence_time'. Only one of the variables is used " 
                 << "in the order shown above" << endl;
        } else if (!no_var_defined) {
            cout << "WARNING: Define one of 'flow_rate', 'mass_flow_rate', "
                 << "'residence_time' within inlet_gas. Using 0 for flow "
                 << "rate, which is equivalent to simulating batch reactor." << endl;
        }
        double mfr{0};
        if (fr_defined) {
            auto flow_rate = rctr_parser.getFlowRate();
            mfr = rctr->mass() * flow_rate / rctr_vol;
        }
        else if (mfr_defined) {
            mfr = rctr_parser.getMassFlowRate();
        } else if (rt_defined) {
            auto rt = rctr_parser.getResidenceTime();
            mfr = rctr->mass()/rt;
        }

        inlet_mfc->setMassFlowRate(mfr);
        outlet->setMaster(inlet_mfc.get());
        outlet->setPressureCoeff(0.00001);
        in_rsrv->insert(*gas);
        exhst->insert(*gas);
        inlet_mfc->install(*in_rsrv, *rctr);
        outlet->install(*rctr, *exhst);
    }

    rctr->setChemistry();

    // Read the reactor mode and set corresponding parameters
    string mode = rctr_parser.getMode();
    cout << "Reactor temperature mode: " << mode << endl;
    gen_info_out << "Reactor Temperature mode: "  << mode << endl;

    // heat_rsrv and wall are nominally defined.
    // They are used if heat transfer is required.
    auto heat_rsrv = make_shared<Reservoir>(); 
    auto wall = make_shared<Wall>();
    if (mode == "isothermal" || mode == "tpd")  
        rctr->setEnergy(0);
    else {
        rctr->setEnergy(1);
        if (mode == "heat"){
            double htc = rctr_parser.getWallHeatTransferCoeff();  // htc
            double wall_abyv = rctr_parser.getWallSpecificArea(); // wall_abyv
            double ext_temp = rctr_parser.getExternalTemp();      // Text
            auto press = gas->pressure();
            gas->setState_TP(ext_temp, press);
            heat_rsrv->insert(*gas);
            wall->setHeatTransferCoeff(htc);
            wall->setArea(wall_abyv * rctr->volume());
            wall->install(*heat_rsrv, *rctr);
        }
    }

    /* Sensitivity Analysis */
    bool sens_on = rctr_parser.isSensitivityAnalysisEnabled();
    vector<string> rxnids;

    if (sens_on) {
        rxnids = rctr_parser.getSensitivityReactions();
        for (auto& id : rxnids){
            // Identify which kinetics the reaction belong to
            // First gas phase kinetics
            bool rxnConsumed = false;
            for (size_t i = 0; i < gas->nReactions(); i++){
                if (gas->reaction(i)->id == id){
                    rctr->addSensitivityReaction(i);
                    rxnConsumed = true;
                    break;
                }
            }
            size_t kin_ind = 0;
            while(!rxnConsumed && kin_ind < surfaces.size()) {
                for (size_t i = 0; i < surfaces[kin_ind]->nReactions(); i++){
                    if (surfaces[kin_ind]->reaction(i)->id == id){
			cat_surfs[kin_ind]->addSensitivityReaction(i);
			rxnConsumed = true;
                        break;
                    }
                }
                kin_ind++;
            }
            if(!rxnConsumed) {
                throw CanteraError("Sensitivity reaction id {} not found in any mechanism", id);
            }
        }

    }
    // Read simulation parameters
    double end_time = 0;
    if (mode == "tpd"){
        auto beg_temp = rctr_parser.T();
        auto end_temp = rctr_parser.getTPDEndTemp();
        auto temp_ramp = rctr_parser.getTPDTempRamp();
        rctr->setBeta(temp_ramp);
        end_time = (end_temp - beg_temp)/temp_ramp;
        cout << "TPD simulation" << endl;
        cout << "Temperature ramp " << temp_ramp << endl;
        cout << "End temperature  " << end_temp << endl;
    } else {
        end_time = rctr_parser.getEndTime();
    }

    vector<double> times;
    bool transient_log = rctr_parser.logTransient();
    if (transient_log && rctr_nos == 1){ // Only for singler CSTR or a batch reactor
        auto step_size = rctr_parser.getInitStep();
        cout << "Simulation Step size " << step_size << endl;

        string stepping = rctr_parser.steppingType();
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
    //ReactorNetHybrid rnet; 
    rnet.addReactor(*rctr);

    // Pass any user defined numerical options to the solver
    if (rctr_parser.tolerancesDefined()){
        auto abs_tol = rctr_parser.get_atol();
        auto rel_tol = rctr_parser.get_rtol();
        rnet.setTolerances(rel_tol, abs_tol);
    }
    /* For these options, changes have to be made in Cantera
    if (rctr_parser.solverInitStepSizeDefined()){
        rnet.setInitialStepSize(rctr_parser.getSolverInitStepSize());
    }
    if (rctr_parser.solverMaxStepsDefined()){
        rnet.setMaxNumSteps(rctr_parser.getSolverMaxSteps());
    }*/

    auto gas_print_specie_header = [&gas, data_format](string ind_var, ostream& out) -> void
    {
        if (data_format == OutputFormat::CSV) {
            out << ind_var;
            for (const auto & sp_name : gas->speciesNames()) {
                out << "," << sp_name;
            }
        } else {
            out << setw(16) << left << ind_var;
            for (const auto & sp_name : gas->speciesNames()) {
                out << setw(16) << left << sp_name;
            }
        }
        out << endl;
    };

    // Steady state condition makes sense For CSTR and PFR_0D 
    vector<double> T_params = rctr_parser.Ts();
    vector<double> P_params = rctr_parser.Ps();
    vector<double> fr_params = rctr_parser.FRs();
    if (!fr_params.size()){
        if (rctr_type != BATCH) {
            auto mfr = inlet_mfc->massFlowRate();
            fr_params.push_back(mfr/gas->density());
        }
        else {
            fr_params.push_back(0.0);
        }
    }

    fs::path curr_dir = ".";
    for (const auto& T : T_params){
        for (const auto& P : P_params){
            for (const auto& fr : fr_params){
                // Set the gas state TPX and sync the appropriate reactors
                auto X = rctr_parser.getGasPhaseComposition();
                gas->setState_TPX(T, P, X);
                if (rctr_type != BATCH){
                    in_rsrv->syncState();
                    exhst->syncState();
                    inlet_mfc->setMassFlowRate(fr * gas->density());
                }
                rctr->syncState();
                rnet.setInitialTime(0);
                rnet.reinitialize();

                string new_dir = "T-"  + to_string(T) + ",P-" + to_string(P) 
                                       + ",fr-" + to_string(fr);
                fs::path out_dir = curr_dir;
                if (T_params.size() > 1 || P_params.size() > 1 || fr_params.size() > 1){
                    out_dir /= new_dir;
                    create_directory(out_dir);
                }
                
                string file_ext;
                if (data_format == OutputFormat::CSV) {
                    file_ext = "csv";
                } else {
                    file_ext = "dat";
                }

                ofstream gas_ss_mole_out ((out_dir / ("gas_mole_ss." + file_ext)).string(), ios::out);
                ofstream gas_ss_mass_out ((out_dir / ("gas_mass_ss." + file_ext)).string(), ios::out);
                ofstream gas_ss_msdot_out ((out_dir / ("gas_msdot_ss." + file_ext)).string(), ios::out);
                ofstream surf_ss_out ((out_dir / ("surf_cov_ss." + file_ext)).string(), ios::out);
                ofstream state_var_ss_out ((out_dir / ("rctr_state_ss." + file_ext)).string(), ios::out);
                ofstream rates_out ((out_dir / "rates_ss.out").string(), ios::out);

                if (data_format == OutputFormat::DAT) {
                    gas_ss_mole_out << "#Gas Mole fractions at Steady State\n";
                    gas_ss_mass_out << "#Gas Mass fractions at Steady State\n";
                    gas_ss_msdot_out << "#Surface Production Rates of  Gas Species at Steady State\n";
                    surf_ss_out << "#Surace Coverages at Steady State\n"; 
                    state_var_ss_out << "#Steady State Reactor State\n";
                }
                
                // Reaction path analysis (RPA) data, consisting of rates of progress.
                // By default these values are written for simulation end time (PFR or CSTR)
                // and for PFR exit. If enabled, RPA data is also written for all z-points
                // of PFR
                auto rpa_flag = rctr_parser.RPA();

                auto state_var_print_hdr = [data_format](ostream& out, const string ind0) -> void
                {
                    if (data_format == OutputFormat::CSV) {
                        out << ind0 << ","//"z(m)"
                            << "Temperature(K)" << "," 
                            << "Pressure(Pa)" << "," 
                            << "Density(kg/m3)" << "," 
                            << "U(J/kg)" 
                            << endl;
                    } else {
                        out << setw(16) << left << ind0 //"z(m)"
                            << setw(16) << left << "Temperature(K)" 
                            << setw(16) << left << "Pressure(Pa)" 
                            << setw(16) << left << "Density(kg/m3)" 
                            << setw(16) << left << "U(J/kg)" 
                            << endl;
                    }
                };

                auto surf_print_hdr = [&](ostream& out, const string ind0) -> void
                {
                    if (data_format == OutputFormat::CSV) {
                        out << ind0;
                        for (const auto surf : surfaces) {
                            for (const auto & sp_name : surf->speciesNames()) {
                                out << "," << sp_name;
                            }
                        }
                    } else {
                        out << setw(16) << left << ind0;
                        for (const auto surf : surfaces) {
                            for (const auto & sp_name : surf->speciesNames()) {
                                out << setw(16) << left << sp_name;
                            }
                        }
                    }
                    out << endl;
                };

                if (rctr_type == PFR_0D) {
                    gas_print_specie_header("z(m)", gas_ss_mole_out);
                    gas_print_specie_header("z(m)", gas_ss_mass_out);
                    gas_print_specie_header("z(m)", gas_ss_msdot_out);
                    surf_print_hdr(surf_ss_out, "z(m)");
                    state_var_print_hdr(state_var_ss_out, "z(m)");
                } else {
                    gas_print_specie_header("t(s)", gas_ss_mole_out);
                    gas_print_specie_header("t(s)", gas_ss_mass_out);
                    gas_print_specie_header("t(s)", gas_ss_msdot_out);
                    surf_print_hdr(surf_ss_out, "t(s)");
                    state_var_print_hdr(state_var_ss_out, "t(s)");
                }
                    /*
                    surf_ss_out 
                        << setw(16) << left << "z(m)";
                    for (const auto surf : surfaces) {
                        for (const auto & sp_name : surf->speciesNames()) {
                            surf_ss_out << setw(16) << left << sp_name;
                        }
                    }
                    surf_ss_out << endl;
                    */
                    
                    // Print the inlet state
                    //rnet.reinitialize();
                if (rctr_type != BATCH) {
                    print_0d_rctr_state(0, rctr.get(), surf_phases, gas_ss_mole_out, 
                                        gas_ss_mass_out, gas_ss_msdot_out, surf_ss_out, 
                                        state_var_ss_out);
                }

                // Transient state makes sense For BATCH and CSTR 
                ofstream gas_tr_mole_out ((out_dir / ("gas_mole_tr." + file_ext)).string(), ios::out);
                ofstream gas_tr_mass_out ((out_dir / ("gas_mass_tr." + file_ext)).string(), ios::out);
                ofstream gas_tr_msdot_out ((out_dir / ("gas_msdot_tr." + file_ext)).string(), ios::out);
                ofstream surf_tr_out ((out_dir / ("surf_cov_tr." + file_ext)).string(), ios::out);
                ofstream state_var_tr_out ((out_dir / ("rctr_state_tr." + file_ext)).string(), ios::out);
                if (rctr_type != PFR_0D && times.size() > 1) {
                    if (data_format == OutputFormat::DAT) {
                        gas_tr_mole_out << "Transient Gas Mole fractions"  << endl;
                        gas_tr_mass_out << "Transient Gas Mass fractions"  << endl;
                        gas_tr_msdot_out << "Transient Surface Production Rates of  Gas Species"  << endl;
                        surf_tr_out << "Transient Surface Coverages"  << endl;
                    }
                    gas_print_specie_header("t(s)", gas_tr_mole_out);
                    gas_print_specie_header("t(s)", gas_tr_mass_out);
                    gas_print_specie_header("t(s)", gas_tr_msdot_out);
                    surf_print_hdr(surf_tr_out, "t(s)");
                    state_var_print_hdr(state_var_tr_out, "t(s)");
                }

                vector<double> gas_X(gas->nSpecies()); // Temporary work array
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

                    /*
                    if (times.size() == 1) { 
                        cout << "Solving with steady state solver" << endl;
                        rnet.setIntegratorEndTime(times[0]);
                        rnet.solve();
                    }
                    else*/ {
                        for (const auto & tm : times) {
                            rnet.advance(tm);
                            if (transient_log) {
                                print_0d_rctr_state(tm, rctr.get(), surf_phases, 
                                                    gas_tr_mole_out, gas_tr_mass_out, 
                                                    gas_tr_msdot_out, surf_tr_out, 
                                                    state_var_tr_out);

                            }
                        }
                    }
                    rctr->restoreState();

                    double ind_val;
                    if (rctr_type == PFR_0D) {
                        ind_val = (i+0.5)*rctr_vol;
                    }
                    else {
                        ind_val = times.back();
                    }
                        
                    print_0d_rctr_state(ind_val, rctr.get(), surf_phases, 
                                        gas_ss_mole_out, gas_ss_mass_out, 
                                        gas_ss_msdot_out, surf_ss_out, state_var_ss_out);

                    if (rpa_flag) {
                        string rpa_file_name = "rates_z-";
                        rpa_file_name += to_string((i+0.5)*rctr_vol);
                        rpa_file_name += ".out";
                        ofstream rates_out((out_dir / rpa_file_name).string(), 
                                           ios::out); // Masks the name
                        print_rxn_rates_hdr(rates_out);

                        rates_out.precision(6);

                        print_rxn_rates(gas.get(), rates_out);
                        for (auto surf : surfaces) {
                            print_rxn_rates(surf.get(), rates_out);
                        }
                    }
                }
                // Print final rpa data
                print_rxn_rates_hdr(rates_out);
                rates_out.precision(6);
                print_rxn_rates(gas.get(), rates_out);
                for (auto surf : surfaces) {
                    print_rxn_rates(surf.get(), rates_out);
                }

            }
        }
    }
}
}
