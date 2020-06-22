#include <string>
#include <vector>
#include <memory>
#include <map>
#include <iostream>
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
namespace fs = boost::filesystem;

namespace OpenMKM 
{

void run_1d_reactor(ReactorParser& rctr_parser,
                    shared_ptr<Solution> gas, 
                    vector<shared_ptr<Solution>>& surfaces,
                    ofstream& gen_info)
{
    //Define the reactor based on the input file

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
    if (mfr_defined) {
        auto mfr = rctr_parser.getMassFlowRate();
        auto flow_rate = mfr / gas->thermo()->density();
        velocity = flow_rate / rctr_xc_area;
    } else if (fr_defined) {
        auto flow_rate = rctr_parser.getFlowRate();
        velocity = flow_rate / rctr_xc_area;
    } else if (rt_defined) {
        auto rt = rctr_parser.getResidenceTime();
        velocity = rctr_len / rt;
    }

    double cat_abyv = 0;
    if(rctr_parser.catalystAreaDefined()) {
        cat_abyv = rctr_parser.getCatalystAbyV();
    }
    cout << "Catalyst loading (Area/Volume): " << cat_abyv << endl;
    if (cat_abyv == 0.0 && surfaces.size() > 0) {
        cout << "WARNING!!!\nCatalyst loading is zero.\n"
             << "Ignoring the surface phases given in the YAML file\n"
             << "--------------\n";
    }

    vector<InterfaceKinetics*> ikin;
    vector<SurfPhase*> surf_ph;
    for (const auto surf_soln: surfaces) {
        ikin.push_back(dynamic_cast<InterfaceKinetics*>(surf_soln->kinetics().get()));
        surf_ph.push_back(dynamic_cast<SurfPhase*> (surf_soln->thermo().get()));
    }

    // Before simulation save the initial coverages of surface species for reuse
    //vector<vector<double>> surf_init_covs;
    vector<double> cov;
    for (const auto surf: surf_ph) {
        cov.resize(surf->nSpecies());
        surf->getCoverages(cov.data());
        //surf_init_covs.push_back(cov);
    }

    // Start the simulation
    gen_info << "Solving for equilibirum surface coverages at PFR inlet" << endl;
    for (size_t i = 0; i < surfaces.size(); i++) {
        cout << "Surface Site Density " << surf_ph[i]->siteDensity() << endl;
        ikin[i]->solvePseudoSteadyStateProblem();
        vector<double> cov(surf_ph[i]->nSpecies());
        surf_ph[i]->getCoverages(cov.data());

        gen_info << "Equilibrium surface coverages on Surface: " 
                 <<  surf_ph[i]->name() << endl;
        for (auto j = 0; j < surf_ph[i]->nSpecies(); j++)
            gen_info << surf_ph[i]->speciesSPName(j) << " coverage: " << cov[j] << endl;
    }

    auto pfr = PFR1d(gas.get(), ikin, surf_ph, rctr_xc_area, cat_abyv, velocity);
    
    //string mode = rctr_node["mode"].as<string>();
    string mode = rctr_parser.getMode();
    cout << "Reactor temperature mode: " << mode << endl;
    gen_info << "Reactor temperature mode: "  << mode << endl;
    if (mode == "isothermal") {
        pfr.setEnergy(0);
    }
    else if (mode == "tprofile") {
        pfr.setEnergy(0);
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
    gen_info << "Energy enabled? "  << pfr.energyEnabled() << endl;

    // Read the sensitivity coefficients 
    bool sens_on = rctr_parser.isSensitivityAnalysisEnabled();
    bool full_sens = rctr_parser.isfullSensitivityAnalysis();
    vector<std::string> sens_ids;
    int nquad;
    if (sens_on) {
        if (!full_sens){ 
            // Read the sensitivity equations and enable them
            sens_ids = rctr_parser.getSensitivityReactions();
            for (auto& id : sens_ids) {
                pfr.addSensitivityReaction(id);
            }
            auto sp_names = rctr_parser.getSensitivitySpecies();
            for (auto& sp : sp_names) {
                pfr.addSensitivitySpecies(sp);
            }
            sens_ids.insert(sens_ids.end(),
                            make_move_iterator(sp_names.begin()),
                            make_move_iterator(sp_names.end()));
        }  else { // Full sens enabled. Here all rxnids are counted and species are ignored
            nquad = gas->kinetics()->nReactions();
            for (int i = 0; i < gas->kinetics()->nReactions(); i++){
                sens_ids.push_back(gas->kinetics()->reaction(i)->id);
            }
            for (auto kin : ikin){
                nquad += kin->nReactions();
                for (int i = 0; i < kin->nReactions(); i++){
                    sens_ids.push_back(kin->reaction(i)->id);
                }
            }
        }
    }

    /*
    vector<double> ydot(25);
    vector<double> y(25);
    pfr.getInitialConditions(0, y.data(), ydot.data());
    for (size_t i = 0; i < 25; i++){
        cout << "i:   " << i << "   y: " << y[i] << "   ydot:   " << ydot[i] << endl;
    }
    */

    PFR1dSolver pfr_solver {make_shared<PFR1d>(pfr)};

    //auto simul_node = tube_node["simulation"];
    if (rctr_parser.tolerancesDefined()){
        auto abs_tol = rctr_parser.get_atol();
        auto rel_tol = rctr_parser.get_rtol();
        pfr_solver.setTolerances(rel_tol, abs_tol);
    }

    // Full sensitivity is set through PFR solver
    if (sens_on && full_sens){
        pfr_solver.setQuadratureSize(nquad);
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
    auto rpa_flag = rctr_parser.RPA();
    vector<double> zvals = get_log10_intervals(rctr_len, simul_init_step); 

    vector<double> T_params = rctr_parser.Ts();
    vector<double> P_params = rctr_parser.Ps();
    vector<double> fr_params = rctr_parser.FRs();
    if (!fr_params.size()){
        auto fr = velocity * rctr_xc_area;
        fr_params.push_back(fr);
    }

    auto get_vel = [&](double fr) -> double {
        return fr / rctr_xc_area;
    };

    // Set the output type
    OutputFormat data_format = rctr_parser.printFormat();
    setOutputFormat(data_format);

    auto surf_init_covs = rctr_parser.getSurfPhaseCompositions();
    fs::path curr_dir = ".";
    for (const auto& T : T_params){
        for (const auto& P : P_params){
            for (const auto& fr : fr_params){
                pfr.setVelocity(get_vel(fr));
                string gas_comp = rctr_parser.getGasPhaseComposition();
                gas->thermo()->setState_TPX(T, P, gas_comp);
                for (size_t i = 0; i < surfaces.size(); i++) {
                    surf_ph[i]->setState_TP(T, P);
                    
                    //cout << "Initial Surface Coverages: " << i << endl;
                    //for (auto cov : surf_init_covs[i])
                    //    cout << cov << " ";
                    //cout << endl;
                    //cout << "Density " << surf->siteDensity() << endl;
                    //surf->setCoveragesByName(surf_init_covs[i++]);
                    ikin[i]->solvePseudoSteadyStateProblem();
                }
                pfr.reinit();
                pfr_solver.reinit();

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

                ofstream gas_mole_out((out_dir / ("gas_mole_ss." + file_ext)).string(), ios::out);
                ofstream gas_mass_out((out_dir / ("gas_mass_ss." + file_ext)).string(), ios::out);
                ofstream gas_sdot_out((out_dir / ("gas_sdot_ss." + file_ext)).string(), ios::out);
                ofstream surf_cov_out((out_dir / ("surf_cov_ss." + file_ext)).string(), ios::out);
                ofstream surf_sdot_out((out_dir / ("surf_sdot_ss." + file_ext)).string(), ios::out);
                ofstream state_var_out((out_dir / ("rctr_state_ss." + file_ext)).string(), ios::out);
                ofstream rates_out((out_dir / "rates_ss.out").string(), ios::out);
                print_rxn_rates_hdr(rates_out);

                if (data_format == OutputFormat::DAT) {
                    gas_mole_out << "#Gas Mole fractions\n";
                    gas_mass_out << "#Gas Mass fractions\n";
                    gas_sdot_out << "#Surface Production Rates of  Gas Species (units of kg/s)\n";
                    surf_cov_out << "#Surace Coverages\n";
                    surf_sdot_out << "#Production Rates of Surface Species (units of kmol/s) \n";
                    state_var_out << "#Steady State Reactor State\n";
                }

                print_gas_species_hdr(gas_mole_out, gas->thermo().get(), "z(m)");
                print_gas_species_hdr(gas_mass_out, gas->thermo().get(), "z(m)");
                print_gas_species_hdr(gas_sdot_out, gas->thermo().get(), "z(m)");
                print_surface_species_hdr(surf_cov_out, surfaces, "z(m)");
                print_surface_species_hdr(surf_sdot_out, surfaces, "z(m)");
                print_pfr_state_hdr(state_var_out);

                gas_mole_out.precision(6);
                gas_mass_out.precision(6);
                gas_sdot_out.precision(6);
                surf_cov_out.precision(6);
                surf_sdot_out.precision(6);
                state_var_out.precision(6);
                rates_out.precision(6);

                for (const auto& z : zvals) {
                    pfr_solver.solve(z);
                    print_pfr_rctr_state(z, &pfr, gas_mole_out, gas_mass_out, 
                                         gas_sdot_out, surf_cov_out, surf_sdot_out, state_var_out);
                    if (rpa_flag) {
                        string rpa_file_name = "rates_z-";
                        rpa_file_name += to_string(z);
                        rpa_file_name += ".out";
                        ofstream rates_out ((out_dir / rpa_file_name).string(), ios::out); // Masks the name
                        print_rxn_rates_hdr(//"Rates (mol/s) and Partial Equilibrium Analysis:",
                                            rates_out);
                        rates_out.precision(6);

                        print_rxn_rates(gas->kinetics().get(), rates_out);
                        for (auto surf : surfaces) {
                            print_rxn_rates(surf->kinetics().get(), rates_out);
                        }
                        rates_out.close();
                    }
                }

                pfr_solver.writeStateData((out_dir / "1d_pfr_state.out").string());
                pfr_solver.writeGasData((out_dir / "1d_pfr_gas.out").string());
                pfr_solver.writeSurfaceData((out_dir / "1d_pfr_surface.out").string());
                if (sens_on){
                    string sep = (file_ext == "csv") ? "," : "\t";
                    if (!full_sens)
                        pfr_solver.writeSensitivityData(
                                (out_dir / ("1d_pfr_sensitivity." + file_ext)).string(), sens_ids, sep);
                    else
                        pfr_solver.writeFisherInformationMatrixDiag(
                                (out_dir / ("1d_pfr_sensitivity." + file_ext)).string(), sens_ids, sep);
                }


                // Print final rpa data
                rates_out.precision(6);
                print_rxn_rates(gas->kinetics().get(), rates_out);
                for (auto surf : surfaces) {
                    print_rxn_rates(surf->kinetics().get(), rates_out);
                }

                gas_mole_out.close();
                gas_mass_out.close();
                gas_sdot_out.close();
                surf_cov_out.close();
                surf_sdot_out.close();
                state_var_out.close();
                rates_out.close();
            }
        }
    }
}

}
