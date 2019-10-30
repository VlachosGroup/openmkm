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
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<InterfaceInteractions>> surfaces,
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
        auto flow_rate = mfr / gas->density();
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
    for (const auto surf: surfaces) {
        ikin.push_back(surf.get());
        surf_ph.push_back(surf.get());
    }

    // Before simulation save the initial coverages of surface species for reuse
    //vector<vector<double>> surf_init_covs;
    vector<double> cov;
    for (const auto surf: surfaces) {
        cov.resize(surf->nSpecies());
        surf->getCoverages(cov.data());
        //surf_init_covs.push_back(cov);
    }

    // Start the simulation
    gen_info << "Solving for equilibirum surface coverages at PFR inlet" << endl;
    for (const auto surf: surfaces) {
        cout << "Surface Site Density " << surf->density() << endl;
        surf->solvePseudoSteadyStateProblem();
        vector<double> cov(surf->nSpecies());
        surf->getCoverages(cov.data());

        gen_info << "Equilibrium surface coverages on Surface: " 
                 <<  surf->name() << endl;
        for (auto i = 0; i < surf->nSpecies(); i++)
            gen_info << surf->speciesSPName(i) << " coverage: " << cov[i] << endl;
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
    vector<double> zvals = get_log10_intervals(rctr_len, simul_init_step); //Use the same function to get z steps

    
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
                gas->setState_TPX(T, P, gas_comp);
                size_t i = 0;
                for (const auto surf : surfaces) {
                    surf->setState_TP(T, P);
                    
                    cout << "Initial Surface Coverages: " << i << endl;
                    for (auto cov : surf_init_covs[i])
                        cout << cov << " ";
                    cout << endl;
                    cout << "Density " << surf->density() << endl;
                    //surf->setCoveragesByName(surf_init_covs[i++]);
                    surf->solvePseudoSteadyStateProblem();
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
                ofstream gas_msdot_out((out_dir / ("gas_msdot_ss." + file_ext)).string(), ios::out);
                ofstream surf_cov_out((out_dir / ("surf_cov_ss." + file_ext)).string(), ios::out);
                ofstream state_var_out((out_dir / ("rctr_state_ss." + file_ext)).string(), ios::out);
                ofstream rates_out((out_dir / "rates_ss.out").string(), ios::out);

                gas_mole_out.precision(6);
                gas_mass_out.precision(6);
                gas_msdot_out.precision(6);
                surf_cov_out.precision(6);
                state_var_out.precision(6);
                rates_out.precision(6);

                for (const auto& z : zvals) {
                    pfr_solver.solve(z);
                    print_pfr_rctr_state(z, &pfr, surf_ph, gas_mole_out, gas_mass_out, 
                                         gas_msdot_out, surf_cov_out, state_var_out);
                    if (rpa_flag) {
                        string rpa_file_name = "rates_z-";
                        rpa_file_name += to_string(z);
                        rpa_file_name += ".out";
                        ofstream rates_out ((out_dir / rpa_file_name).string(), ios::out); // Masks the name
                        print_rxn_rates_hdr(//"Rates (mol/s) and Partial Equilibrium Analysis:",
                                            rates_out);
                        rates_out.precision(6);

                        print_rxn_rates(gas.get(), rates_out);
                        for (auto surf : surfaces) {
                            print_rxn_rates(surf.get(), rates_out);
                        }
                        rates_out.close();
                    }
                }

                pfr_solver.writeStateData((out_dir / "1d_pfr_state.out").string());
                pfr_solver.writeGasData((out_dir / "1d_pfr_gas.out").string());
                pfr_solver.writeSurfaceData((out_dir / "1d_pfr_surface.out").string());

                // Print final rpa data
                rates_out.precision(6);
                print_rxn_rates(gas.get(), rates_out);
                for (auto surf : surfaces) {
                    print_rxn_rates(surf.get(), rates_out);
                }

                gas_mole_out.close();
                gas_mass_out.close();
                gas_msdot_out.close();
                surf_cov_out.close();
                state_var_out.close();
                rates_out.close();
            }
        }
    }
}

}
