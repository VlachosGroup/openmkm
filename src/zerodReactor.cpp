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
#include "reactor.h"


namespace HeteroCt 
{

using namespace std;
using namespace Cantera;

void run_0d_reactor(YAML::Node& tube_node,
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<InterfaceInteractions>> surfaces,
                    ofstream& gen_info, 
                    bool transient_log)
{
    //Define the reactor based on the input file
    auto rctr_node = tube_node["reactor"];
    //auto rctr_type_node = rctr_node["type"];
    //auto rctr_type = rt[rctr_type_node.as<string>()];

    auto rctr = make_shared<Reactor>();
    auto in_rsrv = make_shared<Reservoir>();
    auto exhst = make_shared<Reservoir>();
    rctr->insert(*gas);
    in_rsrv->insert(*gas);
    exhst->insert(*gas);

    // Read the reactor dimensions
    size_t cstr_nos = 1;
    auto nd_node = rctr_node["nodes"];
    if (nd_node)
        if (!nd_node.IsNull())
            cstr_nos = nd_node.as<int>();

    double rctr_vol = strSItoDbl(rctr_node["volume"].as<string>());
    if (!rctr_vol) {
        // TODO: Raise Error
        ;
    }
    rctr_vol /= cstr_nos;

    rctr->setInitialVolume(rctr_vol);
    vector<shared_ptr<ReactorSurface>> cat_surfs;

    double cat_area = strSItoDbl(rctr_node["cat_abyv"].as<string>());
    cat_area *= rctr_vol;
    size_t surf_spec_no = 0;
    for (const auto surf : surfaces) {
        surf_spec_no += surf->nSpecies();
        auto cat_surf = make_shared<ReactorSurface>(); 
        cat_surf->setKinetics(surf.get());
        vector<double> coverages(surf->nSpecies());
        surf->getCoverages(coverages.data());
        cat_surf->setCoverages(coverages.data());
        cat_surf->setArea(cat_area);
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
    outlet->setPressureCoeff(0.00001);
    inlet_mfc->install(*in_rsrv, *rctr);
    outlet->install(*rctr, *exhst);

    cout << "reactor density " << rctr->density() << endl;

    // Start the simulation
    auto simul_node = tube_node["simulation"];
    auto end_time = strSItoDbl(simul_node["end_time"].as<string>());
    auto solver_node = simul_node["solver"];
    auto abs_tol = solver_node["atol"].as<double>();
    auto rel_tol = solver_node["rtol"].as<double>();

    
    ReactorNet rnet; 
    rnet.addReactor(*rctr);
    rnet.setTolerances(rel_tol, abs_tol);
    vector<double> times;

    if (transient_log)
        times = get_log_intervals(end_time);
    else
        times.push_back(end_time);

    vector<double> gas_X(rctr->contents().nSpecies());
    vector<double> surf_cov(surf_spec_no);
    for (size_t i = 0; i < cstr_nos; i++) {
        // The next three lines redefine the TPX of gas for each reactor
        cout << "CSTR #: " << i << endl;
        double temp = rctr->contents().temperature();
        cout << "temp: "  << temp << endl;
        double pressure = rctr->contents().pressure();
        cout << "pressure: "  << pressure << endl;
        rctr->contents().getMoleFractions(gas_X.data());
        cout << "Mole fractions: "  << endl;
        for (auto k = 0; k < gas_X.size(); k++)
            cout << gas_X[k] << " ";
        cout << endl;
        rctr->surface(0)->getCoverages(surf_cov.data());    // Adapt this for multiple surfaces
        cout << "Surface Coverages: "  << endl;
        for (auto k = 0; k < surf_cov.size(); k++)
            cout << surf_cov[k] << " ";
        cout << endl;
        gas->setState_TPX(temp, pressure, gas_X.data());
        in_rsrv->syncState();
        rnet.setInitialTime(0);
        rnet.reinitialize();

        for (const auto & tm : times) {
            rnet.advance(tm);
        }
        rctr->restoreState();

    }
   
    //return rctr;
}

void run_0d_reactor(YAML::Node& tube_node,
                    shared_ptr<IdealGasMix> gas, 
                    vector<shared_ptr<Interface>> surfaces,
                    ofstream& gen_info, 
                    bool transient_log)
{
    //Define the reactor based on the input file
    auto rctr_node = tube_node["reactor"];
    //auto rctr_type_node = rctr_node["type"];
    //auto rctr_type = rt[rctr_type_node.as<string>()];

    auto rctr = make_shared<Reactor>();
    auto in_rsrv = make_shared<Reservoir>();
    auto exhst = make_shared<Reservoir>();
    rctr->insert(*gas);
    in_rsrv->insert(*gas);
    exhst->insert(*gas);

    // Read the reactor dimensions
    size_t cstr_nos = 1;
    auto nd_node = rctr_node["nodes"];
    if (nd_node)
        if (!nd_node.IsNull())
            cstr_nos = nd_node.as<int>();

    double rctr_vol = strSItoDbl(rctr_node["volume"].as<string>());
    if (!rctr_vol) {
        // TODO: Raise Error
        ;
    }
    rctr_vol /= cstr_nos;

    rctr->setInitialVolume(rctr_vol);
    vector<shared_ptr<ReactorSurface>> cat_surfs;

    double cat_area = strSItoDbl(rctr_node["cat_abyv"].as<string>());
    cat_area *= rctr_vol;
    size_t surf_spec_no = 0;
    for (const auto surf : surfaces) {
        surf_spec_no += surf->nSpecies();
        auto cat_surf = make_shared<ReactorSurface>(); 
        cat_surf->setKinetics(surf.get());
        vector<double> coverages(surf->nSpecies());
        surf->getCoverages(coverages.data());
        cat_surf->setCoverages(coverages.data());
        cat_surf->setArea(cat_area);
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
    outlet->setPressureCoeff(0.00001);
    inlet_mfc->install(*in_rsrv, *rctr);
    outlet->install(*rctr, *exhst);

    cout << "reactor density " << rctr->density() << endl;

    // Start the simulation
    auto simul_node = tube_node["simulation"];
    auto end_time = strSItoDbl(simul_node["end_time"].as<string>());
    auto solver_node = simul_node["solver"];
    auto abs_tol = solver_node["atol"].as<double>();
    auto rel_tol = solver_node["rtol"].as<double>();

    
    ReactorNet rnet; 
    rnet.addReactor(*rctr);
    rnet.setTolerances(rel_tol, abs_tol);
    vector<double> times;

    if (transient_log)
        times = get_log_intervals(end_time);
    else
        times.push_back(end_time);

    vector<double> gas_X(rctr->contents().nSpecies());
    vector<double> surf_cov(surf_spec_no);
    for (size_t i = 0; i < cstr_nos; i++) {
        // The next three lines redefine the TPX of gas for each reactor
        double temp = rctr->contents().temperature();
        double pressure = rctr->contents().pressure();
        rctr->contents().getMoleFractions(gas_X.data());
        gas->setState_TPX(temp, pressure, gas_X.data());
        in_rsrv->syncState();
        rnet.setInitialTime(0);
        rnet.reinitialize();

        for (const auto & tm : times) {
            rnet.advance(tm);
        }
        rctr->restoreState();
        cout << "CSTR #: " << i << endl;
        cout << "temp: "  << temp << endl;
        cout << "pressure: "  << pressure << endl;
        cout << "Mole fractions: "  << endl;
        rctr->contents().getMoleFractions(gas_X.data());
        for (auto k = 0; k < gas_X.size(); k++)
            cout << gas_X[k] << " ";
        cout << endl;
        rctr->surface(0)->getCoverages(surf_cov.data());    // Adapt this for multiple surfaces
        cout << "Surface Coverages: "  << endl;
        for (auto k = 0; k < surf_cov.size(); k++)
            cout << surf_cov[k] << " ";
        cout << endl;

    }
   
    //return rctr;
}

}
