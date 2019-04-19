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
                    vector<shared_ptr<InterfaceInteractions>> surfaces)
{
    //Define the reactor based on the input file
    auto rctr_node = tube_node["reactor"];
    //auto rctr_type_node = rctr_node["type"];
    //auto rctr_type = rt[rctr_type_node.as<string>()];

    auto rctr = make_shared<Reactor>();
    shared_ptr<Reservoir> in_rsrv (new Reservoir), exhst (new Reservoir);
    rctr->insert(*gas);
    in_rsrv->insert(*gas);
    exhst->insert(*gas);

    // Read the reactor dimensions
    int nodes = 1;
    auto nd_node = rctr_node["nodes"];
    if (nd_node)
        if (!nd_node.IsNull())
            nodes = nd_node.as<int>();

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
    for (const auto surf : surfaces) {
        cout << "surf address: " << surf << endl;
        shared_ptr<ReactorSurface> cat_surf(new ReactorSurface());
        cat_surf->setKinetics(surf.get());
        vector<double> coverages(surf->nSpecies());
        surf->getCoverages(coverages.data());
        cat_surf->setCoverages(coverages.data());
        cat_surf->setArea(cat_area);
        auto* ip = cat_surf.get();
        cout << "cat_surf.get() address: " << ip << endl;
        cat_surfs.push_back(cat_surf);
        rctr->addSurface(cat_surf.get());
    }
    shared_ptr<ReactorSurface> cat_surf1(new ReactorSurface());
    auto* ip = cat_surf1.get();
    cout << "cat_surf1.get() raw pointer address: " << ip << endl;
    shared_ptr<ReactorSurface> cat_surf2(new ReactorSurface());
    ip = cat_surf2.get();
    cout << "cat_surf2.get() raw pointer address: " << ip << endl;

    auto *s1 = rctr->surface(0);
    cout << "s1 thermo address: " << s1->thermo() << endl;
    auto *s2 = rctr->surface(1);
    cout << "s2 thermo address: " << s2->thermo() << endl;

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
    outlet->setPressureCoeff(0.0001);
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
    auto times = get_times(end_time);
    for (const auto & tm : times) {
        //if (tm < 2e-6) {
            rnet.advance(tm);//(tm);
        //}
    }
   

    //return rctr;
}

}
