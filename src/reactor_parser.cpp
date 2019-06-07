#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfLatIntPhase.h"
#include "cantera/base/stringUtils.h"

#include "omkmexceptions.h"
#include "reactor_parser.h"

using namespace std;
using namespace YAML;
using namespace Cantera;

namespace OpenMKM 
{

std::map<std::string, RctrType> RctrTypeMap = {{"batch", BATCH},
                                               {"cstr", CSTR},
                                               {"pfr_0d", PFR_0D},
                                               {"pfr", PFR}};

// The descendents are specified in reverse order
Node getChildNode(Node& p_nd, string p_name, vector<string> rev_descendants)
{
    if (rev_descendants.size()){
        auto dir_descndnt = rev_descendants.back();
        rev_descendants.pop_back();
        auto d_nd  = p_nd[dir_descndnt];   // Descending node
        auto lineage = p_name + "." + dir_descndnt;
        if (!d_nd || d_nd.IsNull()){
            throw YAMLParserError("ReactorParser", 
                                  lineage, "node not found or null");
        } else {
            return getChildNode(d_nd, lineage, rev_descendants); 
        }
    } else {
        return p_nd;
    }
}

bool IsChildNodeAvailable(Node& p_nd, vector<string> rev_descendants)
{
    if (rev_descendants.size()){
        auto d_nd  = p_nd[rev_descendants.back()];   // Descending node
        rev_descendants.pop_back();
        if (!d_nd || d_nd.IsNull()){
            return false; 
        } else {
            return IsChildNodeAvailable(d_nd, rev_descendants); 
        }
    } else {
        return true;
    }
}

void ReactorParser::read_mandatory_nodes()
{
    m_rctr_nd = getChildNode(m_tube_nd, "tube", 
                             vector<string>{"reactor"});
    m_phase_nd = getChildNode(m_tube_nd, "tube", 
                              vector<string>{"phases"});
    m_simul_nd = getChildNode(m_tube_nd, "tube", 
                              vector<string>{"simulation"});
    // Inlet node could be null for batch reactor
    if (getReactorType() != BATCH){
        m_inlet_nd = getChildNode(m_tube_nd, "tube", 
                                  vector<string>{"inlet_gas"});
    }

    auto T_nd = getChildNode(m_tube_nd, "tube", 
                             vector<string>{"temperature", "reactor"});
    m_T = strSItoDbl(T_nd.as<string>());

    auto P_nd = getChildNode(m_tube_nd, "tube", 
                             vector<string>{"pressure", "reactor"});
    m_P = strSItoDbl(P_nd.as<string>());
}

bool ReactorParser::GasPhaseDefined(string phase_filename)
{
    return IsChildNodeAvailable(m_phase_nd, vector<string>{"name", "gas"});
}


shared_ptr<IdealGasMix> ReactorParser::getGasPhase(string phase_filename)
{
    auto gas_name_nd = getChildNode(m_phase_nd, "tube.phases", 
                                    vector<string>{"name", "gas"});
    auto gas_phase_name = gas_name_nd.as<string>();
    auto gas = make_shared<IdealGasMix>(phase_filename, gas_phase_name);

    auto gas_X_nd = getChildNode(m_phase_nd, "tube.phases", 
                                 vector<string>{"initial_state", "gas"});
    auto gas_X = gas_X_nd.as<string>();
    gas->setState_TPX(m_T, m_P, gas_X);
    return gas;
}

bool ReactorParser::BulkPhaseDefined(string phase_filename)
{
    return IsChildNodeAvailable(m_phase_nd, vector<string>{"name", "bulk"});
}

shared_ptr<StoichSubstance> ReactorParser::getBulkPhase(string phase_filename)
{
    auto blk_name_nd = getChildNode(m_phase_nd, "tube.phases", 
                                    vector<string>{"name", "bulk"});
    auto blk_phase_name = blk_name_nd.as<string>();
    auto blk = make_shared<StoichSubstance>(phase_filename, blk_phase_name);
    blk->setState_TP(m_T, m_P);
    return blk;
}

bool ReactorParser::SurfacePhasesDefined(string phase_filename)
{
    return IsChildNodeAvailable(m_phase_nd, vector<string>{"surfaces"});
}

vector<shared_ptr<InterfaceInteractions>> ReactorParser::getSurfPhases(
        string phase_filename, vector<ThermoPhase*> gb_phases)
{
    auto surf_nds = m_phase_nd["surfaces"];
    if (!surf_nds || surf_nds.IsNull()){
        throw YAMLParserError("ReactorParser::getSurfPhases", 
                              "phases.surfaces", "node not found or null");
    }

    vector<shared_ptr<InterfaceInteractions>> surf_phases;
    vector<SurfPhase*> surf_phases1;

    for (const auto & surf_nd :  surf_nds){
        auto surf_name = surf_nd["name"].as<string>();
        auto surf = make_shared<InterfaceInteractions>(phase_filename, 
                                                       surf_name, gb_phases);
        auto surf_init_state = surf_nd["initial_state"].as<string>();
        surf->setState_TP(m_T, m_P);
        surf->setCoveragesByName(surf_init_state);
        surf_phases.push_back(surf);
        surf_phases1.push_back(surf.get());
    }
    setTotalSiteDensity(surf_phases1);

    return surf_phases;
}

// Get Reactor Type
RctrType ReactorParser::getReactorType()
{
    auto rctr_type_nd = getChildNode(m_rctr_nd, "tube.reactor",
                                     vector<string>{"type"});
    return RctrTypeMap[rctr_type_nd.as<string>()];
}

//! Return reactor (PFR) cross section area
double ReactorParser::getXCArea()
{
    auto area_nd = getChildNode(m_rctr_nd, "tube.reactor",
                                vector<string>{"area"});
    return strSItoDbl(area_nd.as<string>());
}

//! Return reactor (PFR) length 
double ReactorParser::getLength()
{
    auto length_nd = getChildNode(m_rctr_nd, "tube.reactor",
                                  vector<string>{"length"});
    return strSItoDbl(length_nd.as<string>());
}

// Returns reactor volume 
double ReactorParser::getVolume()
{
    auto vol_nd = getChildNode(m_rctr_nd, "tube.reactor",
                               vector<string>{"volume"});
    return strSItoDbl(vol_nd.as<string>());
}

// Returns # of CSTRs to use for PFR_0D
size_t ReactorParser::getNodes()
{
    size_t rctr_nodes = 1;
    auto kids = vector<string>{"nodes"};
    if(!IsChildNodeAvailable(m_rctr_nd, kids)){
        return rctr_nodes;
    }
    auto rctr_node_nd = getChildNode(m_rctr_nd, "tube.reactor", kids);
    return rctr_node_nd.as<size_t>();
}


bool ReactorParser::FlowRateDefined()
{
    return IsChildNodeAvailable(m_inlet_nd, vector<string>{"flow_rate"});
}

bool ReactorParser::MassFlowRateDefined()
{
    return IsChildNodeAvailable(m_inlet_nd, vector<string>{"mass_flow_rate"});
}

bool ReactorParser::ResidenceTimeDefined()
{
    return IsChildNodeAvailable(m_inlet_nd, vector<string>{"residence_time"});
}

double ReactorParser::getFlowRate()
{
    auto fr_nd = getChildNode(m_inlet_nd, "tube.inlet_gas",
                              vector<string>{"flow_rate"});
    return strSItoDbl(fr_nd.as<string>());
}

double ReactorParser::getMassFlowRate()
{
    auto mfr_nd = getChildNode(m_inlet_nd, "tube.inlet_gas",
                               vector<string>{"mass_flow_rate"});
    return strSItoDbl(mfr_nd.as<string>());
}

double ReactorParser::getResidenceTime()
{
    auto rt_nd = getChildNode(m_inlet_nd, "tube.inlet_gas",
                              vector<string>{"residence_time"});
    return strSItoDbl(rt_nd.as<string>());
}

bool ReactorParser::catalystAreaDefined()
{
    return IsChildNodeAvailable(m_rctr_nd, vector<string>{"cat_abyv"});
}

//! Catalyst Area by Reactor Volume
double ReactorParser::getCatalystAbyV()
{
    auto kids = vector<string>{"cat_abyv"};
    if(!IsChildNodeAvailable(m_rctr_nd, kids)){
        return 0.0;
    }
    auto cat_abyv_nd = getChildNode(m_rctr_nd, "tube.reactor", kids);
    return strSItoDbl(cat_abyv_nd.as<string>());
}

//! Get Reactor Operational Modes
//! The implemented modes are 
//! "isothermal" -- isothermal operation, 
//! "Tprofile"   -- Temperature profile along PFR 
//! "TPD"        -- Temperature increased as a function of time 
//! "adiabatic"  -- Adiabatic operation,
//! "heat"       -- Heat conducting Walls
string ReactorParser::getMode()
{
    auto mode_nd = getChildNode(m_rctr_nd, "tube.reactor",
                                vector<string>{"mode"});
    return mode_nd.as<string>();
}

//! Get Temperature Profile imposed on PFR
//! The calling code should call this only if operational mode is "Tprofile".
map<double, double> ReactorParser::getTProfile()
{
    auto tprofile_nd = getChildNode(m_rctr_nd, "tube.reactor",
                                    vector<string>{"TProfile"});
    /*if (!tprofile_nd.IsSequence()){
        throw YAMLParserError("ReactorParser", "tube.reactor.TProfile",
                              "Provide a sequence of dist: T");
    }
    */
    map<double, double> T_profile;
    for(YAML::const_iterator it = tprofile_nd.begin(); it != tprofile_nd.end(); ++it) {
        T_profile.insert(pair<double, double>(strSItoDbl(it->first.as<string>()),
                                              strSItoDbl(it->second.as<string>())));
    }
    return T_profile;
}

//! Get heat transfer coefficient of wall 
//! The calling code should call this only if operational mode is "heat".
double ReactorParser::getWallHeatTransferCoeff()
{
    auto htc_nd = getChildNode(m_rctr_nd, "tube.reactor",
                               vector<string>{"htc"});
    return strSItoDbl(htc_nd.as<string>());
}

//! Get wall specific area (wall area by reactor volume)
//! The calling code should call this only if operational mode is "heat".
double ReactorParser::getWallSpecificArea()
{
    auto wall_nd = getChildNode(m_rctr_nd, "tube.reactor",
                               vector<string>{"wall_abyv"});
    return strSItoDbl(wall_nd.as<string>());
}

//! Get temperature of heat source attached to wall 
//! The calling code should call this only if operational mode is "heat".
double ReactorParser::getExternalTemp()
{
    auto text_nd = getChildNode(m_rctr_nd, "tube.reactor",
                                vector<string>{"Text"});
    return strSItoDbl(text_nd.as<string>());
}


// Parameters to pass to Numerical Solver

//! Checks if atol and rtol are defined
bool ReactorParser::tolerancesDefined()
{
    vector<string> atol_kids {"atol", "solver"};
    bool check_atol = IsChildNodeAvailable(m_simul_nd, atol_kids);
    vector<string> rtol_kids {"rtol", "solver"};
    bool check_rtol = IsChildNodeAvailable(m_simul_nd, rtol_kids);
    return check_atol && check_rtol;
}

//! Parses atol to numerical solver if defined
double ReactorParser::get_atol()
{
    auto atol_nd = getChildNode(m_simul_nd, "tube.simulation",
                                vector<string>{"atol", "solver"});
    return atol_nd.as<double>();
}

//! Parses rtol to numerical solver if defined
double ReactorParser::get_rtol()
{
    auto rtol_nd = getChildNode(m_simul_nd, "tube.simulation",
                                vector<string>{"rtol", "solver"});
    return rtol_nd.as<double>();
}

//! Checks if initial step size to solver is defined
bool ReactorParser::solverInitStepSizeDefined()
{
    return IsChildNodeAvailable(m_simul_nd, 
                                vector<string> {"init_step_size", "solver"});
}

//! Parses initial step size to numerical solver if defined
double ReactorParser::getSolverInitStepSize()
{
    auto initstep_nd = getChildNode(m_simul_nd, "tube.simulation",
                                    vector<string>{"init_step_size", "solver"});
    return initstep_nd.as<double>();
}

//! Checks if maximum no of steps taken by solver is defined
bool ReactorParser::solverMaxStepsDefined()
{
    return IsChildNodeAvailable(m_simul_nd, 
                                vector<string> {"max_steps", "solver"});
}

//! Parses maximum steps taken by numerical solver if defined
double ReactorParser::getSolverMaxSteps()
{
    auto maxsteps_nd = getChildNode(m_simul_nd, "tube.simulation",
                                    vector<string>{"max_steps", "solver"});
    return maxsteps_nd.as<double>();
}

// Simulation parameters that are not part of solver

//! Checks if initial step where simulation ouput is printed is defined
bool ReactorParser::initStepDefined()
{
    return IsChildNodeAvailable(m_simul_nd, 
                                vector<string> {"init_step"});
}

//! Parses initial step where simulation output is printed if defined
double ReactorParser::getInitStep()
{
    auto initstep_nd = getChildNode(m_simul_nd, "tube.simulation",
                                    vector<string>{"init_step"});
    return initstep_nd.as<double>();
}




// IO Flags at various points
//! RPA at intermediate times or distances
bool ReactorParser::RPA()
{
    auto kids = vector<string>{"rpa"};
    if(!IsChildNodeAvailable(m_simul_nd, kids)){
        return false;
    }
    auto rpa_nd = getChildNode(m_simul_nd, "tube.simulation", kids);
    return rpa_nd.as<bool>();
}

// Parametric study related
bool ReactorParser::T_parametric_study_enabled() 
{
    vector<string> kids {"temperature", "parametric_study"};
    return IsChildNodeAvailable(m_simul_nd, kids);
}

bool ReactorParser::P_parametric_study_enabled() 
{
    vector<string> kids {"pressure", "parametric_study"};
    return IsChildNodeAvailable(m_simul_nd, kids);
}

bool ReactorParser::mdot_parametric_study_enabled() 
{
    vector<string> flowrate_kids {"flow_rate", "parametric_study"};
    auto check_fr = IsChildNodeAvailable(m_simul_nd, flowrate_kids);
    vector<string> massflowrate_kids {"mass_flow_rate", "parametric_study"};
    auto check_mfr = IsChildNodeAvailable(m_simul_nd, massflowrate_kids);
    vector<string> residtime_kids {"residence_time", "parametric_study"};
    auto check_rt = IsChildNodeAvailable(m_simul_nd, residtime_kids);
    return (check_fr || check_mfr || check_rt);
}
/*
bool validate() { // TODO: Implement for one shot error checking
    return false;
}
*/

}

