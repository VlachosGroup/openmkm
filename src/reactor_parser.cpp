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

void operator >> (const Node& node, vector_fp& vec)
{
    for (size_t i = 0; i < node.size(); i++){
        vec.push_back(strSItoDbl(node[i].as<string>()));
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

bool ReactorParser::gasPhaseDefined(string phase_filename)
{
    return IsChildNodeAvailable(m_phase_nd, vector<string>{"name", "gas"});
}

shared_ptr<IdealGasMix> ReactorParser::getGasPhase(string phase_filename)
{
    auto gas_name_nd = getChildNode(m_phase_nd, "tube.phases", 
                                    vector<string>{"name", "gas"});
    auto gas_phase_name = gas_name_nd.as<string>();
    auto gas = make_shared<IdealGasMix>(phase_filename, gas_phase_name);

    m_gas_X = getGasPhaseComposition();
    gas->setState_TPX(m_T, m_P, m_gas_X);
    return gas;
}

string ReactorParser::getGasPhaseComposition()
{
    if (!m_gas_X.size()){
        auto gas_X_nd = getChildNode(m_phase_nd, "tube.phases", 
                                     vector<string>{"initial_state", "gas"});
        m_gas_X = gas_X_nd.as<string>();
    }
    return m_gas_X;
}

bool ReactorParser::bulkPhaseDefined(string phase_filename)
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

bool ReactorParser::surfacePhasesDefined(string phase_filename)
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
        m_surf_X.push_back(surf_init_state);
        surf_phases.push_back(surf);
        surf_phases1.push_back(surf.get());
    }
    setTotalSiteDensity(surf_phases1);

    return surf_phases;
}

vector<string> ReactorParser::getSurfPhaseCompositions()
{
    return m_surf_X; 
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

//! Output Format
OutputFormat ReactorParser::printFormat()
{
    /*   The following code may result in error if enums have a default value of 0
    if (m_output_format == OutputFormat::DAT || m_output_format == OutputFormat::CSV)
        return m_output_format;
    */

    auto kids = vector<string>{"output_format"};
    if(!IsChildNodeAvailable(m_simul_nd, kids)){
        m_output_format = OutputFormat::DAT;          // Default is DAT
        return m_output_format;
    }
    auto print_format_nd = getChildNode(m_simul_nd, "simulation", kids);
    auto print_format = print_format_nd.as<string>();
    switch(toupper(print_format[0])){
        case 'C':
            m_output_format = OutputFormat::CSV;
            break;
        case 'D':
            m_output_format = OutputFormat::DAT;
            break;
        default:
            cout << "WARNING: Invalid output format given: Valid options are CSV or DAT. Choosing default DAT" << endl;
            m_output_format = OutputFormat::DAT;
    }

    return m_output_format;
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
    return strSItoDbl(initstep_nd.as<string>());
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

// TPD Params
// TPD end temperature
double ReactorParser::getTPDEndTemp()
{
    auto tend_nd = getChildNode(m_rctr_nd, "tube.reactor",
                                vector<string>{"Tend"});
    return strSItoDbl(tend_nd.as<string>());
}

// TPD temperature ramp
double ReactorParser::getTPDTempRamp()
{
    auto tramp_nd = getChildNode(m_rctr_nd, "tube.reactor",
                                 vector<string>{"Tramp"});
    return strSItoDbl(tramp_nd.as<string>());
}

// Simulation end time for CSTR and batch reactors 
double ReactorParser::getEndTime()
{
    auto end_time_nd = getChildNode(m_simul_nd, "tube.simulation",
                                    vector<string>{"end_time"});
    return strSItoDbl(end_time_nd.as<string>());
}


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
    return strSItoDbl(initstep_nd.as<string>());
}

//! Parses the flag to enable transient output 
bool ReactorParser::logTransient()
{
    auto transient_flag = false;
    auto kids = vector<string>{"transient"};
    if(IsChildNodeAvailable(m_simul_nd, kids)){
        auto transient_nd = getChildNode(m_simul_nd, "tube.simulation", kids);
        transient_flag = transient_nd.as<bool>();
    }
    if (!transient_flag) {// Enable if tpd is enabled
        if(getMode() == "tpd"){
            transient_flag = true;
            cout << "Simulation in TPD Mode: " << endl
                 << "Enabling printing varaibles at transient state" << endl;
        }
    }
    return transient_flag;
}

//! Parses the flag to enable transient output 
string ReactorParser::steppingType()
{
    if(getMode() == "tpd"){
        cout << "Simulation in TPD Mode: Regular stepping." << endl;
        return "regular";
    }
    auto stepping_nd = getChildNode(m_simul_nd, "tube.simulation", 
                                    vector<string>{"stepping"});
    return stepping_nd.as<string>();
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
bool ReactorParser::isT_multi_input() 
{
    vector<string> kids {"temperature", "multi_input"};
    return IsChildNodeAvailable(m_simul_nd, kids);
}

bool ReactorParser::isP_multi_input() 
{
    vector<string> kids {"pressure", "multi_input"};
    return IsChildNodeAvailable(m_simul_nd, kids);
}

bool ReactorParser::isFR_multi_input() 
{
    vector<string> massflowrate_kids {"flow_rate", "multi_input"};
    auto check_mfr = IsChildNodeAvailable(m_simul_nd, massflowrate_kids);
    return check_mfr;
}

/*
bool validate() { // TODO: Implement for one shot error checking
    return false;
}
*/

std::vector<double> ReactorParser::Ts()
{
    vector_fp Ts;
    if (isT_multi_input()){
        vector<string> kids {"temperature", "multi_input"};
        auto parameter_T_nd = getChildNode(m_simul_nd, "simulation", kids);
        parameter_T_nd >> Ts;
    }
    else {
        Ts.push_back(m_T);
    }
    return Ts;
}

std::vector<double> ReactorParser::Ps()
{
    vector_fp Ps;
    if (isP_multi_input()){
        vector<string> kids {"pressure", "multi_input"};
        auto parameter_P_nd = getChildNode(m_simul_nd, "simulation", kids);
        parameter_P_nd >> Ps;
    }
    else {
        Ps.push_back(m_P);
    }
    return Ps;
}

std::vector<double> ReactorParser::FRs()
{
    vector_fp fr;
    if(isFR_multi_input()){
        vector<string> flowrate_kids {"flow_rate", "multi_input"};
        auto parameter_fr_nd = getChildNode(m_simul_nd, "simulation", flowrate_kids);
        parameter_fr_nd >> fr;
    }
    return fr;
}

// Parsing sensitivity related variables
bool ReactorParser::isSensitivityAnalysisEnabled()
{
    vector<string> kids {"sensitivity"};
    return IsChildNodeAvailable(m_simul_nd, kids);
}

bool ReactorParser::isfullSensitivityAnalysis()
{
    vector<string> kids {"full", "sensitivity"};
    if(!IsChildNodeAvailable(m_simul_nd, kids)){
        return false;
    }
    auto fullSA_nd = getChildNode(m_simul_nd, "tube.sumulation", kids);
    return fullSA_nd.as<bool>();
}
double ReactorParser::getSensitivityRtol()
{

    vector<string> kids {"rtol", "sensitivity"};
    if (!IsChildNodeAvailable(m_simul_nd, kids)){
        return 0.0;
    }
    auto rtol_nd = getChildNode(m_simul_nd, "tube.simulation", kids);
    return rtol_nd.as<double>();
}

double ReactorParser::getSensitivityAtol()
{
    vector<string> kids {"atol", "sensitivity"};
    if (!IsChildNodeAvailable(m_simul_nd, kids)){
        return 0.0;
    }
    auto atol_nd = getChildNode(m_simul_nd, "tube.simulation", kids);
    return atol_nd.as<double>();
}

vector<string> ReactorParser::getSensitivityReactions()
{
    vector<string> rxn_ids;
    vector<string> kids {"reactions", "sensitivity"};
    auto SArxn_nd = getChildNode(m_simul_nd, "simulation", kids);
    for(YAML::const_iterator it = SArxn_nd.begin(); it != SArxn_nd.end(); ++it) {
        rxn_ids.push_back(it->as<string>());
    }
    return rxn_ids;
}

vector<string> ReactorParser::getSensitivitySpecies()
{
    vector<string> sp_names;
    vector<string> kids {"species", "sensitivity"};
    auto SArxn_nd = getChildNode(m_simul_nd, "simulation", kids);
    for(YAML::const_iterator it = SArxn_nd.begin(); it != SArxn_nd.end(); ++it) {
        sp_names.push_back(it->as<string>());
    }
    return sp_names;
}


}
