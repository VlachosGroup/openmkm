#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfLatIntPhase.h"

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
    auto rctr_type_node = getChildNode(m_rctr_nd, "reactor",
                                       vector<string>{"type"});
    return RctrTypeMap[rctr_type_node.as<string>()];
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

