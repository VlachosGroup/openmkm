#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfLatIntPhase.h"

#include "omkmexceptions.h"
#include "reactor_parser.h"

using namespace std;
using namespace YAML;
using namespace Cantera;

namespace OpenMKM 
{

Node getChildNode(Node& p_nd, string p_name, vector<string> descendants)
{
    string lineage = p_name;
    Node d_nd (p_nd);   // Descending node
    for (const auto& child : descendants)  {
        d_nd = d_nd[child];
        lineage += "." + child;
        if (!d_nd || d_nd.IsNull()){
            throw YAMLParserError("ReactorParser", 
                                  lineage, "node not found or null");
        }
    }
    return d_nd;
}

void ReactorParser::read_mandatory_nodes(){
    m_rctr_nd = m_tube_nd["reactor"];
    m_phase_nd = m_tube_nd["phases"];
    m_simul_nd = m_tube_nd["simulation"];
    m_inlet_nd = m_tube_nd["inlet_gas"];

    if (!m_rctr_nd || m_rctr_nd.IsNull()){
        throw YAMLParserError("ReactorParser::read_mandatory_nodes", 
                              "reactor", "node not found or null");
    }
    // Inlet node could be null for batch reactor
    //
    if (!m_simul_nd || m_simul_nd.IsNull()){
        throw YAMLParserError("ReactorParser::read_mandatory_nodes", 
                              "simulation", "node not found or null");
    }
    if (!m_phase_nd || m_phase_nd.IsNull()){
        throw YAMLParserError("ReactorParser::read_mandatory_nodes", 
                              "phases", "node not found or null");
    }
}

shared_ptr<IdealGasMix> getGasPhase(){
    auto gas_nd = m_phase_nd["gas"];
    if (!gas_nd || gas_nd.IsNull()){
        throw YAMLParserError("ReactorParser::getGasPhase", 
                              "phases.gas", "node not found or null");
    }
    auto gas_name_nd = gas_nd["name"];
    if (!gas_name_nd || gas_name_nd.IsNull()){
        throw YAMLParserError("ReactorParser::getGasPhase", 
                              "phases.gas.name", "node not found or null");
    }
    auto gas_phase_name = phase_node["gas"]["name"].as<string>();
    return make_shared<IdealGasMix>(phase_file_name, gas_phase_name);
}



/*
bool validate() { // TODO: Implement for one shot error checking
    return false;
}
*/

}

#endif
