#ifndef OMKM_RCTR_PARSER_H
#define OMKM_RCTR_PARSER_H

#include <vector>
#include <string>
#include <memory>
#include <yaml-cpp/yaml.h>

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfLatIntPhase.h"

#include "omkmexceptions.h"


namespace OpenMKM 
{

YAML::Node getChildNode(YAML::Node& p_nd, 
                        string p_name, 
                        std::vector<string> descendants);

class ReactorParser {
public:
    ReactorParser() {}
    ReactorParser(std::string rctr_file) {
        LoadFile(rctr_file);
    }

    void LoadFile(std::string rctr_file) {
        m_tube_nd = YAML::LoadFile(rctr_file);
        read_mandatory_nodes();
    }

    void read_mandatory_nodes();
    
    std::shared_ptr<Cantera::IdealGasMix> getGasPhase();

    bool validate() { // TODO: Implement for one shot error checking
        return false;
    }

private:
    YAML::Node m_tube_nd;
    // The below YAML nodes are child nodes of tube_nd, but 
    // defined explicitly for convenience
    YAML::Node m_rctr_nd;       
    YAML::Node m_simul_nd;
    YAML::Node m_phase_nd;
    YAML::Node m_inlet_nd;
};

}

#endif
