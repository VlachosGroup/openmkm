/**
 *  @file ZerodReactor.h
 *  Header for all functions that read the reactor parameters 
 *  from YAML input file and run the given ZeroD Reactor type
 *  with the supplied parameters
 */

// This file is part of hetero_ct. See License.txt in the top level directory

#ifndef OMKM_ZEROD_REACTOR_H
#define OMKM_ZEROD_REACTOR_H

#include <vector>
#include <memory>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include "cantera/IdealGasMix.h"
#include "cantera/InterfaceLatInt.h"

#include "reactor_parser.h"

namespace OpenMKM
{

/*
enum RctrType {
    BATCH,
    CSTR,
    PFR_0D,
    PFR
};
*/


void run_0d_reactor(RctrType rctr_type, 
                    YAML::Node& tube_node,
                    ReactorParser& rctr_parser,
                    std::shared_ptr<Cantera::IdealGasMix> gas, 
                    std::vector<std::shared_ptr<Cantera::InterfaceInteractions>> surfaces,
                    std::ofstream& gen_info);

void run_0d_reactor(RctrType rctr_type, 
                    YAML::Node& tube_node,
                    ReactorParser& rctr_parser,
                    std::shared_ptr<Cantera::IdealGasMix> gas, 
                    std::vector<std::shared_ptr<Cantera::Interface>> surfaces,
                    std::ofstream& gen_info);


//void run_1d_reactor(YAML::Node& tube_node, ReactorParser& rctr_parser,
void run_1d_reactor(ReactorParser& rctr_parser,
                    std::shared_ptr<Cantera::IdealGasMix> gas, 
                    std::vector<std::shared_ptr<Cantera::Interface>> surfaces,
                    std::ofstream& gen_info);

void run_1d_reactor(ReactorParser& rctr_parser,
                    std::shared_ptr<Cantera::IdealGasMix> gas, 
                    std::vector<std::shared_ptr<Cantera::InterfaceInteractions>> surfaces,
                    std::ofstream& gen_info);
                    

}

#endif 
