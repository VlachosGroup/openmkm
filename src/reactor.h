/**
 *  @file ZerodReactor.h
 *  Header for all functions that read the reactor parameters 
 *  from YAML input file and run the given ZeroD Reactor type
 *  with the supplied parameters
 */

// This file is part of hetero_ct. See License.txt in the top level directory

#ifndef HTRCT_ZEROD_REACTOR_H
#define HTRCT_ZEROD_REACTOR_H

#include <vector>
#include <memory>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include "cantera/IdealGasMix.h"
#include "cantera/InterfaceLatInt.h"

namespace HeteroCt
{

enum RctrType {
    BATCH,
    CSTR,
    PFR_0D,
    PFR
};


void run_0d_reactor(RctrType rctr_type, 
                    YAML::Node& tube_node,
                    std::shared_ptr<Cantera::IdealGasMix> gas, 
                    std::vector<std::shared_ptr<Cantera::InterfaceInteractions>> surfaces,
                    std::ofstream& gen_info);

/*
void run_1d_reactor(YAML::Node& tube_node,
                    std::shared_ptr<Cantera::IdealGasMix> gas, 
                    std::vector<std::shared_ptr<Cantera::Interface>> surfaces,
                    std::ofstream& gen_info);
                    */

void run_1d_reactor(YAML::Node& tube_node,
                    std::shared_ptr<Cantera::IdealGasMix> gas, 
                    std::vector<std::shared_ptr<Cantera::InterfaceInteractions>> surfaces,
                    std::ofstream& gen_info);
                    

}

#endif 
