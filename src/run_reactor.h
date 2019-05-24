/**
 *  @file run_reactor.h
 *  Header for all functions that read the reactor parameters 
 *  from YAML input file and run the given Reactor type
 *  with the supplied parameters
 */

// This file is part of hetero_ct. See License.txt in the top level directory

#ifndef OMKM_REACTOR_RUN_H
#define OMKM_REACTOR_RUN_H

#include <vector>
#include <memory>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include "cantera/IdealGasMix.h"
#include "cantera/InterfaceLatInt.h"

#include "reactor_parser.h"

namespace OpenMKM
{

void run_0d_reactor(ReactorParser& rctr_parser,
                    std::shared_ptr<Cantera::IdealGasMix> gas, 
                    std::vector<std::shared_ptr<Cantera::InterfaceInteractions>> surfaces,
                    std::ofstream& gen_info);

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
