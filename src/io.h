/**
 * @file io.h
 */

// This file is part of Hetero_ct.

#ifndef HTRCT_IO_H
#define HTRCT_IO_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>

#include "cantera/base/config.h"
#include "cantera/kinetics.h"
#include "cantera/thermo/ThermoPhase.h"

namespace HeteroCt
{

/**
 * Utility function to print formation enthalpies of species
 */
void print_formation_enthalpy(std::vector<std::shared_ptr<Cantera::ThermoPhase> > phases, 
                              std::string output_file);

/**
 * Utility function to print formation entropies of species
 */
void print_formation_entropy(std::vector<std::shared_ptr<Cantera::ThermoPhase> > phases, 
                             std::string output_file);

/**
 * Utility function to print reaction enthalpies
 */
void print_rxn_enthalpy(std::vector<Cantera::Kinetics*> kinetic_mgrs, 
                        doublereal T, 
                        std::string output_file);

/**
 * Utility function to print reaction entropies
 */
void print_rxn_entropy(std::vector<Cantera::Kinetics*> kinetic_mgrs, 
                       std::string output_file);

/**
 * Utility function to print equilibrium constants of reactions 
 */
void print_rxn_eq_consts(std::vector<Cantera::Kinetics*> kinetic_mgrs, 
                         std::string output_file);

/**
 * Utility function to print Gibbs energy of reactions 
 */
void print_rxn_gibbs(std::vector<Cantera::Kinetics*> kinetic_mgrs, 
                     doublereal T,  
                     std::string output_file);

/**
 * Utility function to print equilibrium constants of reactions 
 */
void print_rxn_kc(std::vector<Cantera::Kinetics*> kinetic_mgrs, 
                  std::string output_file);

/**
 * Utility function to print forward rate constants of reactions 
 */
void print_rxn_kf(std::vector<Cantera::Kinetics*> kinetic_mgrs, 
                  std::string output_file);

/**
 * Utility function to print reverse rate constants of reactions 
 */
void print_rxn_kr(std::vector<Cantera::Kinetics*> kinetic_mgrs, 
                  std::string output_file);

/**
 * Utility function to print hetero_ct header 
 */
void print_omkm_header(std::ostream& out); 


}
#endif
