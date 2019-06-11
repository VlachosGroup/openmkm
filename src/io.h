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
#include "cantera/thermo/SurfPhase.h"
#include "cantera/zeroD/Reactor.h"

#include "pfr1d.h"

namespace OpenMKM
{
/**
 * Utility function to print all species in the given phases
 */
void print_species(std::vector<std::shared_ptr<Cantera::ThermoPhase> > phases, 
                   std::string output_file);

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
 * Utility function to print all given reactions 
 */
void print_rxns(std::vector<Cantera::Kinetics*> kinetic_mgrs, 
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
 * Utility function to print OpenMKM header 
 */
void print_omkm_header(std::ostream& out);

/**
 * Utility function to print reaction rates
 */
void print_rxn_rates(Cantera::Kinetics* kin, int rxn_index, std::ofstream& out);

/**
 * Utility function to print header before printing reaction rates
 */
void print_rxn_rates_hdr(std::string hdr, std::ofstream& out);

/**
 * Utility function to print (1d) PFR reactor state at given distance from inlet
 */
void print_pfr_rctr_state(double z, PFR1d* rctr, 
                          std::vector<Cantera::SurfPhase*> surfaces,
                          std::ofstream& gas_mole_out, std::ofstream& gas_mass_out,
                          std::ofstream& gas_msdot_out, std::ofstream& surf_cov_out,
                          std::ofstream& state_var_out);
/**
 * Utility function to print 0d reactor state at given distance from inlet
 */
void print_0d_rctr_state(double z, Cantera::Reactor* rctr, 
                          std::vector<Cantera::SurfPhase*> surfaces,
                          std::ofstream& gas_mole_out, std::ofstream& gas_mass_out,
                          std::ofstream& gas_msdot_out, std::ofstream& surf_cov_out,
                          std::ofstream& state_var_out);

}
#endif
