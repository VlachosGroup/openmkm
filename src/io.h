#ifndef HTRCT_IO_H
#define HTRCT_IO_H

#include <vector>
#include <string>

#include "cantera/base/config.h"
#include "cantera/kinetics.h"
#include "cantera/thermo/ThermoPhase.h"

void print_formation_enthalpy(std::vector<Cantera::ThermoPhase*> phases, std::string output_file);
void print_formation_entropy(std::vector<Cantera::ThermoPhase*> phases, std::string output_file);
void print_rxn_enthalpy(std::vector<Cantera::Kinetics*> kinetic_mgrs, doublereal T, std::string output_file);
void print_rxn_eq_consts(std::vector<Cantera::Kinetics*> kinetic_mgrs, std::string output_file);

#endif
