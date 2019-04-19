#include <iostream>
#include <fstream>

#include "io.h"

using namespace std;
using namespace Cantera;

namespace HeteroCt
{

void print_formation_enthalpy(vector<ThermoPhase*> phases, string output_file) 
{
    vector<doublereal> hform; 
    ofstream out;
    out.open(output_file);
    out << "#Dimensionless formation enthalpies of species (H/RT)\n" << endl;
    for  (const auto phase: phases){
        hform.resize(phase->nSpecies());
        phase->getEnthalpy_RT(hform.data());
        for (size_t k = 0; k < phase->nSpecies(); k++) {
            out.width(12); 
            out << std::left << phase->speciesName(k); 
            out.width(12); 
            out << std::right << hform[k] << endl;
        }   
    }   
    out.close();
}

void print_formation_entropy(vector<ThermoPhase*> phases, string output_file) 
{
    vector<doublereal> sform; 
    ofstream out;
    out.open(output_file);
    out << "#Dimensionless formation entropies of species (S/R)\n" << endl;
    for  (const auto phase: phases){
        sform.resize(phase->nSpecies());
        phase->getEntropy_R(sform.data());
        for (size_t k = 0; k < phase->nSpecies(); k++) {
            out.width(12); 
            out << std::left << phase->speciesName(k); 
            out.width(12); 
            out << std::right << sform[k] << endl;
        }   
    }   
    out.close();
}


void print_rxn_enthalpy(vector<Kinetics*> kinetic_mgrs, doublereal T, string output_file)
{
    vector<doublereal> hrxn;
    ofstream out;
    out.open(output_file);
    out << "#Dimensionless enthalpies of reactions (H/RT)\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            hrxn.resize(size);
            mgr->getDeltaEnthalpy(hrxn.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12);
                out << std::left << hrxn[k]/(GasConstant*T) ;
                out.width(12);
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
    out.close();
}

void print_rxn_eq_consts(vector<Kinetics*> kinetic_mgrs, string output_file)
{
    vector<doublereal> kc;
    ofstream out;
    out.open(output_file);
    out << "#Equilibrium constants of reactions\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            kc.resize(size);
            mgr->getEquilibriumConstants(kc.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12);
                out << std::left << kc[k];
                out.width(12);
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
    out.close();
}

}
