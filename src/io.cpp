#include <iostream>
#include <fstream>
#include <iomanip>

#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"
#include "io.h"
#include "reactor_parser.h"

using namespace std;
using namespace Cantera;

namespace OpenMKM
{

static OutputFormat data_format;

void setOutputFormat(OutputFormat output_type)
{
    data_format = output_type;
}


// Prints the number of species in the mechanism
void print_species_number(vector<shared_ptr<ThermoPhase>> phases) 
{
    int total_no_species = 0;
    for  (const auto phase: phases){
        total_no_species += phase->nSpecies();
    }
    cout << "Total # of species: " << total_no_species << endl;

    int total_surf_species = 0;
    for  (const auto phase: phases){
        auto phase_type = phase->type();
        if (phase_type == "IdealGas") {
            cout << "Total # of  gas phase species: " << phase->nSpecies() 
                 << endl;
        }
        if (phase_type == "SurfCoverage" || phase_type == "Surf"){
            total_surf_species += phase->nSpecies();
        }
    }
    cout << "Total # of surface species: " << total_surf_species << endl;

}

//! Prints the data to conform to species input of reaction path 
//! visualisation code, RenView
void print_species(vector<shared_ptr<ThermoPhase>> phases, string output_file) 
{
    vector<string> species; 
    ofstream out (output_file);
    int sp_indx = 1;

    // Get all the elements in the system
    auto add_element = [] (string el, vector<string>& elements) {
        auto it = find(elements.begin(), elements.end(), el);
        if (it == elements.end()){
            elements.push_back(el);
        }
    };

    auto expandedComposition = [] (const compositionMap& species_comp, 
                                    const vector<string>& elements) {
        vector<double> compositions;
        for (const auto& el : elements) {
            auto it = species_comp.find(el);
            if (it == species_comp.end()) {
                compositions.push_back(0.0);
            } else {
                compositions.push_back(it->second);
            }
        }
        return compositions;
    };

    vector<string> allElements; 
    for  (const auto phase: phases){
        for (const auto& el : phase->elementNames()){
            add_element(el, allElements);
        }
    }

    // Print the header
    out << setw(16) << left << "Species_name" 
        << setw(10) << left << "Phase" 
        << setw(16) << left << "Surf_cov";
    for (const auto& el : allElements) 
        out << setw(4) <<  el;
    out << endl;

    vector_fp coverages;
    for  (auto phase: phases){
        // Get the phase type
        auto phase_type = phase->type();
        if (phase_type == "StoichSubstance")
            phase_type = "Bulk";
        if (phase_type == "IdealGas")
            phase_type = "Gas";
        if (phase_type == "SurfCoverage" || phase_type == "Surf")
            phase_type = "Surface";
        if (phase_type == "Surface") {
            coverages.resize(phase->nSpecies());
            dynamic_pointer_cast<SurfPhase>(phase)->getCoverages(coverages.data());
        }

        //out.width(16); 
        for (size_t k = 0; k < phase->nSpecies(); k++) {
            auto species = phase->species(k);
            out << setw(16) << left << species->name 
                << setw(10) << left  << phase_type;
            if (phase_type == "Surface"){
                out << setw(16) << coverages[k];
            } else {
                out << setw(16) << 0.0;
            }
            vector<double> comps = expandedComposition(species->composition, 
                                                       allElements);
            for (const auto& comp : comps)
                out << setw(4) << comp;
            out << endl;
        } 
    }   
}
void print_formation_enthalpy(vector<shared_ptr<ThermoPhase>> phases, string output_file) 
{
    vector<doublereal> hform; 
    ofstream out (output_file);
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
}

void print_formation_entropy(vector<shared_ptr<ThermoPhase>> phases, string output_file) 
{
    vector<doublereal> sform; 
    ofstream out (output_file);
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
}

void print_rxns(vector<Kinetics*> kinetic_mgrs, string output_file)
{
    ofstream out (output_file);
    //out << "#Dimensionless enthalpies of reactions (H/RT)\n" << endl;
    int rxn_indx = 1;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            for (size_t k = 0; k < size; k++) {
                //out.width(12);
                //out << std::left << rxn_indx++;
                out.width(12);
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
}

void print_rxn_enthalpy(vector<Kinetics*> kinetic_mgrs, doublereal T, string output_file)
{
    vector<doublereal> hrxn;
    ofstream out (output_file);
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
}

void print_rxn_entropy(vector<Kinetics*> kinetic_mgrs, string output_file)
{
    vector<doublereal> sRxn;
    ofstream out (output_file);
    out << "#Dimensionless entropies of reactions (S/R)\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            sRxn.resize(size);
            mgr->getDeltaSSEntropy(sRxn.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12);
                out << scientific << std::left << sRxn[k]/GasConstant ;
                out.width(12);
                out << scientific << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
}


void print_rxn_eq_consts(vector<Kinetics*> kinetic_mgrs, string output_file)
{
    vector<doublereal> kc;
    ofstream out (output_file);
    out << "#Equilibrium constants of reactions\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            kc.resize(size);
            mgr->getEquilibriumConstants(kc.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12);
                out << scientific << std::left << kc[k];
                out.width(12);
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
}

void print_rxn_gibbs(vector<Kinetics*> kinetic_mgrs, doublereal T,  
        string output_file)
{
    vector<doublereal> muRxn; 
    ofstream out (output_file);
    out << "#Dimensionless Gibbs Energies of reactions (G/RT)\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            muRxn.resize(size);
            mgr->getDeltaSSGibbs(muRxn.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12); 
                out << scientific << std::left << muRxn[k]/(GasConstant*T) ;
                out.width(12); 
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
}


void print_rxn_kc(vector<Kinetics*> kinetic_mgrs, string output_file)
{
    vector<doublereal> kc; 
    ofstream out (output_file);
    out << "#Equilibrium constants of reactions\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            kc.resize(size);
            mgr->getEquilibriumConstants(kc.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12); 
                out << scientific << std::left << kc[k];
                out.width(12); 
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
}

void print_rxn_kf(vector<Kinetics*> kinetic_mgrs, string output_file)
{
    vector<doublereal> kf; 
    ofstream out (output_file);
    out << "#Forward rate constants of reactions\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            kf.resize(size);
            mgr->getFwdRateConstants(kf.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12); 
                out << scientific << std::left << kf[k];
                out.width(12); 
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
}

void print_rxn_kr(vector<Kinetics*> kinetic_mgrs, string output_file)
{
    vector<doublereal> kr; 
    ofstream out (output_file);
    out << "#Reverse rate constants of reactions\n" << endl;
    for  (const auto mgr: kinetic_mgrs){
        size_t size = mgr->nReactions();
        if (size > 0) {
            kr.resize(size);
            mgr->getRevRateConstants(kr.data());
            for (size_t k = 0; k < size; k++) {
                out.width(12); 
                out << scientific << std::left << kr[k];
                out.width(12); 
                out << std::left << mgr->reactionString(k) << endl;
            }
        }
    }
}

void print_omkm_header(std::ostream& out) {
    out << "-----------------------------------------------------------\n" 
        << "OpenMKM: version 0.3.0\n" 
        << "-----------------------------------------------------------\n\n" 
        << "OpenMKM is a multiphysics, multiscale, and open source software \n"
        << "aimed at modelng chemical kinetics. It can run pure gas phase \n"
        << "as well as surface mechanisms for heterogeneous catalysis.\n"
        << "OpenMKM is developed at Delaware Energy Institute, University\n"
        << "of Delaware. The development of OpenMKM is funded by RAPID.\n\n\n"; 
}

/**
 * Utility function to print reaction rates
 */
void print_rxn_rates(Kinetics* kin, ofstream& out)
{
    vector<double> fwd_rts;
    vector<double> rev_rts;
    vector<double> net_rts;
    auto nRxns = kin->nReactions();
    if (nRxns) {
        fwd_rts.resize(nRxns);
        rev_rts.resize(nRxns);
        net_rts.resize(nRxns);
        kin->getFwdRatesOfProgress(fwd_rts.data());
        kin->getRevRatesOfProgress(rev_rts.data());
        kin->getNetRatesOfProgress(net_rts.data());
        for (size_t i = 0; i < nRxns; i++) {
            auto pe = fwd_rts[i]/(fwd_rts[i] + rev_rts[i]);
            out << scientific
                << setw(16) << left << fwd_rts[i]
                << setw(16) << left << rev_rts[i]
                << setw(16) << left << net_rts[i]
                << setw(16) << left << pe
                << kin->reactionString(i) << endl;
        }
    }
}

/**
 * Utility function to print header before printing reaction rates
 */
void print_rxn_rates_hdr(ofstream& out)
{
    out //<< hdr << endl
        //<< setw(16) << left << "Reaction No."
        << setw(16) << left << "Fwd Rate"
        << setw(16) << left << "Rev Rate"
        << setw(16) << left << "Net Rate"
        << setw(16) << left << "Partial Equil."
        << setw(32) << left << "Reaction String"
        << endl;
}


/**
 * Utility function to print PFR state at given distance z from inlet
 */
void print_pfr_rctr_state(double z, PFR1d* rctr, vector<SurfPhase*> surfaces,
                          ofstream& gas_mole_out, ofstream& gas_mass_out,
                          ofstream& gas_msdot_out, ofstream& surf_cov_out,
                          ofstream& state_var_out)
{

    vector<double> work(rctr->contents().nSpecies());

    state_var_out << scientific;
    if (data_format == OutputFormat::CSV) {
        state_var_out << z << ","
                      << rctr->contents().temperature() << ","
                      << rctr->contents().pressure() << ","
                      << rctr->contents().density() 
                      << endl;
    } else {
        state_var_out << setw(16) << left  << z
                      << setw(16) << left  << rctr->contents().temperature()
                      << setw(16) << left  << rctr->contents().pressure()
                      << setw(16) << left  << rctr->contents().density()
                      << endl;
    }

    auto gas_var_print = [&work, &z](ostream& out) -> void
    {
        out << scientific;
        if (data_format == OutputFormat::CSV) {
            out << z;
            for (size_t k = 0; k < work.size(); k++) {
                out << "," << work[k];
            }
        } else {
            out << setw(16) << left << z;
            for (size_t k = 0; k < work.size(); k++) {
                out << setw(16) << left << work[k];
            }
        }
        out << endl;
    };

    rctr->contents().getMoleFractions(work.data());
    gas_var_print(gas_mole_out);

    rctr->contents().getMassFractions(work.data());
    gas_var_print(gas_mass_out);

    rctr->getSurfaceProductionRates(work.data());
    gas_var_print(gas_msdot_out);
    
    surf_cov_out << scientific;
    if (data_format == OutputFormat::CSV) {
        surf_cov_out << z;
        for (size_t j = 0;  j <  surfaces.size(); j++) {
            work.resize(surfaces[j]->nSpecies());
            rctr->surface(j)->getCoverages(work.data());
            for (size_t k = 0; k < work.size(); k++) {
                surf_cov_out << "," << work[k];
            }
        }
    } else {
        surf_cov_out << setw(16) << left << z;
        for (size_t j = 0;  j <  surfaces.size(); j++) {
            work.resize(surfaces[j]->nSpecies());
            rctr->surface(j)->getCoverages(work.data());
            for (size_t k = 0; k < work.size(); k++) {
                surf_cov_out << setw(16) << left << work[k];
            }
        }
    }
    surf_cov_out << endl;
}

void print_0d_rctr_state(double z, Reactor* rctr, vector<SurfPhase*> surfaces,
                         ofstream& gas_mole_out, ofstream& gas_mass_out,
                         ofstream& gas_msdot_out, ofstream& surf_cov_out,
                         ofstream& state_var_out)
{

    vector<double> work(rctr->contents().nSpecies());

    state_var_out << scientific;
    if (data_format == OutputFormat::CSV) {
        state_var_out << z << ","
                      << rctr->temperature() << ","
                      << rctr->pressure()  << ","
                      << rctr->density()  << ","
                      << rctr->intEnergy_mass() 
                      << endl;
    }  else { //if (data_format == OutputFormat::DAT) {
        state_var_out << setw(16) << left  << z
                      << setw(16) << left  << rctr->temperature()
                      << setw(16) << left  << rctr->pressure() 
                      << setw(16) << left  << rctr->density() 
                      << setw(16) << left  << rctr->intEnergy_mass() 
                      << endl;
    } 

    auto gas_var_print = [&work, &z](ostream& out) -> void
    {
        out << scientific;
        if (data_format == OutputFormat::CSV) {
            out << z;
            for (size_t k = 0; k < work.size(); k++) {
                out << "," << work[k];
            }
        } else { //if (data_format == OutputFormat::DAT) {
            out << setw(16) << left << z;
            for (size_t k = 0; k < work.size(); k++) {
                out << setw(16) << left << work[k];
            }
        } 
        out << endl;
    };

    rctr->contents().getMoleFractions(work.data());
    gas_var_print(gas_mole_out);

    rctr->contents().getMassFractions(work.data());
    gas_var_print(gas_mass_out);

    rctr->getSurfaceProductionRates(work.data());
    gas_var_print(gas_msdot_out);

    surf_cov_out << scientific;
    if (data_format == OutputFormat::CSV) {
        surf_cov_out << z;
        for (size_t j = 0;  j <  surfaces.size(); j++) {
            work.resize(surfaces[j]->nSpecies());
            rctr->surface(j)->getCoverages(work.data());    
            for (size_t k = 0; k < work.size(); k++) {
                surf_cov_out << "," << work[k];
            }
        }
    } else { //if (data_format == OutputFormat::DAT) {
        surf_cov_out << setw(16) << left << z;
        for (size_t j = 0;  j <  surfaces.size(); j++) {
            work.resize(surfaces[j]->nSpecies());
            rctr->surface(j)->getCoverages(work.data());    
            for (size_t k = 0; k < work.size(); k++) {
                surf_cov_out << setw(16) << left << work[k];
            }
        }
    }
    surf_cov_out << endl;
}


}
