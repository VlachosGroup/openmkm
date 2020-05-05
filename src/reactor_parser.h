//! @file ReactorParser.h

// This file is part of OpenMKM. See License.txt in the top-level directory 
// for license and copyright information.

#ifndef OMKM_RCTR_PARSER_H
#define OMKM_RCTR_PARSER_H

#include <vector>
#include <string>
#include <memory>
#include <yaml-cpp/yaml.h>

#include "cantera/IdealGasMix.h"
#include "cantera/InterfaceLatInt.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/SurfLatIntPhase.h"

#include "omkmexceptions.h"


namespace OpenMKM 
{

/*
YAML::Node getChildNode(YAML::Node& p_nd, 
                        std::string p_name, 
                        std::vector<std::string> rev_descendants);
*/

enum RctrType {
    BATCH,
    CSTR,
    PFR_0D,
    PFR
};

enum OutputFormat {
    DAT,
    CSV
};

class ReactorParser {
public:
    ReactorParser() :
        m_T(0), m_P(0)
    {}

    ReactorParser(std::string rctr_file) {
        LoadFile(rctr_file);
    }

    void LoadFile(std::string rctr_file) {
        m_tube_nd = YAML::LoadFile(rctr_file);
        read_mandatory_nodes();
    }

    double temperature() {
        return m_T;
    }

    double T() {
        return m_T;
    }

    double pressure() {
        return m_P;
    }
    
    double P() {
        return m_P;
    }
    
    // Check for phases and read them
    bool gasPhaseDefined(std::string phase_filename);

    std::string getGasPhaseComposition();

    std::shared_ptr<Cantera::IdealGasMix> getGasPhase(
            std::string phase_filename);

    std::shared_ptr<Cantera::Solution> getGasSolution(
            std::string phase_filename);

    bool bulkPhaseDefined(std::string phase_filename);

    std::shared_ptr<Cantera::StoichSubstance> getBulkPhase(
            std::string phase_filename);

    std::shared_ptr<Cantera::Solution> getBulkSolution(
            std::string phase_filename);

    bool surfacePhasesDefined(std::string phase_filename);

    std::vector<std::shared_ptr<Cantera::InterfaceInteractions>> getSurfPhases(
            std::string phase_filename, 
            std::vector<Cantera::ThermoPhase*> gb_phases);

    std::vector<std::shared_ptr<Cantera::Solution>> getSurfaceSolutions(
            std::string phase_filename, 
            std::vector<std::shared_ptr<Cantera::Solution>>& gb_phases);

    std::vector<std::string> getSurfPhaseCompositions();

    //! Get Reactor Type
    RctrType getReactorType();

    // Get Reactor Dimensions

    //! PFR Length 
    double getLength();

    //! PFR Cross Section Area
    double getXCArea();

    //! Reactor volume for 0d reactors 
    double getVolume();

    //! # CSTRs to use for 0d PFR
    size_t getNodes();

    bool catalystAreaDefined();
    double getCatalystAbyV();

    // Inlet parameters
    //! 
    bool FlowRateDefined();
    bool MassFlowRateDefined();
    bool ResidenceTimeDefined();
    double getFlowRate();
    double getMassFlowRate();
    double getResidenceTime();

    // Reactor Operational Modes

    std::string getMode();
    std::map<double, double> getTProfile();
    double getWallHeatTransferCoeff();
    double getWallSpecificArea();
    double getExternalTemp();

    // User defined numerical solver options 
    bool tolerancesDefined();
    double get_atol();
    double get_rtol();
    bool solverInitStepSizeDefined();
    double getSolverInitStepSize();
    bool solverMaxStepsDefined();
    double getSolverMaxSteps();

    // Simulation output flags & parameters
    
    bool initStepDefined();
    bool logTransient();
    double getInitStep();
    std::string steppingType();
    double getEndTime();
    bool RPA();

    // Simulation flags for TPD
    double getTPDTempRamp();
    double getTPDEndTemp();

    // Sensitivity Analysis
    bool SAEnabled();
    double get_satol();
    double get_srtol();
    std::vector<std::string> getSAReactions();


    //! Parametric study 
    bool parametric_study_enabled(){
        return (isT_multi_input() || 
                isP_multi_input() || 
                isFR_multi_input());
    }

    bool isT_multi_input(); 
    bool isP_multi_input(); 
    bool isFR_multi_input(); 

    std::vector<double> Ts();
    std::vector<double> Ps();
    std::vector<double> FRs();

    OutputFormat printFormat();

    /*
    bool validate() { // TODO: Implement for one shot error checking
        return true;
    }
    */
protected:
    void read_mandatory_nodes();

private:
    YAML::Node m_tube_nd;
    // The below YAML nodes are child nodes of tube_nd, but 
    // defined explicitly for convenience
    YAML::Node m_rctr_nd;       
    YAML::Node m_simul_nd;
    YAML::Node m_phase_nd;
    YAML::Node m_inlet_nd;
    double m_T;         // Temperature
    double m_P;         // Pressure
    std::string m_gas_X;         // Gas Composition
    std::vector<std::string> m_surf_X;    // Surfaces Compositions
    OutputFormat m_output_format;
};

}

#endif
