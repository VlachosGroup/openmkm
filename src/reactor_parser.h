#ifndef OMKM_RCTR_PARSER_H
#define OMKM_RCTR_PARSER_H

#include <vector>
#include <string>
#include <memory>
#include <yaml-cpp/yaml.h>

#include "cantera/IdealGasMix.h"
#include "cantera/InterfaceLatInt.h"
//#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/SurfLatIntPhase.h"

#include "omkmexceptions.h"


namespace OpenMKM 
{

YAML::Node getChildNode(YAML::Node& p_nd, 
                        std::string p_name, 
                        std::vector<std::string> rev_descendants);

enum RctrType {
    BATCH,
    CSTR,
    PFR_0D,
    PFR
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

    void read_mandatory_nodes();
    
    // Check for phases and read them
    bool GasPhaseDefined(std::string phase_filename);

    std::shared_ptr<Cantera::IdealGasMix> getGasPhase(
            std::string phase_filename);

    bool BulkPhaseDefined(std::string phase_filename);

    std::shared_ptr<Cantera::StoichSubstance> getBulkPhase(
            std::string phase_filename);

    bool SurfacePhasesDefined(std::string phase_filename);

    std::vector<std::shared_ptr<Cantera::InterfaceInteractions>> getSurfPhases(
            std::string phase_filename, 
            std::vector<Cantera::ThermoPhase*> gb_phases);

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

    double getInitStep();

    bool RPA();



    //! Parametric study 
    bool parametric_study_enabled(){
        return (T_parametric_study_enabled() || 
                P_parametric_study_enabled() || 
                mdot_parametric_study_enabled());
    }

    bool T_parametric_study_enabled(); 
    bool P_parametric_study_enabled(); 
    bool mdot_parametric_study_enabled(); 

    std::vector<double> T_parameter_samples();
    std::vector<double> P_parameter_samples();
    std::vector<double> mdot_parameter_samples();

    /*
    bool validate() { // TODO: Implement for one shot error checking
        return true;
    }
    */

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
};

}

#endif
