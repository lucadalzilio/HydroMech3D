//
// Created by Federico Ciardo on 30.08.21. All rights reserved.
//

// Import from IL library
#include <armadillo>
#include <JSON/nlohmann/json.hpp>

#include "FaultInSituStress.h"
#include "FaultProperties.h"
#include "FluidProperties.h"
#include "Injection.h"
#include "Mesh.h"
#include "SolidMatrixProperties.h"
#include "SolverParameters.h"

// Using declaration for json library to lighten the code
using json = nlohmann::json;

#ifndef INC_3DEQSIM_SRC_SOLUTION_H
#define INC_3DEQSIM_SRC_SOLUTION_H

namespace EQSim {

// Forward declaration to avoid circular dependency between header files
class FrictionProperties;
class PermeabilityProperties;

class Solution {
 protected:
  double current_time_;
  double time_step_;
  EQSim::SolverParameters solver_parameters_;
  arma::vec insitu_tractions_;
  arma::vec total_tractions_;
  arma::vec DDs_;
  arma::vec DDs_rates_;
  arma::vec pressure_;
  arma::vec pressure_o_;
  EQSim::Injection injection_;
  arma::vec hydraulic_aperture_;
  EQSim::Mesh mesh_;

 public:
  Solution() = default;

  Solution(double &current_time, double &time_step,
           EQSim::SolverParameters &solver_parameters,
           arma::vec &insitu_tractions, arma::vec &DDs,
           arma::vec &DDs_rates, arma::vec &pressure,
           arma::vec &init_pressure,
           arma::vec &hydraulic_aperture, EQSim::Mesh &mesh) {
    mesh_ = mesh;
    current_time_ = current_time;
    time_step_ = time_step;
    solver_parameters_ = solver_parameters;
    insitu_tractions_ = insitu_tractions;
    DDs_rates_ = DDs_rates;
    DDs_ = DDs;
    pressure_ = pressure;
    pressure_o_ = init_pressure;
    hydraulic_aperture_ = hydraulic_aperture;
  }

  // Methods
  json createJsonObjectToExport();
  void exportBaseSolutionToJson(json &j_obj_to_export, std::string &filename);
  void exportBaseSolutionToUBJson(json &j_obj_to_export, std::string &filename);

  // Getter methods
  EQSim::Mesh getMesh() const { return mesh_; };
  EQSim::SolverParameters getSolverParametersStructure() {
    return solver_parameters_;
  };
  arma::vec getDDs() const { return DDs_; };
  double getDDs(arma::uword &i) const { return DDs_[i]; };
  arma::vec getSlip1() const {
    arma::vec slip1(DDs_.size() / 3, arma::fill::zeros);
    for (arma::uword I = 0; I < slip1.n_elem; ++I) {
      slip1[I] = DDs_[3 * I];
    }
    return slip1;
  };
  arma::vec getSlip2() const {
    arma::vec slip2(DDs_.size() / 3, arma::fill::zeros);
    for (arma::uword I = 0; I < slip2.n_elem; ++I) {
      slip2[I] = DDs_[3 * I + 1];
    }
    return slip2;
  };
  arma::vec getOpening() const {
    arma::vec opening(DDs_.size() / 3, arma::fill::zeros);
    for (arma::uword I = 0; I < opening.n_elem; ++I) {
      opening[I] = DDs_[3 * I + 2];
    }
    return opening;
  };
  arma::vec getDDsRates() const { return DDs_rates_; };
  double getDDsRates(arma::uword index) const { return DDs_rates_[index]; };
  arma::vec getSlipRate1() const {
    arma::vec slip_rate1(DDs_rates_.size() / 3, arma::fill::zeros);
    for (arma::uword I = 0; I < slip_rate1.n_elem; ++I) {
      slip_rate1[I] = DDs_rates_[3 * I];
    }
    return slip_rate1;
  };
  arma::vec getSlipRate2() const {
    arma::vec slip_rate2(DDs_rates_.size() / 3, arma::fill::zeros);
    for (arma::uword I = 0; I < slip_rate2.n_elem; ++I) {
      slip_rate2[I] = DDs_rates_[3 * I + 1];
    }
    return slip_rate2;
  };
  arma::vec getPressure() const { return pressure_; };
  double getPressure(arma::uword i) const { return pressure_[i]; };
  arma::vec getInitPressure() const { return pressure_o_; };
  double getInitPressure(arma::uword &i) const { return pressure_o_[i]; };
  arma::vec getInsituTractions() const { return insitu_tractions_; };
  double getInsituTractions(arma::uword &i) const {
    return insitu_tractions_[i];
  };
  arma::vec getNormalInSituTractions() const {
    arma::vec normal_tractions(insitu_tractions_.n_elem / 3, arma::fill::zeros);
    for (arma::uword I = 0; I < normal_tractions.n_elem; ++I) {
      normal_tractions[I] = insitu_tractions_[3 * I + 2];
    }
    return normal_tractions;
  };
  double getNormalInSituTractions(arma::uword &i) const {
    return insitu_tractions_[3 * i + 2];
  };
  double getNormalTotalTractions(arma::uword &i) const {
    return total_tractions_[3 * i + 2];
  };
  double getCurrentTime() const { return current_time_; };
  double getTimeStep() const { return time_step_; };
  arma::vec getHydraulicAperture() const {
    return hydraulic_aperture_;
  };
  double getHydraulicAperture(arma::uword &i) const {
    return hydraulic_aperture_[i];
  };
  double getSlipRateAlongSlipDirection(arma::uword &elm_i,
                                       arma::uword &slip_direction) {
    return DDs_rates_[3 * elm_i + slip_direction];
  }

  // Setter methods
  void setSlip1(double &slip1, arma::uword &i) { DDs_[3 * i] = slip1; }
  void setSlip2(double &slip2, arma::uword &i) { DDs_[3 * i + 1] = slip2; }
  void setOpening(double &opening, arma::uword &i) { DDs_[3 * i + 2] = opening; }
  void setSlipRate1(double &slip_rate1, arma::uword &i) {
    DDs_rates_[3 * i] = slip_rate1;
  }
  void setSlipRate2(double &slip_rate2, arma::uword &i) {
    DDs_rates_[3 * i + 1] = slip_rate2;
  }
  void setSlipAlongSlipDirection(arma::vec &slip,
                                 arma::uword &slip_direction) {
    for (arma::uword I = 0; I < DDs_.n_elem / 3; ++I) {
      DDs_[3 * I + slip_direction] = slip[I];
    }
  }
  void setSlipAlongSlipDirection(arma::uword &i, const double &slip,
                                 arma::uword &slip_direction) {
    DDs_[3 * i + slip_direction] = slip;
  }
  void setSlipRateAlongSlipDirection(arma::vec &slip_rate,
                                     arma::uword &slip_direction) {
    for (arma::uword I = 0; I < DDs_rates_.n_elem / 3; ++I) {
      DDs_rates_[3 * I + slip_direction] = slip_rate[I];
    }
  }
  void setSlipRateAlongSlipDirection(arma::uword &i, const double &slip_rate,
                                     arma::uword &slip_direction) {
    DDs_rates_[3 * i + slip_direction] = slip_rate;
  }
  void setTotalTractions(arma::vec &tot_tractions) {
    total_tractions_ = tot_tractions;
  }

  void setTotalTractions(arma::uword &elm_i, double &shear_stress_1,
                         double &shear_stress_2, double &normal_stress) {
    total_tractions_[3 * elm_i] = shear_stress_1;
    total_tractions_[3 * elm_i + 1] = shear_stress_2;
    total_tractions_[3 * elm_i + 2] = normal_stress;
  }
  void setCurrentTime(const double &curr_time) { current_time_ = curr_time; }
  void setTimeStep(const double &time_step) { time_step_ = time_step; }
  void setPressure(const arma::vec &pressure) { pressure_ = pressure; }
  void setPressure(arma::uword &i, const double &value) { pressure_[i] = value;} 

  void setInitialPressure(const arma::vec &init_pressure) {
    pressure_o_ = init_pressure;
  }
};

/// RK45 Solution -> subclass of Solution class
class SolutionRK45 : public Solution {
 protected:
  arma::vec state_variables_;
  arma::vec plastic_fault_porosity_;
  arma::vec fault_hydraulic_aperture_;
  EQSim::Injection injection_;
  EQSim::FluidProperties fluid_properties_;
  EQSim::FaultProperties fault_properties_;
  EQSim::SolidMatrixProperties solidMatrix_Properties_;
  FrictionProperties *fric_properties_;
  PermeabilityProperties *permeab_properties_;
  std::string results_path_, baseFileName_;
  arma::mat ElastMatrix_;
  arma::imat neigh_elts_;

 public:
  SolutionRK45() = default;
  SolutionRK45(double &current_time, double &time_step,
               EQSim::SolverParameters &solver_parameters,
               EQSim::FluidProperties &fluid_prop,
               EQSim::SolidMatrixProperties &solidMatrix_Properties,
               arma::vec &insitu_tractions,
               arma::vec &total_tractions, arma::vec &DDs,
               arma::vec &DDs_rates, arma::vec &state_variables,
               arma::vec &plastic_fault_porosity,
               arma::vec &fault_hydraulic_aperture,
               arma::vec &pressure_o, arma::vec &pressure,
               EQSim::Injection &InjectionObj,
               EQSim::FaultProperties &fault_properties, EQSim::Mesh &mesh,
               FrictionProperties *fric_prop,
               PermeabilityProperties *permeab_prop, std::string &res_path,
               std::string &baseFileName, arma::mat const &elast_matrix,
               arma::imat &neigh_elts) {
    mesh_ = mesh;
    current_time_ = current_time;
    time_step_ = time_step;
    solver_parameters_ = solver_parameters;
    fluid_properties_ = fluid_prop;
    solidMatrix_Properties_ = solidMatrix_Properties;
    insitu_tractions_ = insitu_tractions;
    total_tractions_ = total_tractions;
    DDs_ = DDs;
    DDs_rates_ = DDs_rates;
    state_variables_ = state_variables;
    plastic_fault_porosity_ = plastic_fault_porosity;
    fault_hydraulic_aperture_ = fault_hydraulic_aperture;
    pressure_o_ = pressure_o;
    pressure_ = pressure;
    injection_ = InjectionObj;
    fault_properties_ = fault_properties;
    fric_properties_ = fric_prop;
    permeab_properties_ = permeab_prop;
    results_path_ = res_path;
    baseFileName_ = baseFileName;
    ElastMatrix_ = elast_matrix;
    neigh_elts_ = neigh_elts;
  }

  // Methods
  json createJsonObjectToExport(bool &export_background_stresses);
  void exportToJson(json &j_obj_to_export, std::string &filename);
  void exportToUBJson(json &j_obj_to_export, std::string &filename);

  // Getter methods
  Injection getInjectionObj() { return injection_; };
  FluidProperties getFluidProperties() { return fluid_properties_; };
  SolidMatrixProperties getSolidMatrixProperties() {
    return solidMatrix_Properties_;
  };
  FaultProperties getFaultProperties() { return fault_properties_; };
  FrictionProperties *getFrictionPtr() { return fric_properties_; };
  PermeabilityProperties *getPermeabPtr() { return permeab_properties_; };
  std::string getResPath() { return results_path_; };
  std::string getBaseFileName() { return baseFileName_; };
  arma::mat getElastMatrix() { return ElastMatrix_; };
  arma::vec getStateVariables() { return state_variables_; };
  arma::vec getPlasticFaultPorosity() {
    return plastic_fault_porosity_;
  };
  arma::vec getFaultHydraulicAperture() {
    return fault_hydraulic_aperture_;
  };
  double getInsituTractions(arma::uword i) const { return insitu_tractions_[i]; };
  arma::imat getNeighElts() const { return neigh_elts_; };

  // Setter methods
  void setStateVariables(arma::vec &state_variables) {
    state_variables_ = state_variables;
  };
  void setPlasticFaultPorosity(arma::vec &plastic_fault_porosity) {
    plastic_fault_porosity_ = plastic_fault_porosity;
  };
  void setFaultHydraulicAperture(arma::vec &fault_hydraulic_aperture) {
    fault_hydraulic_aperture_ = fault_hydraulic_aperture;
  };
};

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_SOLUTION_H
