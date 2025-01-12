//
// Created by Federico Ciardo on 30.09.21. All rights reserved.
//
#include <armadillo>

#ifndef INC_3DEQSIM_SRC_FAULTPROPERTIES_H
#define INC_3DEQSIM_SRC_FAULTPROPERTIES_H

namespace EQSim {

class FrictionProperties;
class DilatancyProperties;
class PermeabilityProperties;

// This class encapsulates all the fault properties
// Class members:
// - fault gauge porosity [-]
// - void volume compressibility [1/Pa]
// - fault hydraulic aperture [m]

class FaultProperties {
 private:
  double porosity_;
  double void_compressibility_;
  arma::vec initial_hydraulic_aperture_;
  EQSim::PermeabilityProperties *permeabilityProperties_;
  EQSim::FrictionProperties *friction_properties_;
  arma::vec initial_DDs_;
  arma::vec initial_DDs_rates_;
  arma::vec initial_state_variables_;
  arma::vec initial_plastic_porosity_;
  arma::vec internal_loads_;

  // Constructors
 public:
  FaultProperties() = default;
  FaultProperties(double &porosity, double &void_compressibility,
                  arma::vec &initial_hydraulic_aperture,
                  EQSim::FrictionProperties *friction_properties,
                  EQSim::PermeabilityProperties *permeabilityProperties,
                  arma::vec &initial_DDs,
                  arma::vec &initial_DDs_rates,
                  arma::vec &initial_state_variables,
                  arma::vec &initial_plastic_porosity,
                  arma::vec &internal_loads) {
    porosity_ = porosity;
    void_compressibility_ = void_compressibility;
    initial_hydraulic_aperture_ = initial_hydraulic_aperture;
    friction_properties_ = friction_properties;
    permeabilityProperties_ = permeabilityProperties;
    initial_DDs_ = initial_DDs;
    initial_DDs_rates_ = initial_DDs_rates;
    initial_state_variables_ = initial_state_variables;
    initial_plastic_porosity_ = initial_plastic_porosity;
    internal_loads_ = internal_loads;
  }

  // Getter methods
  double getFaultPorosity() const { return porosity_; };
  double getInitialPlasticFaultPorosity(int index) const {
    return initial_plastic_porosity_[index];
  };
  double getFaultVoidCompressibility() const { return void_compressibility_; };
  double getInitialFaultHydraulicAperture(int i) const {
    return initial_hydraulic_aperture_[i];
  };
  PermeabilityProperties *getPermeabilityProperties() {
    return permeabilityProperties_;
  };
  FrictionProperties *getFrictionProperties() { return friction_properties_; };
  arma::vec getInitialDDs() { return initial_DDs_; };
  arma::vec getInitialDDsRates() { return initial_DDs_rates_; };
  arma::vec getInitialStateVariables() {
    return initial_state_variables_;
  };
  arma::vec getInternalLoads() { return internal_loads_; };
  arma::vec getAmbientPressureDistribution() {
    arma::vec ambient_pressure(internal_loads_.n_elem / 3, arma::fill::zeros);
    for (arma::uword I = 0; I < internal_loads_.n_elem / 3; ++I) {
      ambient_pressure[I] = internal_loads_[3 * I + 2];
    }
    return ambient_pressure;
  };

  // Setter methods
  void setFaultPorosity(double &porosity) { porosity_ = porosity; };
  void setFaultVoidCompressibility(double &void_compressibility) {
    void_compressibility_ = void_compressibility;
  };
};
}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_FAULTPROPERTIES_H
