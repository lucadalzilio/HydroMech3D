//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

#ifndef INC_3DEQSIM_SRC_FLUIDPROPERTIES_H
#define INC_3DEQSIM_SRC_FLUIDPROPERTIES_H

namespace EQSim {

// This class encapsulates all the fluid properties
// Class members:
// - fluid density [kg/m^3]
// - fluid kinematic viscosity [Pa s]
// - fluid compressibility [1/Pa]

class FluidProperties {
 private:
  double density_;
  double viscosity_;
  double compressibility_;

  // Constructors
 public:
  FluidProperties() = default;
  FluidProperties(double &density, double &compressibility, double &viscosity) {
    density_ = density;
    viscosity_ = viscosity;
    compressibility_ = compressibility;
  }

  // Getter methods
  double getFluidDensity() const { return density_; };
  double getFluidViscosity() const { return viscosity_; };
  double getFluidCompressibility() const { return compressibility_; };

  // Setter methods
  void setFluidDensity(double &density) { density_ = density; };
  void setFluidViscosity(double &viscosity) { viscosity_ = viscosity; };
  void setFluidCompressibility(double &compressibility) {
    compressibility_ = compressibility;
  };
};
}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_FLUIDPROPERTIES_H
