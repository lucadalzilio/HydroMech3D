//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

#ifndef INC_3DEQSIM_SRC_SOLIDMATRIXPROPERTIES_H
#define INC_3DEQSIM_SRC_SOLIDMATRIXPROPERTIES_H

namespace EQSim {

// This class encapsulates all the solid matrix properties
// Class members:
// - Young's modulus [Pa]
// - Shear modulus [Pa]
// - Poisson's ratio [-]
// - Rock density [Kg/m^3]

class SolidMatrixProperties {
 private:
  double young_;
  double shear_;
  double poiss_;
  double density_;

 public:
  // Constructors
  SolidMatrixProperties() = default;

  SolidMatrixProperties(double& YoungModulus, double& PoissonRatio,
                        double& RockDensity) {
    young_ = YoungModulus;
    poiss_ = PoissonRatio;
    shear_ = young_ / (2.0 * (1.0 + poiss_));
    density_ = RockDensity;
  };

  // Getter methods
  double getYoungModulus() const { return young_; }
  double getShearModulus() const { return shear_; }
  double getPoissonRatio() const { return poiss_; }
  double getRockDensity() const { return density_; }
  double getBulkModulus() const { return young_ / (3. * (1. - 2. * poiss_)); }
  };
}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_SOLIDMATRIXPROPERTIES_H
