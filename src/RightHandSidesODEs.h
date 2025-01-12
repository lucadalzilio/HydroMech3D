//
// Created by Federico Ciardo on 31.08.21. All rights reserved.
//

// Import from IL library
#include <armadillo>

// Import from the project
#include "Solution.h"

#ifndef INC_3DEQSIM_SRC_RIGHTHANDSIDESODES_H
#define INC_3DEQSIM_SRC_RIGHTHANDSIDESODES_H

namespace EQSim {
class RightHandSideODEs {
  friend class RK45;

 private:
  double (*RhsFunc1_)(double &permeab,
                      double &fluid_compress, double &viscosity,
                      double &hydraulic_aperture, double &void_compress,
                      double &porosity, double &d2_p, double &Q);
  double (*RhsFunc2_)(double &theta, double &dc, double &slip_rate);
  long double (*RhsFunc3_)(double &ref_slip_velocity, double &a_value,
                           double &b_value, double &theta, double &dc,
                           double &ref_fric_coeff, double &sigma_n, double &p,
                           double &slip_rate, double &GtimesSlipRate,
                           double &permeab, double &fluid_compress,
                           double &viscosity, double &hydraulic_aperture,
                           double &void_compress, double &porosity,
                           double &d2_p, double &Q);

 public:
  RightHandSideODEs(
      double (*rhsFunc1)(double &, double &, double &, double &, double &,
                         double &, double &, double &),
      double (*rhsFunc2)(double &, double &, double &),
      long double (*rhsFunc3)(double &, double &, double &, double &, double &,
                              double &, double &, double &, double &, double &,
                              double &, double &, double &, double &, double &,
                              double &, double &, double &)) {
    RhsFunc1_ = rhsFunc1;
    RhsFunc2_ = rhsFunc2;
    RhsFunc3_ = rhsFunc3;
  }
};
}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_RIGHTHANDSIDESODES_H
