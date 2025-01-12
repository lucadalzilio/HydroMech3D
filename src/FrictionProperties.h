//
// Created by Federico Ciardo on 30.08.21. All rights reserved.
//

// Inclusion from Standard library
#include <iostream>
#include <cassert>

// Import from the project
#include "Mesh.h"
#include "Solution.h"

#ifndef INC_3DEQSIM_SRC_FRICTION_H
#define INC_3DEQSIM_SRC_FRICTION_H

namespace EQSim {

// This class encapsulates all the friction properties
// It includes the virtual function (interface) that sets the friction
// parameters of a given friction law to each element of the mesh (based on the
// local MatID)

class FrictionProperties {
 protected:
  arma::vec friction_coefficients_;
  arma::vec a_values_;
  arma::vec b_values_;
  arma::vec reference_velocities_;
  arma::vec state_evolution_distances_;
  arma::vec reference_friction_coefficients_;
  arma::vec peak_friction_coefficients_;
  arma::vec residual_friction_coefficients_;
  arma::vec residual_slips_;

 public:
  // Constructor
  FrictionProperties() = default;

  // Interface
  virtual void setFrictionParameters(json &j_friction_param,
                                     EQSim::Mesh &Mesh) = 0;

  // Getter methods
  arma::vec getFrictionCoefficient() { return friction_coefficients_; }
  double getFrictionCoefficient(arma::uword i) {
    return friction_coefficients_[i];
  }
  arma::vec get_a_values() const { return a_values_; };
  double get_a_values(arma::uword &i) const { return a_values_[i]; };
  arma::vec get_b_values() const { return b_values_; };
  double get_b_values(arma::uword &i) const { return b_values_[i]; };
  arma::vec get_reference_velocities() {
    return reference_velocities_;
  };
  double get_reference_velocities(arma::uword &i) {
    return reference_velocities_[i];
  };
  arma::vec get_state_evolution_distances_() const {
    return state_evolution_distances_;
  };
  double get_state_evolution_distances_(arma::uword &i) const {
    return state_evolution_distances_[i];
  };
  arma::vec get_reference_friction_coefficients_() const {
    return reference_friction_coefficients_;
  };
  double get_reference_friction_coefficients_(arma::uword i) const {
    return reference_friction_coefficients_[i];
  };
};

// This sub-class inherits from FrictionProperties class and encapsulates all
// the parameters of the Rate and State friction law
// The first five members are used to import the input parameters from json
// input file.
// The last five members, instead, include the input parameters for
// all the mesh.

class RateAndStateFriction : public FrictionProperties {
 private:
  arma::vec a_;
  arma::vec b_;
  arma::vec Vo_;
  arma::vec Dc_;
  arma::vec fo_;

 public:
  void setFrictionParameters(json &j_friction_param,
                             EQSim::Mesh &Mesh) override {
    // Initial check of keywords
    if (j_friction_param.count("Reference friction coefficient") != 1) {
      std::cout
          << "Reference friction coefficient keyword is wrong in input file! "
          << std::endl;
      std::abort();
    }

    if (j_friction_param.count("a-values") != 1) {
      std::cout << "a-values keyword is wrong in input file! " << std::endl;
      std::abort();
    }

    if (j_friction_param.count("b-values") != 1) {
      std::cout << "b-values keyword is wrong in input file! " << std::endl;
      std::abort();
    }

    if (j_friction_param.count("State evolution distances") != 1) {
      std::cout << "State evolution distances keyword is wrong in input file! "
                << std::endl;
      std::abort();
    }

    if (j_friction_param.count("Reference velocities") != 1) {
      std::cout << "Reference velocities keyword is wrong in input file! "
                << std::endl;
      std::abort();
    }

    arma::vec a(1, arma::fill::zeros);
    arma::vec b(1, arma::fill::zeros);
    arma::vec Dc(1, arma::fill::zeros);
    arma::vec Vo(1, arma::fill::zeros);
    arma::vec fo(1, arma::fill::zeros);

    assert(fo.size() ==
           j_friction_param["Reference friction coefficient"].size());
    assert(a.size() == j_friction_param["a-values"].size());
    assert(b.size() == j_friction_param["b-values"].size());
    assert(Dc.size() ==
           j_friction_param["State evolution distances"].size());
    assert(Vo.size() == j_friction_param["Reference velocities"].size());

    for (arma::uword i = 0; i < 1; ++i) {
      fo[i] = j_friction_param["Reference friction coefficient"][i];
      a[i] = j_friction_param["a-values"][i];
      b[i] = j_friction_param["b-values"][i];
      Dc[i] = j_friction_param["State evolution distances"][i];
      Vo[i] = j_friction_param["Reference velocities"][i];
    }

    fo_ = fo;
    a_ = a;
    b_ = b;
    Dc_ = Dc;
    Vo_ = Vo;

    // Assign all the material parameters to each element in the mesh

    arma::uword Nelts = Mesh.getNumberOfElts();
    arma::vec foValues(Nelts, arma::fill::zeros);
    arma::vec aValues(Nelts, arma::fill::zeros);
    arma::vec bValues(Nelts, arma::fill::zeros);
    arma::vec reference_velocities(Nelts, arma::fill::zeros);
    arma::vec state_evolution_distances(Nelts, arma::fill::zeros);

    for (arma::uword I = 0; I < Nelts; ++I) {
      foValues[I] = fo_[0];
      aValues[I] = a_[0];
      bValues[I] = b_[0];
      reference_velocities[I] = Vo_[0];
      state_evolution_distances[I] = Dc_[0];
    }

    reference_friction_coefficients_ = foValues;
    a_values_ = aValues;
    b_values_ = bValues;
    reference_velocities_ = reference_velocities;
    state_evolution_distances_ = state_evolution_distances;
  };
};

}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_FRICTION_H
