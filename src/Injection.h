//
// Created by Federico Ciardo on 01.10.21. All rights reserved.
//

// Inclusion from Standard library
#include <iostream>
#include <cassert>
// Import from the project
#include <JSON/nlohmann/json.hpp>

#include "ElementData.h"
#include "FluidProperties.h"
#include "Mesh.h"
#include <armadillo>

#ifndef INC_3DEQSIM_SRC_INJECTION_H
#define INC_3DEQSIM_SRC_INJECTION_H

// Using declaration for json library to lighten the code
using json = nlohmann::json;

namespace EQSim {

// This class encapsulates all the injection parameters
// Class members:
// - mass rates [Kg/s]
// - constant overpressure [Pa]
// - time at which injection change [s]
// - source elements, i.e. elements in which fluid is injected into
// - type of injection

class Injection {
 private:
  arma::vec mass_rates_source_elm_;
  arma::vec constant_overpressures_source_elm_;
  arma::vec time_at_which_injection_change_source_elm_;
  arma::uword source_element_;

  std::string type_of_injection_;

 public:
  // Constructor
  Injection() = default;

  Injection(json &j_injection) {
    if (j_injection.count("Type") == 1) {
      type_of_injection_ = j_injection["Type"].get<std::string>();
    } else {
      std::cout << "No injection type in input data!";
      std::abort();
    }

    if (j_injection.count("Source element") == 1) {
      source_element_ = j_injection["Source element"].get<arma::uword>();
    } else {
      std::cout << "No Source element in input data!";
      std::abort();
    }

    if (j_injection.count("Time at which injection change - Source Elem") ==
        1) {
      time_at_which_injection_change_source_elm_.resize(
          j_injection["Time at which injection change - Source Elem"].size());
      for (arma::uword I = 0;
           I < time_at_which_injection_change_source_elm_.n_elem; ++I) {
        time_at_which_injection_change_source_elm_[I] =
            j_injection["Time at which injection change - Source Elem"][I];
      }
    } else {
      std::cout << "No Time at which injection change - Source Elem in "
                   "input data!";
      std::abort();
    }

    if (type_of_injection_ == "Rate control") {
      if (j_injection.count("Mass rate - Source Element") == 1) {
        mass_rates_source_elm_.resize(
            j_injection["Mass rate - Source Element"].size());
        for (arma::uword I = 0; I < mass_rates_source_elm_.n_elem; ++I) {
          mass_rates_source_elm_[I] =
              j_injection["Mass rate - Source Element"][I];
        }
      } else {
        std::cout << "No Mass rate - Source Element in input data!";
        std::abort();
      }

      assert(time_at_which_injection_change_source_elm_.n_elem ==
                mass_rates_source_elm_.n_elem);

    } else if (type_of_injection_ == "Pressure control") {
      if (j_injection.count("Constant overpressure - Source Element") == 1) {
        constant_overpressures_source_elm_.resize(
            j_injection["Constant overpressure - Source Element"].size());
        for (arma::uword I = 0; I < constant_overpressures_source_elm_.n_elem;
             ++I) {
          constant_overpressures_source_elm_[I] =
              j_injection["Constant overpressure - Source Element"][I];
        }
      } else {
        std::cout
            << "No constant overpressure - Source Element in input data!";
        std::abort();
      }

      assert(time_at_which_injection_change_source_elm_.n_elem ==
                constant_overpressures_source_elm_.n_elem);
    }
  }

  // Getter methods
  arma::vec getTimeAtWhichInjectionChangeSourceElem() {
    return time_at_which_injection_change_source_elm_;
  };
  double getTimeAtWhichInjectionChangeSourceElem(arma::uword &i) {
    return time_at_which_injection_change_source_elm_[i];
  };
  arma::uword getSourceElement1() const { return source_element_; };
  arma::vec getMassRatesSourceElem() {
    return mass_rates_source_elm_;
  };
  arma::vec getConstantOverpressuresSourceElem() {
    return constant_overpressures_source_elm_;
  };

  // Setter methods
  void setTimeAtWhichInjectionChangeSourceElem(arma::uword &i, double &time) {
    time_at_which_injection_change_source_elm_[i] = time;
  };
  void setMassRateSourceElem(arma::uword &i, double &mass_rate) {
    mass_rates_source_elm_[i] = mass_rate;
  };
  void setConstantOverpressureSourceElem(arma::uword &i, double &const_dp) {
    constant_overpressures_source_elm_[i] = const_dp;
  };

  // Methods
  arma::vec calculate_Q_vector(EQSim::Mesh &Mesh,
                                       EQSim::FluidProperties &FluidProperties,
                                       double &curr_time) {
    arma::uword Nelts = Mesh.getNumberOfElts();
    arma::vec Q_vector(Nelts, arma::fill::zeros);
    arma::uword source_elt;
    arma::vec new_time_at_which_injection_change{};
    arma::vec new_mass_rate{};

    // If the type of injection is rate control, then calculates the Q vector
    // with the proper mass rates, otherwise keep the Q vector null
    if (type_of_injection_ == "Rate control") {
      if ((time_at_which_injection_change_source_elm_.n_elem >= 2) &&
          (curr_time >= time_at_which_injection_change_source_elm_[0])) {
        new_time_at_which_injection_change.resize(
            time_at_which_injection_change_source_elm_.n_elem - 1);
        new_mass_rate.resize(
            time_at_which_injection_change_source_elm_.n_elem - 1);
        for (arma::uword I = 0; I < new_time_at_which_injection_change.n_elem;
             ++I) {
          new_time_at_which_injection_change[I] =
              time_at_which_injection_change_source_elm_[I + 1];
          new_mass_rate[I] = mass_rates_source_elm_[I + 1];
        }

        time_at_which_injection_change_source_elm_ =
            new_time_at_which_injection_change;
        mass_rates_source_elm_ = new_mass_rate;
      }

      EQSim::ElementData DataSourceElem1 =
          Mesh.getElementData(source_element_);
      for (arma::uword I = 0; I < mass_rates_source_elm_.n_elem; ++I) {
        Q_vector[source_element_] = (mass_rates_source_elm_[I] /
                                      (FluidProperties.getFluidDensity() *
                                       DataSourceElem1.getAreaElt()));  // m / s
      }
    }

    return Q_vector;
  }
};

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_INJECTION_H
