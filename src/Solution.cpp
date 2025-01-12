//
// Created by Federico Ciardo on 22.09.21. All rights reserved.
//

// Inclusion from the standard library
#include <fstream>
#include <iomanip>

// Inclusion from the project
#include <src/FrictionProperties.h>
#include <src/PermeabilityProperties.h>

#include "Injection.h"

namespace EQSim {

// Implementation of methods that export to Json/Ubjson file the current
// solution

/////////////////////////    BaseSolution     /////////////////////////

json Solution::createJsonObjectToExport() {
  json json_coord = json::array();
  for (arma::uword m = 0; m < mesh_.getCoordinates().n_rows; ++m) {
    json_coord[m] = {mesh_.getCoordinates(m, 0), mesh_.getCoordinates(m, 1),
                     mesh_.getCoordinates(m, 2)};
  }
  //  connectivity, and dofs array
  json json_connectivity = json::array();
  for (arma::uword m = 0; m < mesh_.getNumberOfElts(); ++m) {
    json_connectivity[m] = {
        mesh_.getNodesConnettivity(m, 0), mesh_.getNodesConnettivity(m, 1),
        mesh_.getNodesConnettivity(m, 2), mesh_.getNodesConnettivity(m, 3)};
  }

  json json_source_elt1 = json::array();
  auto source_elt1 = injection_.getSourceElement1();
  json_source_elt1[0] = source_elt1;

  json json_mass_rates = json::array();
  auto mass_rates = injection_.getMassRatesSourceElem();
  for (arma::uword I = 0; I < mass_rates.n_elem; ++I) {
    json_mass_rates[I] = mass_rates[I];
  }

  json json_const_overpressure = json::array();
  auto const_overpress = injection_.getConstantOverpressuresSourceElem();
  for (arma::uword I = 0; I < const_overpress.n_elem; ++I) {
    json_const_overpressure[I] = const_overpress[I];
  }

  //  json json_tau1 = json::array();
  //  json json_tau2 = json::array();
  //  json json_sigma_n = json::array();
  //  json json_slip1 = json::array();
  //  json json_slip2 = json::array();
  //  json json_slip_rate1 = json::array();
  //  json json_slip_rate2 = json::array();
  //  json json_opening = json::array();
  json json_pressure = json::array();
  json json_pressure_o = json::array();

  for (arma::uword m = 0; m < mesh_.getNumberOfElts(); ++m) {
    //    json_tau1[m] = tau1_[m];
    //    json_tau2[m] = tau2_[m];
    //    json_slip1[m] = slip1_[m];
    //    json_slip2[m] = slip2_[m];
    //    json_slip_rate1[m] = slip_rate1_[m];
    //    json_slip_rate2[m] = slip_rate2_[m];
    //    json_opening[m] = opening_[m];
    json_pressure[m] = pressure_[m];
    json_pressure_o[m] = pressure_o_[m];
  }

  json j_mesh = {{"Interpolation order", mesh_.getInterpolationOrder()},
                 {"Node coordinates", json_coord},
                 {"Connectivity", json_connectivity}};

  json j_obj = {{"Current Time", current_time_},
                {"Time step", time_step_},
                {"Mesh", j_mesh},
                {"Source element 1", json_source_elt1},
                {"Mass rate", json_mass_rates},
                {"Constant overpressure", json_const_overpressure},
                //                {"Tau1", json_tau2},
                //                {"Tau2", json_tau2},
                //                {"Slip1", json_slip1},
                //                {"Slip2", json_slip2},
                //                {"Slip_Rate1", json_slip_rate1},
                //                {"Slip_Rate2", json_slip_rate2},
                //                {"Opening", json_opening},
                {"Pressure", json_pressure},
                {"Initial_Pressure", json_pressure_o}};

  return j_obj;
};

void Solution::exportBaseSolutionToJson(json &j_obj_to_export,
                                        std::string &filename) {
  std::ofstream output(filename + std::string{".json"});
  output << std::setw(4) << j_obj_to_export << std::endl;
};

void Solution::exportBaseSolutionToUBJson(json &j_obj_to_export,
                                          std::string &filename) {
  // Move to ubjson object
  std::vector<std::uint8_t> j_ubjson_obj = json::to_ubjson(j_obj_to_export);
  // write prettified binary file
  std::ofstream output;
  output.open(filename + std::string{".ubj"}, std::ios::out | std::ios::binary);
  output.write((char *)&j_ubjson_obj[0],
               j_ubjson_obj.size() * sizeof(j_ubjson_obj[0]));
  output.close();
};

// Implementation of methods that export to Json/Ubjson file the current
// solution object

/////////////////////////    RK45 Solution     /////////////////////////

json SolutionRK45::createJsonObjectToExport(bool &export_background_stresses) {
  json json_source_elt1 = json::array();
  arma::uword source_elt1 = injection_.getSourceElement1();
  json_source_elt1 = source_elt1;

  json json_mass_rates = json::array();
  auto mass_rates = injection_.getMassRatesSourceElem();
  for (arma::uword I = 0; I < mass_rates.n_elem; ++I) {
    json_mass_rates[I] = mass_rates[I];
  }

  json json_tau1 = json::array();
  json json_tau2 = json::array();
  json json_sigma_n = json::array();
  json json_slip1 = json::array();
  json json_slip2 = json::array();
  json json_slip_rate1 = json::array();
  json json_slip_rate2 = json::array();
  json json_opening = json::array();
  json json_cumul_incrm_opening = json::array();
  json json_pressure = json::array();
  json json_state_variable = json::array();
  json json_plastic_fault_porosity = json::array();

  arma::uword Nelts = mesh_.getNumberOfElts();

  for (arma::uword m = 0; m < Nelts; ++m) {
    json_tau1[m] = total_tractions_[3 * m];
    json_tau2[m] = total_tractions_[3 * m + 1];
    json_sigma_n[m] = total_tractions_[3 * m + 2];
    json_slip1[m] = DDs_[3 * m];
    json_slip2[m] = DDs_[3 * m + 1];
    json_cumul_incrm_opening[m] = DDs_[3 * m + 2];
    json_opening[m] = fault_hydraulic_aperture_[m];
    json_slip_rate1[m] = DDs_rates_[3 * m];
    json_slip_rate2[m] = DDs_rates_[3 * m + 1];
    json_pressure[m] = pressure_[m];
    json_state_variable[m] = state_variables_[m];
    json_plastic_fault_porosity[m] = plastic_fault_porosity_[m];
  }

  json j_obj = {{"Current Time", current_time_},
                {"Time step", time_step_},
                {"Source element 1", json_source_elt1},
                {"Mass rate", json_mass_rates},
                {"Tau1", json_tau1},
                {"Tau2", json_tau2},
                {"SigmaN", json_sigma_n},
                {"Slip1", json_slip1},
                {"Slip2", json_slip2},
                {"Opening", json_opening},
                {"Cumul. increment of opening", json_opening},
                {"Slip_Rate1", json_slip_rate1},
                {"Slip_Rate2", json_slip_rate2},
                {"Pressure", json_pressure},
                {"State variable", json_state_variable},
                {"Plastic fault porosity", json_plastic_fault_porosity}};

  return j_obj;
};

void SolutionRK45::exportToJson(json &j_obj_to_export, std::string &filename) {
  std::ofstream output(filename + std::string{".json"});
  output << std::setw(4) << j_obj_to_export << std::endl;
};

void SolutionRK45::exportToUBJson(json &j_obj_to_export,
                                  std::string &filename) {
  // Move to ubjson object
  std::vector<std::uint8_t> j_ubjson_obj = json::to_ubjson(j_obj_to_export);
  // write prettified binary file
  std::ofstream output;
  output.open(filename + std::string{".ubj"}, std::ios::out | std::ios::binary);
  output.write((char *)&j_ubjson_obj[0],
               j_ubjson_obj.size() * sizeof(j_ubjson_obj[0]));
  output.close();
};

}  // namespace EQSim
