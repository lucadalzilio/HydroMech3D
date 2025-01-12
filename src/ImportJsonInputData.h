//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

#ifndef INC_3DEQSIM_SRC_IMPORTJSONINPUTDATA_H
#define INC_3DEQSIM_SRC_IMPORTJSONINPUTDATA_H

// Inclusion from Standard library
#include <iostream>

// Inclusion from InsideLoop library
#include <armadillo>

// Import from the project
#include <JSON/nlohmann/json.hpp>

#include "FaultInSituStress.h"
#include "FaultProperties.h"
#include "FluidProperties.h"
#include "Mesh.h"
#include "SolidMatrixProperties.h"
#include "SolverParameters.h"

// Using declaration for json library to lighten the code
using json = nlohmann::json;

namespace EQSim {

void ImportFirstLayerJsonInput(json &js, bool &checkRestart, 
                               std::string &solver_description,
                               std::string &basefilename, std::string &date,
                               std::string &program_name, std::string &res_path,
                               json &j_mesh, json &j_model_params,
                               json &j_solver_params);

void ImportSecondLayerJsonInput(json &j_model_params, bool &checkRestart,
                                json &j_injection,
                                json &j_fluid_params, json &j_fault_properties,
                                json &j_fault_insitu_params,
                                json &j_initial_conditions,
                                json &j_rock_properties, json &j_fluid_flow);

void ImportThirdLayerJsonInput(bool &checkRestart, json &j_fault_properties,
                               json &j_friction_properties,
                               json &j_permeability_properties);

EQSim::FluidProperties LoadFluidProperties(json &j_fluid_params);

EQSim::FaultProperties LoadFaultProperties(EQSim::Mesh &Mesh,
                                           json &j_fault_properties,
                                           json &j_friction_properties,
                                           json &j_permeability_properties,
                                           json &j_initial_conditions);

EQSim::SolidMatrixProperties LoadMatrixProperties(json &j_rock_properties);

EQSim::SolverParameters LoadSolverParameters(json &j_solver_params);

EQSim::Mesh LoadMeshData(json &j_mesh);

EQSim::FaultInSituStress LoadFaultInSituStressComponents(json &j_fault_insitu,
                                                         EQSim::Mesh &Mesh);

}  // namespace EQSim
#endif  // INC_3DEQSIM_SRC_IMPORTJSONINPUTDATA_H
