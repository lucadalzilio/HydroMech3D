//
// Created by Federico Ciardo on 28.07.21. All rights reserved.
//

// Inclusion from the project
#include <nlohmann/json.hpp>

// Inclusion from InsideLoop library
#include <armadillo>

#ifndef INC_3DEQSIM_SRC_EQSOLVER_H
#define INC_3DEQSIM_SRC_EQSOLVER_H

// Using declaration for json library to lighten the code
using json = nlohmann::json;

namespace EQSim {
void EQsolver(json &js, bool &checkRestart);
}

#endif  // INC_3DEQSIM_SRC_EQSOLVER_H
