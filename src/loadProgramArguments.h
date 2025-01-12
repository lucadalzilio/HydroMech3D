//
// Created by Federico Ciardo on 26.07.21. All rights reserved.
//

// Inclusion from Standard Library
#include <fstream>

// Inclusion from InsideLoop library
#include <armadillo>

// Inclusion from the project
#include <nlohmann/json.hpp>

#ifndef INC_3DEQSIM_SRC_LOADPROGRAMARGUMENTS_H
#define INC_3DEQSIM_SRC_LOADPROGRAMARGUMENTS_H

// Using declaration for json library to lighten the code
using json = nlohmann::json;

namespace EQSim {
void LoadProgramArguments(int argc, const char *const *argv,
                          std::string &argFileName, bool &checkRestart, json &js);
}
#endif  // INC_3DEQSIM_SRC_LOADPROGRAMARGUMENTS_H
