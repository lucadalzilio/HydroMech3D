
// Inclusion from Standard library
#include <iostream>

// Inclusion from the project
#include <src/EQsolver.h>
#include <src/loadProgramArguments.h>
int main(const int argc, const char* const* argv) {
  std::string argFileName;
  bool checkRestart = false;
  json js;

  EQSim::LoadProgramArguments(argc, argv, argFileName, checkRestart,
                              js);

  EQSim::EQsolver(js, checkRestart);

  std::cout << std::endl;
  std::cout << "All good! The simulation is ended! " << std::endl;

  return 0;
}
