// Created by Federico Ciardo on 26.07.21. All rights reserved.

// Inclusion from Standard Library
#include <unistd.h>

// Inclusion from the project
#include "loadProgramArguments.h"

namespace EQSim {

void LoadProgramArguments(int argc, const char *const *argv, 
                    std::string &argFileName, bool &checkRestart, json &js) {
    // First check the number of passed arguments to the program
    if (argc != 3) {
        // If the number of arguments is not equal to 3, throw an error
        // and close the program.
        // Remember: the first argument is always the program name!
        std::cout << "Program arguments must be passed in this way: " << std::endl;
        std::cout << "-i inputfile ---> for new analyses;" << std::endl;
        std::cout << "-r restartfile ---> for restarted analyses." << std::endl;
        std::cout << "-- Press ENTER to exit...";
        std::cin.get();
        exit(1);
    } else {
        // Initialize status variable to evaluate the existence of files
        int statusOfAccess;
        // If we got the exact number of program arguments, load them
        for (int i = 1; i < argc; i++) {
            /* We iterate over argv[] to get the parameters stored inside.
            Note that we're starting on 1 because we don't need to know
            the program name, which is stored in argv[0] */
            if (std::string(argv[i]) == "-i") {
                /// NEW ANALYSIS
                // We know the next argument *should* be the input filename
                argFileName = std::string(argv[i + 1]);
                // Check for existence and readability of the input file
                statusOfAccess = access(argFileName.c_str(), F_OK | R_OK);
                if (statusOfAccess != 0) {
                    std::cerr << "Error: impossible to access the input file: "
                              << argFileName << std::endl;
                    std::cerr << strerror(errno) << std::endl;
                    exit(1);
                } else {
                    // All good with new analysis, so check restart set to false!
                    checkRestart = false;
                    i++;
                }
            } else if (std::string(argv[i]) == "-r") {
                /// RESTART ANALYSIS
                // We know the next argument *should* be the restart filename:
                argFileName = std::string(argv[i + 1]);
                // Check for existence and readability of the input file
                statusOfAccess = access(argFileName.c_str(), F_OK | R_OK);
                if (statusOfAccess != 0) {
                    std::cerr << "Error: impossible to access the restart file: "
                              << argFileName << std::endl;
                    std::cerr << strerror(errno) << std::endl;
                    exit(1);
                } else {
                    // All good with restart analysis
                    checkRestart = true;
                    i++;
                }
            } else {
                /// ERROR
                // There is an invalid argument, so send an error.
                std::cerr << "Invalid argument " << argv[i] << std::endl;
                std::cerr << "-- Press ENTER to exit...";
                std::cin.get();
                exit(1);
            }
        }
        // Now that we are sure the program arguments are passed correctly,
        // extract json object from input or restart file
        std::ifstream inputRestart(argFileName);
        inputRestart >> js;
    }
}

} // namespace EQSim
