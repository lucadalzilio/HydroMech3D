//
// Created by Federico Ciardo on 30.08.21. All rights reserved.
//

// Import from the project
#include "RightHandSidesODEs.h"

#ifndef INC_3DEQSIM_SRC_RK45_H
#define INC_3DEQSIM_SRC_RK45_H

namespace EQSim {
class RK45 {
 private:
  RightHandSideODEs* rhsODEs_;
  SolutionRK45 Solution_;

 public:
  RK45(RightHandSideODEs* rhsODEs, SolutionRK45 &SolutionRK45){
    rhsODEs_ = rhsODEs;
    Solution_ = SolutionRK45;
  }

  void Solve();
};
}

#endif  // INC_3DEQSIM_SRC_RK45_H
