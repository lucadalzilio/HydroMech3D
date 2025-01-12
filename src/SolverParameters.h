//
// Created by Federico Ciardo on 03.10.21. All rights reserved.
//

#ifndef INC_3DEQSIM_SRC_SOLVERPARAMETERS_H
#define INC_3DEQSIM_SRC_SOLVERPARAMETERS_H

namespace EQSim {

struct SolverParameters {
  double initial_time;
  double time_Step;
  double maximum_time;
  double tolerance_RK = 10e-6;  // default value
  double time_step_amplification_factor;
  double time_step_reduction_factor;
  int export_current_solution_every_i_time_steps_;
};

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_SOLVERPARAMETERS_H
