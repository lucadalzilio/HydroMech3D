//
// Created by Federico Ciardo on 11.08.21. All rights reserved.
//

// Inclusion from IL library
#include <armadillo>
// Inclusion from the project
#include "SolidMatrixProperties.h"
#include "src/StressKernelsP0/StressKernelsDxP0.h"

#ifndef INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H
#define INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H

namespace EQSim {

arma::mat TractionsDueToDDsOnSingleEltP0(
    double &a, double &b, arma::vec &coor_centroid_elt,
    arma::vec &shear1_vector, arma::vec &shear2_vector,
    arma::vec &normal_vector, SolidMatrixProperties Matrix_Prop);

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_FULLSPACEELASTICITY_H
