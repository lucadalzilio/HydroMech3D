//
// Created by Federico Ciardo on 11.10.22. All rights reserved.
//
#include <armadillo>
#include <src/Mesh.h>
#include "FullSpaceElasticity.h"

#ifndef INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H
#define INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H

namespace EQSim {

arma::mat AssembleElastMat(Mesh &mesh, EQSim::SolidMatrixProperties &matrix_prop);

arma::mat RotationMatrix3D(arma::vec &normal_vector,
                                     double &theta);

}
#endif  // INC_3DEQSIM_SRC_ASSEMBLEELASTICITYMATRIX_H
