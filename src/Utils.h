//
// Created by Federico Ciardo on 11.08.21. All rights reserved.
//

// Inclusion from IL library
#include <armadillo>

// Inclusion from the project
#include "ElementData.h"
#include "Mesh.h"

#ifndef INC_3DEQSIM_SRC_UTILS_H
#define INC_3DEQSIM_SRC_UTILS_H

namespace EQSim {

inline arma::mat RotationMatrix3D(arma::vec &normal_vector,
                                            double &theta);

inline double euclideanDistance(double x1, double x2, double x3, double y1,
                                double y2, double y3);

template <typename T>
inline arma::Mat<T> position_2d_array(const arma::Mat<T> &arr2D, T seek);

template <class T>
inline arma::Col<T> deleteDuplicates(const arma::Col<T> &arr);

arma::vec CalculatePressureLaplacianUniformMeshOnly(
    arma::imat &neigh_elts, EQSim::Mesh &Mesh,
    const arma::vec &my_vector);

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_UTILS_H
