//
// Created by Federico Ciardo on 11.08.21. All rights reserved.
//

#include "Utils.h"

namespace EQSim {

inline arma::mat RotationMatrix3D(arma::vec &normal_vector,
                                            double &theta) {
  arma::mat R(3, 3, arma::fill::zeros);

  double n1, n2, n3, n1square, n2square, n3square;
  n1 = normal_vector[0];
  n2 = normal_vector[1];
  n3 = normal_vector[2];
  n1square = n1 * n1;
  n2square = n2 * n2;
  n3square = n3 * n3;

  R(0, 0) = cos(theta) + (n1square * (1 - cos(theta)));
  R(0, 1) = (n1 * n2 * (1 - cos(theta))) - (n3 * sin(theta));
  R(0, 2) = (n1 * n3 * (1 - cos(theta))) + (n2 * sin(theta));

  R(1, 0) = (n1 * n2 * (1 - cos(theta))) + n3 * sin(theta);
  R(1, 1) = cos(theta) + (n2square * (1 - cos(theta)));
  R(1, 2) = (n2 * n3 * (1 - cos(theta))) - (n1 * sin(theta));

  R(2, 0) = (n1 * n3 * (1 - cos(theta))) - (n2 * sin(theta));
  R(2, 1) = (n2 * n3 * (1 - cos(theta))) + (n1 * sin(theta));
  R(2, 2) = cos(theta) + (n3square * (1 - cos(theta)));

  return R;
}

/////////
// This function calculates the 3D euclidean distance between two points;
// Point 1=(x1 = x1, x2 = y1, x3 = z1) -- Point 2=(y1 = x2, y2 = y2, y3 = z2)
inline double euclideanDistance(double x1, double x2, double x3, double y1,
                                double y2, double y3) {
  double res = 0.;

  res = sqrt(pow(x1 - y1, 2) + pow(x2 - y2, 2) + pow(x3 - y3, 2));

  return res;
}

/////////
template <typename T>
inline arma::Mat<T> position_2d_array(const arma::Mat<T> &arr2D, T seek) {
  arma::uword arr2D_size1 = arr2D.n_cols;
  arma::uword arr2D_size0 = arr2D.n_rows;
  arma::Mat<T> M(arr2D_size1 * arr2D_size0, 2, arma::fill::value(-1));
  arma::uword k = 0;

  for (arma::uword i = 0; i < arr2D_size0; ++i) {
    for (arma::uword j = 0; j < arr2D_size1; ++j) {
      if (arr2D(i, j) == seek) {
        M(k, 0) = i;
        M(k, 1) = j;
        k = k + 1;
      }
    }
  }

  arma::imat outp(k, 2, arma::fill::zeros);

  for (arma::uword l = 0; l < k; ++l) {
    for (arma::uword k2 = 0; k2 < 2; ++k2) {
      outp(l, k2) = M(l, k2);
    }
  }

  return outp;
}

/////////
template <class T>
inline arma::Col<T> deleteDuplicates(const arma::Col<T> &arr) {
  arma::Col<T> res;

  for (arma::uword i = 0; i < arr.n_elem; ++i) {
    bool already_there = false;
    for (arma::uword j = 0; j < res.n_elem; ++j) {
      if (arr[i] == res[j]) {
        already_there = true;
      }
    }
    if (!already_there) {
    	res.insert_rows(res.n_elem, 1);
	res(res.n_elem - 1) = arr[i];
    }
  }

  return res;
}

////////
inline arma::vec CalculatePressureLaplacianUniformMeshOnly(
    arma::imat &neigh_elts, EQSim::Mesh &Mesh,
    const arma::vec &my_vector) {
  arma::vec res(my_vector.n_elem, arma::fill::zeros);
  arma::vec resX(my_vector.n_elem, arma::fill::zeros);
  arma::vec resY(my_vector.n_elem, arma::fill::zeros);
  arma::ivec neighElts;
  arma::uword t;
  EQSim::ElementData elemI;
  arma::uword N_eltsX, N_eltsY;

  for (arma::uword I = 0; I < neigh_elts.n_rows; ++I) {
    elemI = Mesh.getElementData(I);

    // Find the real neighbour elements, i.e. the ones different than -1.
    t = 0;
    for (arma::uword J = 0; J < neigh_elts.n_cols; ++J) {
      if (neigh_elts(I, J) != -1) {
        neighElts.resize(t + 1);
        neighElts[t] = neigh_elts(I, J);
        ++t;
      }
    }

    N_eltsX = 0;
    N_eltsY = 0;
    // Loop over the real neighbour elements
    for (arma::uword J = 0; J < neighElts.n_elem; ++J) {
      if ((Mesh.getCentroids(neighElts[J], 1) == Mesh.getCentroids(I, 1)) &&
          (neighElts[J] != static_cast<long long int>(I))) {
        resX[I] += my_vector[neighElts[J]];
        ++N_eltsX;
      }

      if ((Mesh.getCentroids(neighElts[J], 0) == Mesh.getCentroids(I, 0)) &&
          (neighElts[J] != static_cast<long long int>(I))) {
        resY[I] += my_vector[neighElts[J]];
        ++N_eltsY;
      }
    }

    if (N_eltsX == 2) {
      resX[I] -= 2. * my_vector[I];
      resX[I] = (resX[I] / (pow((2. * elemI.getA()), 2)));
    } else {
      resX[I] -= my_vector[I];
      resX[I] = (resX[I] / (pow((2. * elemI.getA()), 2)));
    }

    if (N_eltsY == 2) {
      resY[I] -= 2. * my_vector[I];
      resY[I] = (resY[I] / (pow((2. * elemI.getB()), 2)));
    } else {
      resY[I] -= my_vector[I];
      resY[I] = (resY[I] / (pow((2. * elemI.getB()), 2)));
    }

    res[I] = resX[I] + resY[I];
  }

  return res;
}

}  // namespace EQSim
