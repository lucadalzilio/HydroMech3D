//
// Created by Federico Ciardo on 10.08.21. All rights reserved.
//

// Inclusion from armadillo library
#include <armadillo>

#ifndef INC_3DEQSIM_SRC_ELEMENTDATA_H
#define INC_3DEQSIM_SRC_ELEMENTDATA_H

namespace EQSim {

class ElementData {
 private:
  double theta_;                // Angle of s1 with respect to e1 = (1,0,0)
  double a_;                    // Half-length of element on one direction
  double b_;                    // Half-length of element on the other direction
  double areaElt_;     // Area of the element
  arma::vec n_;  // Normal vector
  arma::vec s1_;  // Shear vector 1
  arma::vec s2_;  // Shear vector 2
  arma::vec centroid_;
  
 public:
  ElementData() = default;
  ElementData(const arma::mat &CoorElt, const arma::vec &CentroidElt) {
    if (CoorElt.n_rows < 3 || CoorElt.n_cols < 3) {
      throw std::invalid_argument("CoorElt must be at least 3x3");
    }
    if (CentroidElt.n_elem < 3) {
      throw std::invalid_argument("CentroidElt must have at least 3 elements");
    }
    a_ = (std::sqrt(pow(fabs(CoorElt(1, 0) - CoorElt(2, 0)), 2) +
                    pow(fabs(CoorElt(1, 1) - CoorElt(2, 1)), 2) +
                    pow(fabs(CoorElt(1, 2) - CoorElt(2, 2)), 2))) /
         2.;

    b_ = (std::sqrt(pow(fabs(CoorElt(0, 0) - CoorElt(1, 0)), 2) +
                    pow(fabs(CoorElt(0, 1) - CoorElt(1, 1)), 2) +
                    pow(fabs(CoorElt(0, 2) - CoorElt(1, 2)), 2))) /
         2.;

    // The sign of the othonormal vectors follow the Right Hand rule
    arma::mat DiffCoorElt(3, 3, arma::fill::zeros);
    for (arma::uword I = 0; I < DiffCoorElt.n_rows; ++I) {
      DiffCoorElt(I, 0) = -1. * CoorElt(I, 0) + CoorElt(I + 1, 0);
      DiffCoorElt(I, 1) = -1. * CoorElt(I, 1) + CoorElt(I + 1, 1);
      DiffCoorElt(I, 2) = -1. * CoorElt(I, 2) + CoorElt(I + 1, 2);
    }

    arma::vec aux(3, arma::fill::zeros);
    double normaux;
    aux[0] = (-1. * DiffCoorElt(0, 2) * DiffCoorElt(1, 1)) +
             DiffCoorElt(0, 1) * DiffCoorElt(1, 2);
    aux[1] = DiffCoorElt(0, 2) * DiffCoorElt(1, 0) -
             (DiffCoorElt(0, 0) * DiffCoorElt(1, 2));
    aux[2] = (-1. * DiffCoorElt(0, 1) * DiffCoorElt(1, 0)) +
             DiffCoorElt(0, 0) * DiffCoorElt(1, 1);
    normaux = arma::norm(aux, 1);

    n_ = aux / normaux;

    arma::vec aux1(3, arma::fill::zeros);
    double normaux1;
    aux1[0] = (-1. * (-1. * n_[2]) * DiffCoorElt(0, 1)) +
              (-1. * n_[1]) * DiffCoorElt(0, 2);
    aux1[1] =
        (-1. * n_[2]) * DiffCoorElt(0, 0) - ((-1. * n_[0]) * DiffCoorElt(0, 2));
    aux1[2] = (-1. * (-1. * n_[1]) * DiffCoorElt(0, 0)) +
              (-1. * n_[0]) * DiffCoorElt(0, 1);

    normaux1 = arma::norm(aux1, 1);

    s1_ = aux1 / normaux1;

    arma::vec aux2(3, arma::fill::zeros);
    double normaux2;
    aux2[0] = (-1. * (-1. * n_[2]) * DiffCoorElt(1, 1)) +
              (-1. * n_[1]) * DiffCoorElt(1, 2);
    aux2[1] =
        (-1. * n_[2]) * DiffCoorElt(1, 0) - ((-1. * n_[0]) * DiffCoorElt(1, 2));
    aux2[2] = (-1. * (-1. * n_[1]) * DiffCoorElt(1, 0)) +
              (-1. * n_[0]) * DiffCoorElt(1, 1);
    normaux2 = arma::norm(aux2, 1);

    s2_ = aux2 / normaux2;

    arma::vec e1(3, arma::fill::zeros);
    e1[0] = 1.;
    double costheta = arma::dot(e1, s1_) /
                      (arma::norm(e1, 1) * arma::norm(s1_, 1));
    theta_ = std::acos(costheta);
    if (s1_[1]<0){ theta_ = -1.* theta_; }

    centroid_ = CentroidElt;
    areaElt_ = ((2*a_) * (2*b_));
  }

  // Getter methods
  arma::vec getN() { return n_; };
  arma::vec getS1() { return s1_; };
  arma::vec getS2() { return s2_; };
  double getTheta() const { return theta_; };
  double getAreaElt() const { return areaElt_; };
  arma::vec getCentroidElt() const { return centroid_; };
  double getCentroidElt(arma::uword i) const { return centroid_[i]; };
  double getA() const { return a_; };
  double getB() const { return b_; };

};

}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_ELEMENTDATA_H
