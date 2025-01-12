//
// Created by Federico Ciardo on 10.08.21. All rights reserved.
//

// Import from the project
#include "Mesh.h"

#ifndef INC_3DEQSIM_SRC_INSITUSTRESS_H
#define INC_3DEQSIM_SRC_INSITUSTRESS_H

namespace EQSim {

class FaultInSituStress {
 private:
  arma::vec sxx_;  // sigma_xx
  arma::vec syy_;  // sigma_yy
  arma::vec szz_;
  arma::vec sxy_;
  arma::vec sxz_;
  arma::vec syz_;

 public:
  // Constructors
  FaultInSituStress() = default;

  FaultInSituStress(arma::vec &sigmaxx, arma::vec &sigma_yy,
                    arma::vec &sigma_zz, arma::vec &sigma_xy,
                    arma::vec &sigma_xz, arma::vec &sigma_yz) {
    sxx_ = sigmaxx;
    syy_ = sigma_yy;
    szz_ = sigma_zz;
    sxy_ = sigma_xy;
    sxz_ = sigma_xz;
    syz_ = sigma_yz;
  };

  // Getter methods
  arma::vec getSxx() const { return sxx_; };
  arma::vec getSyy() const { return syy_; };
  arma::vec getSzz() const { return szz_; };
  arma::vec getSxy() const { return sxy_; };
  arma::vec getSxz() const { return sxz_; };
  arma::vec getSyz() const { return syz_; };
  double getSxx(arma::uword &i) const { return sxx_[i]; };
  double getSyy(arma::uword &i) const { return syy_[i]; };
  double getSzz(arma::uword &i) const { return szz_[i]; };
  double getSxy(arma::uword &i) const { return sxy_[i]; };
  double getSxz(arma::uword &i) const { return sxz_[i]; };
  double getSyz(arma::uword &i) const { return syz_[i]; };

  // Methods
  arma::vec InsituTractionsElt(EQSim::ElementData &my_elt,
                                       arma::uword &elt_idx) const {
    arma::vec n = my_elt.getN();
    arma::vec s1 = my_elt.getS1();
    arma::vec s2 = my_elt.getS2();

    arma::vec traction(3, arma::fill::zeros);

    // s1.Sig.n
    traction[0] = (s1[0] * sxx_[elt_idx] * n[0] + s1[0] * sxy_[elt_idx] * n[1] +
                   s1[0] * sxz_[elt_idx] * n[2]) +
                  (s1[1] * sxy_[elt_idx] * n[0] + s1[1] * syy_[elt_idx] * n[1] +
                   s1[1] * syz_[elt_idx] * n[2]) +
                  (s1[2] * sxz_[elt_idx] * n[0] + s1[2] * syz_[elt_idx] * n[1] +
                   s1[2] * szz_[elt_idx] * n[2]);  // ts1

    // s2.Sig.n
    traction[1] = (s2[0] * sxx_[elt_idx] * n[0] + s2[0] * sxy_[elt_idx] * n[1] +
                   s2[0] * sxz_[elt_idx] * n[2]) +
                  (s2[1] * sxy_[elt_idx] * n[0] + s2[1] * syy_[elt_idx] * n[1] +
                   s2[1] * syz_[elt_idx] * n[2]) +
                  (s2[2] * sxz_[elt_idx] * n[0] + s2[2] * syz_[elt_idx] * n[1] +
                   s2[2] * szz_[elt_idx] * n[2]);  // ts2

    // n.Sig.n
    traction[2] = (n[0] * sxx_[elt_idx] * n[0] + n[0] * sxy_[elt_idx] * n[1] +
                   n[0] * sxz_[elt_idx] * n[2]) +
                  (n[1] * sxy_[elt_idx] * n[0] + n[1] * syy_[elt_idx] * n[1] +
                   n[1] * syz_[elt_idx] * n[2]) +
                  (n[2] * sxz_[elt_idx] * n[0] + n[2] * syz_[elt_idx] * n[1] +
                   n[2] * szz_[elt_idx] * n[2]);  // tn

    return traction;
  };

  arma::vec AllInSituTractions(EQSim::Mesh &my_mesh) const {
    arma::vec inSituTractionsEltE;
    arma::vec all_tractions(my_mesh.getNumberOfDofs(), arma::fill::zeros);

    // Number of Dofs per element
    arma::uword NDofsPerElt = 3;
    // Number of total elements
    arma::uword Nelts = my_mesh.getNumberOfElts();

    // Loop over all the elements
    for (arma::uword e = 0; e < Nelts; e++) {
      // Get data of element e
      EQSim::ElementData my_elt = my_mesh.getElementData(e);

      // Get tractions on element e
      inSituTractionsEltE = InsituTractionsElt(my_elt, e);

      // Fill the vector
      all_tractions[e * NDofsPerElt] = inSituTractionsEltE[0];
      all_tractions[e * NDofsPerElt + 1] = inSituTractionsEltE[1];
      all_tractions[e * NDofsPerElt + 2] = inSituTractionsEltE[2];
    }

    return all_tractions;
  }
};
}  // namespace EQSim

#endif  // INC_3DEQSIM_SRC_INSITUSTRESS_H
