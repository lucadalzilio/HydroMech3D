//
// Created by Federico Ciardo on 17.08.21. All rights reserved.
//

#include "StressKernelsDzP0.h"

namespace EQSim {

arma::mat StressTensorDueToDDzOnSingleEltP0(double &a, double &b,
                                                      double &x1, double &x2,
                                                      double &x3, double &Nu,
                                                      double &G) {
  double Sigma_xx_Dz, Sigma_yy_Dz, Sigma_zz_Dz, Sigma_xy_Dz, Sigma_xz_Dz,
      Sigma_yz_Dz;

  arma::mat StressTensorDueToDz(3, 3, arma::fill::zeros);

  Sigma_xx_Dz =
      (G *
       (((-a + x1) * (-b + x2) *
         (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        ((a + x1) * (-b + x2) *
         (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        ((-a + x1) * (b + x2) *
         (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) +
        ((a + x1) * (b + x2) *
         (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        x3 * (((-a + x1) * (-b + x2) * x3 *
               (3 * pow(a, 2) + 2 * pow(b, 2) - 6 * a * x1 + 3 * pow(x1, 2) -
                4 * b * x2 + 2 * pow(x2, 2) + 3 * pow(x3, 2))) /
                  (pow(pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) -
              ((a + x1) * (-b + x2) * x3 *
               (3 * pow(a, 2) + 2 * pow(b, 2) + 6 * a * x1 + 3 * pow(x1, 2) -
                4 * b * x2 + 2 * pow(x2, 2) + 3 * pow(x3, 2))) /
                  (pow(pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) -
              ((-a + x1) * (b + x2) * x3 *
               (3 * pow(a, 2) + 2 * pow(b, 2) - 6 * a * x1 + 3 * pow(x1, 2) +
                4 * b * x2 + 2 * pow(x2, 2) + 3 * pow(x3, 2))) /
                  (pow(pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) +
              ((a + x1) * (b + x2) * x3 *
               (3 * pow(a, 2) + 2 * pow(b, 2) + 6 * a * x1 + 3 * pow(x1, 2) +
                4 * b * x2 + 2 * pow(x2, 2) + 3 * pow(x3, 2))) /
                  (pow(pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5))) +
        (1 - 2 * Nu) *
            ((-b + x2) /
                 (sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                  (-a + x1 +
                   sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             (-b + x2) /
                 (sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                  (a + x1 +
                   sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             (b + x2) /
                 (sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                  (-a + x1 +
                   sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) +
             (b + x2) /
                 (sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                  (a + x1 +
                   sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_yy_Dz =
      (G *
       (((-a + x1) * (-b + x2) *
         (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        ((a + x1) * (-b + x2) *
         (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        ((-a + x1) * (b + x2) *
         (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) +
        ((a + x1) * (b + x2) *
         (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        x3 * (((-a + x1) * (-b + x2) * x3 *
               (2 * pow(a, 2) + 3 * pow(b, 2) - 4 * a * x1 + 2 * pow(x1, 2) -
                6 * b * x2 + 3 * pow(x2, 2) + 3 * pow(x3, 2))) /
                  (pow(pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) -
              ((a + x1) * (-b + x2) * x3 *
               (2 * pow(a, 2) + 3 * pow(b, 2) + 4 * a * x1 + 2 * pow(x1, 2) -
                6 * b * x2 + 3 * pow(x2, 2) + 3 * pow(x3, 2))) /
                  (pow(pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) -
              ((-a + x1) * (b + x2) * x3 *
               (2 * pow(a, 2) + 3 * pow(b, 2) - 4 * a * x1 + 2 * pow(x1, 2) +
                6 * b * x2 + 3 * pow(x2, 2) + 3 * pow(x3, 2))) /
                  (pow(pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) +
              ((a + x1) * (b + x2) * x3 *
               (2 * pow(a, 2) + 3 * pow(b, 2) + 4 * a * x1 + 2 * pow(x1, 2) +
                6 * b * x2 + 3 * pow(x2, 2) + 3 * pow(x3, 2))) /
                  (pow(pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5))) +
        (1 - 2 * Nu) *
            ((-a + x1) /
                 (sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                  (-b + x2 +
                   sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             (a + x1) /
                 (sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                  (-b + x2 +
                   sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             (-a + x1) /
                 (sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                  (b + x2 +
                   sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) +
             (a + x1) /
                 (sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                  (b + x2 +
                   sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_zz_Dz =
      (G *
       (((-a + x1) * (-b + x2) *
         (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        ((a + x1) * (-b + x2) *
         (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        ((-a + x1) * (b + x2) *
         (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) +
        ((a + x1) * (b + x2) *
         (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
          pow(x2, 2) + 2 * pow(x3, 2))) /
            ((pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
             (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
             sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                  pow(x2, 2) + pow(x3, 2))) -
        x3 * (((-a + x1) * (-b + x2) * x3 *
               (4 * (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) -
                (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)) -
                2 * (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)) -
                2 * (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)))) /
                  (pow(pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
                   pow(pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) -
              ((a + x1) * (-b + x2) * x3 *
               (4 * (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) -
                (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)) -
                2 * (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)) -
                2 * (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)))) /
                  (pow(pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
                   pow(pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) -
              ((-a + x1) * (b + x2) * x3 *
               (4 * (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) -
                (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)) -
                2 * (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)) -
                2 * (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)))) /
                  (pow(pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
                   pow(pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5)) +
              ((a + x1) * (b + x2) * x3 *
               (4 * (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) -
                (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)) -
                2 * (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)) -
                2 * (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
                    (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                     2 * b * x2 + pow(x2, 2) + 2 * pow(x3, 2)))) /
                  (pow(pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
                   pow(pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
                   pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                       1.5))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_xy_Dz =
      (G * (-(x3 * (-(x3 / pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                                   2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                               1.5)) +
                    x3 / pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                                 2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                             1.5) +
                    x3 / pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                                 2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                             1.5) -
                    x3 / pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                                 2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                             1.5))) +
            (-1 + 2 * Nu) *
                (1 / sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2)) -
                 1 / sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2)) -
                 1 / sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2)) +
                 1 / sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_xz_Dz =
      -0.25 *
      (G * x3 *
       (((b - x2) *
         (-(pow(x3, 2) * (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2))) -
          2 * pow(x3, 2) *
              (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)) +
          (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
              (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)))) /
            (pow(pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
             pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                     pow(x2, 2) + pow(x3, 2),
                 1.5)) -
        ((b - x2) *
         (-(pow(x3, 2) * (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2))) -
          2 * pow(x3, 2) *
              (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)) +
          (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
              (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)))) /
            (pow(pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
             pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                     pow(x2, 2) + pow(x3, 2),
                 1.5)) -
        ((-b - x2) *
         (-(pow(x3, 2) * (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2))) -
          2 * pow(x3, 2) *
              (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)) +
          (pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
              (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)))) /
            (pow(pow(a, 2) - 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
             pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                     pow(x2, 2) + pow(x3, 2),
                 1.5)) +
        ((-b - x2) *
         (-(pow(x3, 2) * (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2))) -
          2 * pow(x3, 2) *
              (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)) +
          (pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2)) *
              (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)))) /
            (pow(pow(a, 2) + 2 * a * x1 + pow(x1, 2) + pow(x3, 2), 2) *
             pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                     pow(x2, 2) + pow(x3, 2),
                 1.5)))) /
      ((1 - Nu) * arma::datum::pi);

  Sigma_yz_Dz =
      -0.25 *
      (G * x3 *
       (((a - x1) *
         (-(pow(x3, 2) * (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2))) -
          2 * pow(x3, 2) *
              (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)) +
          (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
              (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)))) /
            (pow(pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
             pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                     pow(x2, 2) + pow(x3, 2),
                 1.5)) -
        ((-a - x1) *
         (-(pow(x3, 2) * (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2))) -
          2 * pow(x3, 2) *
              (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)) +
          (pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
              (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)))) /
            (pow(pow(b, 2) - 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
             pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) - 2 * b * x2 +
                     pow(x2, 2) + pow(x3, 2),
                 1.5)) -
        ((a - x1) *
         (-(pow(x3, 2) * (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2))) -
          2 * pow(x3, 2) *
              (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)) +
          (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
              (pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)))) /
            (pow(pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
             pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                     pow(x2, 2) + pow(x3, 2),
                 1.5)) +
        ((-a - x1) *
         (-(pow(x3, 2) * (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2))) -
          2 * pow(x3, 2) *
              (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)) +
          (pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2)) *
              (pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
               pow(x2, 2) + pow(x3, 2)))) /
            (pow(pow(b, 2) + 2 * b * x2 + pow(x2, 2) + pow(x3, 2), 2) *
             pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) + 2 * b * x2 +
                     pow(x2, 2) + pow(x3, 2),
                 1.5)))) /
      ((1 - Nu) * arma::datum::pi);

  StressTensorDueToDz(0, 0) = Sigma_xx_Dz;
  StressTensorDueToDz(0, 1) = Sigma_xy_Dz;
  StressTensorDueToDz(0, 2) = Sigma_xz_Dz;
  StressTensorDueToDz(1, 0) = Sigma_xy_Dz;
  StressTensorDueToDz(1, 1) = Sigma_yy_Dz;
  StressTensorDueToDz(1, 2) = Sigma_yz_Dz;
  StressTensorDueToDz(2, 0) = Sigma_xz_Dz;
  StressTensorDueToDz(2, 1) = Sigma_yz_Dz;
  StressTensorDueToDz(2, 2) = Sigma_zz_Dz;

  return StressTensorDueToDz;
}

}
