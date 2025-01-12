//
// Created by Federico Ciardo on 17.08.21. All rights reserved.
//

#include "StressKernelsDyP0.h"

namespace EQSim {

arma::mat StressTensorDueToDDyOnSingleEltP0(double &a, double &b,
                                                      double &x1, double &x2,
                                                      double &x3, double &Nu,
                                                      double &G) {
  double Sigma_xx_Dy, Sigma_yy_Dy, Sigma_zz_Dy, Sigma_xy_Dy, Sigma_xz_Dy,
      Sigma_yz_Dy;

  arma::mat StressTensorDueToDy(3, 3, arma::fill::zeros);

  Sigma_xx_Dy =
      (G *
       (-(x3 *
          ((a - x1) / pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                              2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                          1.5) -
           (-a - x1) / pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                           1.5) -
           (a - x1) / pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                              2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                          1.5) +
           (-a - x1) / pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                           1.5))) +
        2 * Nu *
            (x3 / (sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (-a + x1 +
                    sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             x3 / (sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (a + x1 +
                    sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             x3 / (sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (-a + x1 +
                    sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) +
             x3 / (sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (a + x1 +
                    sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_yy_Dy =
      (G *
       (2 * (x3 / (sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (-a + x1 +
                    sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             x3 / (sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (a + x1 +
                    sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             x3 / (sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (-a + x1 +
                    sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) +
             x3 / (sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (a + x1 +
                    sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))))) -
        x3 * ((-(pow(-b + x2, 2) *
                 sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) -
               pow(-b + x2, 2) *
                   (-a + x1 +
                    sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) +
               (pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (-a + x1 +
                    sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) /
                  (pow(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2), 1.5) *
                   pow(-a + x1 +
                           sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)),
                       2)) -
              (-(pow(-b + x2, 2) *
                 sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) -
               pow(-b + x2, 2) *
                   (a + x1 +
                    sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) +
               (pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (a + x1 +
                    sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) /
                  (pow(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2), 1.5) *
                   pow(a + x1 +
                           sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)),
                       2)) -
              (-(pow(b + x2, 2) *
                 sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) -
               pow(b + x2, 2) *
                   (-a + x1 +
                    sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) +
               (pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (-a + x1 +
                    sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) /
                  (pow(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2), 1.5) *
                   pow(-a + x1 +
                           sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)),
                       2)) +
              (-(pow(b + x2, 2) *
                 sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) -
               pow(b + x2, 2) *
                   (a + x1 +
                    sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) +
               (pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (a + x1 +
                    sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) /
                  (pow(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2), 1.5) *
                   pow(a + x1 +
                           sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)),
                       2))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_zz_Dy =
      -0.25 *
      (G * x3 *
       ((-(pow(x3, 2) * sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) -
         pow(x3, 2) *
             (-a + x1 + sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) +
         (pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
             (-a + x1 + sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) /
            (pow(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2), 1.5) *
             pow(-a + x1 + sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)),
                 2)) -
        (-(pow(x3, 2) * sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) -
         pow(x3, 2) *
             (a + x1 + sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) +
         (pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
             (a + x1 + sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) /
            (pow(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2), 1.5) *
             pow(a + x1 + sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)),
                 2)) -
        (-(pow(x3, 2) * sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) -
         pow(x3, 2) *
             (-a + x1 + sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) +
         (pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
             (-a + x1 + sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) /
            (pow(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2), 1.5) *
             pow(-a + x1 + sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)),
                 2)) +
        (-(pow(x3, 2) * sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) -
         pow(x3, 2) *
             (a + x1 + sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) +
         (pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
             (a + x1 + sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) /
            (pow(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2), 1.5) *
             pow(a + x1 + sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)),
                 2)))) /
      ((1 - Nu) * arma::datum::pi);

  Sigma_xy_Dy =
      (G *
       (-(x3 *
          ((b - x2) / pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                              2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                          1.5) -
           (b - x2) / pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                              2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                          1.5) -
           (-b - x2) / pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                           1.5) +
           (-b - x2) / pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                           1.5))) +
        (1 - Nu) *
            (x3 / (sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (-b + x2 +
                    sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             x3 / (sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (-b + x2 +
                    sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) -
             x3 / (sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (b + x2 +
                    sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) +
             x3 / (sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (b + x2 +
                    sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_xz_Dy =
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
                             1.5))) -
            Nu * (1 / sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2)) -
                  1 / sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2)) -
                  1 / sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2)) +
                  1 / sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                           2 * b * x2 + pow(x2, 2) + pow(x3, 2))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_yz_Dy =
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
        x3 *
            (-(((-b + x2) * x3 *
                (-a + x1 +
                 2 * sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2)))) /
               (pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                        2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                    1.5) *
                pow(-a + x1 +
                        sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                             2 * b * x2 + pow(x2, 2) + pow(x3, 2)),
                    2))) +
             ((-b + x2) * x3 *
              (a + x1 +
               2 * sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                        2 * b * x2 + pow(x2, 2) + pow(x3, 2)))) /
                 (pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                      1.5) *
                  pow(a + x1 +
                          sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2)),
                      2)) +
             ((b + x2) * x3 *
              (-a + x1 +
               2 * sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                        2 * b * x2 + pow(x2, 2) + pow(x3, 2)))) /
                 (pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                      1.5) *
                  pow(-a + x1 +
                          sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2)),
                      2)) -
             ((b + x2) * x3 *
              (a + x1 +
               2 * sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                        2 * b * x2 + pow(x2, 2) + pow(x3, 2)))) /
                 (pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                      1.5) *
                  pow(a + x1 +
                          sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2)),
                      2))) +
        Nu * ((-a + x1) /
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

  StressTensorDueToDy(0, 0) = Sigma_xx_Dy;
  StressTensorDueToDy(0, 1) = Sigma_xy_Dy;
  StressTensorDueToDy(0, 2) = Sigma_xz_Dy;
  StressTensorDueToDy(1, 0) = Sigma_xy_Dy;
  StressTensorDueToDy(1, 1) = Sigma_yy_Dy;
  StressTensorDueToDy(1, 2) = Sigma_yz_Dy;
  StressTensorDueToDy(2, 0) = Sigma_xz_Dy;
  StressTensorDueToDy(2, 1) = Sigma_yz_Dy;
  StressTensorDueToDy(2, 2) = Sigma_zz_Dy;

  return StressTensorDueToDy;
}

}
