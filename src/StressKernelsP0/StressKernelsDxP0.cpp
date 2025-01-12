//
// Created by Federico Ciardo on 16.08.21. All rights reserved.
//

#include "StressKernelsDxP0.h"

namespace EQSim {

arma::mat StressTensorDueToDDxOnSingleEltP0(double &a, double &b,
                                                      double &x1, double &x2,
                                                      double &x3, double &Nu,
                                                      double &G) {
  double Sigma_xx_Dx, Sigma_yy_Dx, Sigma_zz_Dx, Sigma_xy_Dx, Sigma_xz_Dx,
      Sigma_yz_Dx;

  arma::mat StressTensorDueToDx(3, 3, arma::fill::zeros);

  Sigma_xx_Dx =
      (G *
       (2 * (x3 / (sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
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
                    sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))))) -
        x3 * ((-(pow(-a + x1, 2) *
                 sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) -
               pow(-a + x1, 2) *
                   (-b + x2 +
                    sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) +
               (pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (-b + x2 +
                    sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) /
                  (pow(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2), 1.5) *
                   pow(-b + x2 +
                           sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)),
                       2)) -
              (-(pow(a + x1, 2) *
                 sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) -
               pow(a + x1, 2) *
                   (-b + x2 +
                    sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) +
               (pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
                   (-b + x2 +
                    sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) /
                  (pow(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2), 1.5) *
                   pow(-b + x2 +
                           sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)),
                       2)) -
              (-(pow(-a + x1, 2) *
                 sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) -
               pow(-a + x1, 2) *
                   (b + x2 +
                    sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) +
               (pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (b + x2 +
                    sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) /
                  (pow(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2), 1.5) *
                   pow(b + x2 +
                           sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)),
                       2)) +
              (-(pow(a + x1, 2) *
                 sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) -
               pow(a + x1, 2) *
                   (b + x2 +
                    sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) +
               (pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
                   (b + x2 +
                    sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) /
                  (pow(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2), 1.5) *
                   pow(b + x2 +
                           sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)),
                       2))))) /
      (4. * (1 - Nu) * arma::datum::pi);

  Sigma_yy_Dx =
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
        2 * Nu *
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

  Sigma_zz_Dx =
      -0.25 *
      (G * x3 *
       ((-(pow(x3, 2) * sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) -
         pow(x3, 2) *
             (-b + x2 + sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) +
         (pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
             (-b + x2 + sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) /
            (pow(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2), 1.5) *
             pow(-b + x2 + sqrt(pow(-a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)),
                 2)) -
        (-(pow(x3, 2) * sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) -
         pow(x3, 2) *
             (-b + x2 + sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2))) +
         (pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)) *
             (-b + x2 + sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)))) /
            (pow(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2), 1.5) *
             pow(-b + x2 + sqrt(pow(a + x1, 2) + pow(-b + x2, 2) + pow(x3, 2)),
                 2)) -
        (-(pow(x3, 2) * sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) -
         pow(x3, 2) *
             (b + x2 + sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) +
         (pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
             (b + x2 + sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) /
            (pow(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2), 1.5) *
             pow(b + x2 + sqrt(pow(-a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)),
                 2)) +
        (-(pow(x3, 2) * sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) -
         pow(x3, 2) *
             (b + x2 + sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2))) +
         (pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)) *
             (b + x2 + sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)))) /
            (pow(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2), 1.5) *
             pow(b + x2 + sqrt(pow(a + x1, 2) + pow(b + x2, 2) + pow(x3, 2)),
                 2)))) /
      ((1 - Nu) * arma::datum::pi);

  Sigma_xy_Dx =
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
        (1 - Nu) *
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

  Sigma_xz_Dx =
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
            (-(((-a + x1) * x3 *
                (-b + x2 +
                 2 * sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2)))) /
               (pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                        2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                    1.5) *
                pow(-b + x2 +
                        sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) -
                             2 * b * x2 + pow(x2, 2) + pow(x3, 2)),
                    2))) +
             ((a + x1) * x3 *
              (-b + x2 +
               2 * sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                        2 * b * x2 + pow(x2, 2) + pow(x3, 2)))) /
                 (pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                      1.5) *
                  pow(-b + x2 +
                          sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) -
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2)),
                      2)) +
             ((-a + x1) * x3 *
              (b + x2 +
               2 * sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                        2 * b * x2 + pow(x2, 2) + pow(x3, 2)))) /
                 (pow(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                      1.5) *
                  pow(b + x2 +
                          sqrt(pow(a, 2) + pow(b, 2) - 2 * a * x1 + pow(x1, 2) +
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2)),
                      2)) -
             ((a + x1) * x3 *
              (b + x2 +
               2 * sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                        2 * b * x2 + pow(x2, 2) + pow(x3, 2)))) /
                 (pow(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                          2 * b * x2 + pow(x2, 2) + pow(x3, 2),
                      1.5) *
                  pow(b + x2 +
                          sqrt(pow(a, 2) + pow(b, 2) + 2 * a * x1 + pow(x1, 2) +
                               2 * b * x2 + pow(x2, 2) + pow(x3, 2)),
                      2))) +
        Nu * ((-b + x2) /
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

  Sigma_yz_Dx =
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

  StressTensorDueToDx(0, 0) = Sigma_xx_Dx;
  StressTensorDueToDx(0, 1) = Sigma_xy_Dx;
  StressTensorDueToDx(0, 2) = Sigma_xz_Dx;
  StressTensorDueToDx(1, 0) = Sigma_xy_Dx;
  StressTensorDueToDx(1, 1) = Sigma_yy_Dx;
  StressTensorDueToDx(1, 2) = Sigma_yz_Dx;
  StressTensorDueToDx(2, 0) = Sigma_xz_Dx;
  StressTensorDueToDx(2, 1) = Sigma_yz_Dx;
  StressTensorDueToDx(2, 2) = Sigma_zz_Dx;

  return StressTensorDueToDx;
}

}
