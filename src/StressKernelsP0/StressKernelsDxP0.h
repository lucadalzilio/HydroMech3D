//
// Created by Federico Ciardo on 16.08.21. All rights reserved.
//

// Inclusion from IL library
#include <armadillo>

#ifndef INC_3DEQSIM_SRC_ELASTICKERNELSP0_H
#define INC_3DEQSIM_SRC_ELASTICKERNELSP0_H

namespace EQSim {

arma::mat StressTensorDueToDDxOnSingleEltP0(double &a, double &b,
                                            double &x1, double &x2,
                                            double &x3, double &Nu,
                                            double &G);

}
#endif  // INC_3DEQSIM_SRC_ELASTICKERNELSP0_H
