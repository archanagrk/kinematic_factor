#ifndef __POLVEC__
#define __POLVEC__

#include "rot_matrx.h"


namespace PolVec {
  Eigen::MatrixXcd getPolarization(double& mom_sq, const int& two_helicity, double& mass_sq, double& phi, double& theta, double& psi);
}


#endif
