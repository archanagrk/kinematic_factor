#ifndef __POLVEC__
#define __POLVEC__

#include "kfac_utils.h"
#include "rot_matrx.h"
#include "little_group.h"
#include "subduction.h"

  //**********************************************************************************************************************

namespace PolVec {

  /* Pol vector along the z-axis from helicity ops paper */
  Eigen::MatrixXcd getPolz4(double& mom_sq, const int& two_helicity, double& mass_sq);
  Eigen::MatrixXcd getPolz3(double& mom_sq, const int& two_helicity, double& mass_sq);

  /* Get the polarization by rotating from the z-axis */
  Eigen::MatrixXcd getPol4(double& mom_sq, const int& two_helicity, double& mass_sq, double& phi, double& theta, double& psi);
  Eigen::MatrixXcd getPol3(double& mom_sq, const int& two_helicity, double& mass_sq, double& phi, double& theta, double& psi);
}

  //**********************************************************************************************************************
  
#endif
