#ifndef __KFACTOR__
#define __KFACTOR__

//std lib
#include <Eigen/Dense>

#include "levi_civita.h"
#include "kfac_utils.h"
#include "phase.h"


namespace KFac {

  std::complex<double> KinematicFactor(Eigen::VectorXd& qp, Eigen::VectorXd& qm, Eigen::MatrixXcd& Sub1 , Eigen::MatrixXcd& SubCurr , Eigen::MatrixXcd& Sub3 );
    Eigen::MatrixXcd subPhSum(map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2, Ph::phChars& phase);
  std::complex<double> KinematicFactorwithPhase(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub3, Ph::phChars& phase);

}
#endif


