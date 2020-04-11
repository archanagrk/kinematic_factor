#ifndef __KFACTOR__
#define __KFACTOR__

//std lib
#include </usr/local/Eigen/Dense>
#include <float.h>

#include "levi_civita.h"
#include "kfac_utils.h"
#include "phase.h"
#include "pol_vec.h"

#define ULP_N 8

namespace KFac {

  bool nearlyequal(complex<double>& a, complex<double>& b);

  std::complex<double> KFacSVV(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub3);
    
  std::complex<double> KFacSSV(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >&  Sub1 , map< int, Eigen::MatrixXcd >&  SubCurr , map< int, Eigen::MatrixXcd >&  Sub3 );

  std::complex<double> KFacold(Eigen::VectorXd& qp, Eigen::VectorXd& qm, Eigen::MatrixXcd& Sub1 , Eigen::MatrixXcd& SubCurr , Eigen::MatrixXcd& Sub3 );
    
  Eigen::MatrixXcd subPhSum(map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2, Ph::phChars& phase);
    
  Eigen::MatrixXcd subSum(map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2);

  std::complex<double> KinematicFactorwithPhase(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub3, Ph::phChars& phase);
    
  std::complex<double> KinematicFactorwithPhase_j0(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub3, Ph::phChars& phase);
  
  std::complex<double> KFacwoS(Eigen::VectorXd& qp, Eigen::VectorXd& qm, std::vector<double>& r3, double& mom3_sq, const int& two_lam3, double& m3_sq, std::vector<double>& r_curr, double& mom_curr_sq, const int& two_lam_curr, double& m_curr_sq);
  
  std::complex<double> KFacwoSwP(Eigen::VectorXd& qp, Eigen::VectorXd& qm, const int& two_lam1, std::vector<double>& r3, double& mom3_sq, const int& two_lam3, double& m3_sq, std::vector<double>& r_curr, double& mom_curr_sq, const int& two_lam_curr, double& m_curr_sq, Ph::phChars& phase);

}
#endif


