#ifndef __PHASE__
#define __PHASE__

#include "rot_matrx.h"
#include "little_group.h"
#include "kfac_utils.h"
#include "subduction.h"

  //**********************************************************************************************************************


namespace Ph {

  typedef std::tuple<int, int, int> tripKey;

  //**********************************************************************************************************************

  //struct

  struct phChars{
    Eigen::Vector3d mom2;   /* new mom */
    Eigen::Vector3d mom1;        /* new mom */
    map < Ph::tripKey , complex<double>> lam_phase; /* helicity1/2 and phase */
    Eigen::MatrixXcd r; /* R rot matrix */
        bool operator<(const phChars &rhs) const; /*usual map label problems */
        //bool operator!=(const irrep_label &rhs) const; /*usual map label problems */
    };

  //**********************************************************************************************************************

  /* Compute the phase factor after rotation of a helicity operator or state from Appendix A of Shulz paper */

  Ph::phChars phaseFactor(int twoJ1, int twoJ2, int twoJCurr, Eigen::Vector3d mom1, Eigen::Vector3d mom2, bool compute, bool right);



  /* Composition of two Wigner-D matrices R'R, R'^-1R, R'R^-1 */

  std::complex <double> comp_Wigner_d(int twoJ, int twolam1,int twolam2, double a1, double b1, double c1, double a2, double b2, double c2, int n);
    
  //**********************************************************************************************************************
 
}


#endif
