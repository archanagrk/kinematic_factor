#include "rot_matrx.h"

namespace Rot {

  //**********************************************************************************************************************

    /* Rotation Matrix using euler angles z-y-z convention */

  //**********************************************************************************************************************

Eigen::MatrixXd eulerRotMat(double alpha, double beta, double gamma)
{
  // to act upon cartesian 4-vectors
  // R(a,b,g) = Rz(a) Ry(b) Rz(g)


  Eigen::MatrixXd Rzg = Eigen::MatrixXd::Zero(4,4);
  Rzg(0,0) = 1;
  Rzg(1,1) = cos(gamma); Rzg(1,2) = -1 * sin(gamma);
  Rzg(2,1) = -1 * (Rzg(1,2)); Rzg(2,2) = Rzg(1,1);
  Rzg(3,3) = 1;

  Eigen::MatrixXd Ryb = Eigen::MatrixXd::Zero(4,4);
  Ryb(0,0) = 1;
  Ryb(1,1) = cos(beta); Ryb(1,3) = sin(beta);
  Ryb(2,2) = 1;
  Ryb(3,1) = -1 * (Ryb(1,3)); Ryb(3,3) = Ryb(1,1);

  Eigen::MatrixXd Rza = Eigen::MatrixXd::Zero(4,4);
  Rza(0,0) = 1;
  Rza(1,1) = cos(alpha); Rza(1,2) = -1 * sin(alpha);
  Rza(2,1) = -1 * (Rza(1,2)); Rza(2,2) = Rza(1,1);
  Rza(3,3) = 1;

  return Rza * Ryb * Rzg;
 }

  
};


