#ifndef __ROTMATRX__
#define __ROTMATRX__

#include "kfac_utils.h"

//std lib
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "math.h"
#include <stdio.h>
#include <Eigen/Dense>
#include <vector>


using namespace std;

namespace Rot {
  Eigen::MatrixXd eulerRotMat(double alpha, double beta, double gamma);
  Eigen::MatrixXd getRotMat(Eigen::Vector3d mom2, Eigen::Vector3d mom1);
  std::vector<double> euAng(Eigen::MatrixXd R);
}


#endif
