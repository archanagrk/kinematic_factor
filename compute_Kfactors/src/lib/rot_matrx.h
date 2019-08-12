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
#include </usr/local/Eigen/Dense>
#include <vector>


using namespace std;



namespace Rot {
  
  /*  Euler matrix - converted from itpp to Eigen3 - z-y-z convention */
  Eigen::MatrixXd eulerRotMat(double alpha, double beta, double gamma);

}


#endif
