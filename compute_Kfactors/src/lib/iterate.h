#ifndef __ITER__
#define __ITER__


//std lib
#include <cmath>
#include "math.h"
#include </usr/local/Eigen/Dense>
#include <float.h>
#include <vector>

using namespace Eigen;


  //**********************************************************************************************************************

    /* Naming Scheme */

  //**********************************************************************************************************************

namespace iter {
  std::vector<Vector3d> itermom(double max_mom);
}


#endif