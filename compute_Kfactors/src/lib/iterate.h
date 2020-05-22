#ifndef __ITER__
#define __ITER__


//std lib
#include <cmath>
#include "math.h"
#include </usr/local/Eigen/Dense>
#include <float.h>
#include <vector>

#include "io/adat_io.h"
#include "io/adat_xml_group_reader.h"

#include "hadron/irrep_util.h"

using namespace Eigen;


  //**********************************************************************************************************************

    /* function to iterate over each index of three-momenta */

  //**********************************************************************************************************************

namespace iter {
  std::vector<Vector3d> itermom(double max_mom, double min_mom, ADAT::Array1dO<ADAT::Array1dO<int>> omit_mom, bool canonical);
}


#endif