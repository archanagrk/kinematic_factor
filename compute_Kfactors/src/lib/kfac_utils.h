#ifndef __KFACUTILS__
#define __KFACUTILS__

//std lib
#include <cmath>
#include "math.h"
#include </usr/local/Eigen/Dense>

//some kfac libs
#include "phase.h"


//adat
#include <adat/handle.h>
#include "hadron/irreps_cubic_factory.h"
#include "hadron/irreps_cubic_helicity_factory.h"
#include "hadron/irrep_util.h" //adat_devel
#include "ensem/ensem.h"



using namespace std;
using namespace ADAT;
using namespace ADATXML;
using namespace Hadron;
using namespace Eigen;


  //**********************************************************************************************************************


namespace KfUt{

  /* Truncate small numbers  */

  double truncate(double num,int precision);

  /* Gmunu matrix */

  Eigen::MatrixXcd Gmunu();

  /* Convert Eigen and Array1dO to array data type */
  
  class ToArray {
  public:
    static XMLArray::Array<int> toArray(Eigen::Vector3d input);
    static XMLArray::Array<int> toArray(Array1dO<int> input);
  };

  
}

 /* Define π */

namespace { const double PI = (atan(double(1)) * double(4.0));}


#endif
