#ifndef __KFACUTILS__
#define __KFACUTILS__

//std lib
#include <cmath>
#include "math.h"
#include <Eigen/Dense>


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


using namespace std;


namespace KfUt{
  double truncate(double num,int precision);
  //XMLArray::Array<int> toArray(Eigen::Vector3d input);
  
  class ToArray {
  public:
    static XMLArray::Array<int> toArray(Eigen::Vector3d input);
    static XMLArray::Array<int> toArray(Array1dO<int> input);
  };
  
}


#endif
