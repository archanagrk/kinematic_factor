#ifndef __LITTLEGROUP__
#define __LITTLEGROUP__



//std lib
#include <vector>
#include <map>
#include <string>
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
using namespace ENSEM;
using namespace Hadron;
using namespace Eigen;


namespace LittleGrp{
  string generateLittleGroup(Eigen::Vector3d& mom_);
  //std::vector<double> refAngles(string little_group);
  std::vector<double> refAngles(Eigen::Vector3d mom1);
}


namespace{ bool is_equal(double a, double b); }
namespace{ const double PI = (atan(double(1)) * double(4.0)); }


#endif
