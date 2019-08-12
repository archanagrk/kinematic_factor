#ifndef __LITTLEGROUP__
#define __LITTLEGROUP__



//std lib
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include "math.h"
#include </usr/local/Eigen/Dense>


//kfac lib
#include "kfac_utils.h"

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


  //**********************************************************************************************************************


namespace LittleGrp{

  /* Find the little group based on the momentum */
  string generateLittleGroup(Eigen::Vector3d& mom_);

/* Angles that take pz to the given p */
/*  Get the reference angles for each LG - Appendix E Table VI - Helicity ops for Mesons - z-y-z  Jacob-Wick convention */
  std::vector<double> refAngles(Eigen::Vector3d mom1);
}

/* function to check equality */
namespace{ bool is_equal(double a, double b); }


  //**********************************************************************************************************************

#endif
