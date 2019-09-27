#ifndef __NAMING__
#define __NAMING__


//std lib
#include "phase.h"

  //**********************************************************************************************************************

    /* Naming Scheme */

  //**********************************************************************************************************************

namespace naming {
  string name(int npt, Ph::tripKey two_abs_lam, Vector3d mom1 , Vector3d mom_curr, Vector3d mom3, irrep_label rep1,
                                             irrep_label rep_curr, irrep_label rep3, string LG1, string LG_curr, string LG3);
}


#endif