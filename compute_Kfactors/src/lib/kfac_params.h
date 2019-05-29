#ifndef __KFPARAMS__
#define __KFPARAMS__

//std lib
#include <cmath>
#include "math.h"
#include <Eigen/Dense>
#include <float.h>

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


  /* Struct to store the subductions and moms */
  
  class KFacParams
  {
    public:
      virtual ~KFacParams() {};                                    /* virtual destructor */

      map< int, Eigen::MatrixXcd > Sub1;                  /* source subduction */

      map< int, Eigen::MatrixXcd > SubCurr;                /* current subdcution */

      map< int, Eigen::MatrixXcd > Sub3;                   /* sink subduction */

      Ph::phChars phase;                               /* Phase = 1 by default  */

      VectorXd  qp;                                       /* source mom + sink mom */

      VectorXd  qm;                                       /* sink mom - source mom */

      KFacParams(map< int, Eigen::MatrixXcd >, map< int, Eigen::MatrixXcd >, map< int, Eigen::MatrixXcd >, Ph::phChars, VectorXd, VectorXd);

      virtual MatrixXcd subPhSum() const;  /* returns the phase sum */

      virtual Ph::tripKey two_abs_lam() const; 

  };

  


#endif
