#ifndef __K_FACTOR_PIGAMMAPI_H__
#define __K_FACTOR_PIGAMMAPI_H__

#include "kfac_params.h"
#include "levi_civita.h"
#include "kfactor_pigammarho.h"




/* build a factory of KFactor types */
/* many will be specific to the matrix elements */

/*
  the kinematic factor constructor can depend on fixed parameters
  so an xml control is required
*/


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// some other Kinematic Factors
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class KfacSVS : public KFactor{
 public:
  vector<complex<double>> operator()( const KFacParams& params ) const;
  string name() const { return "kin_fac_scalar_vector_scalar"; } 
};

class KfacSSS : public KFactor{
 public:
  vector<complex<double>> operator()( const KFacParams& params ) const;
  string name() const { return "kin_fac_scalar_scalar_scalar"; } 
};



#endif
