#ifndef __K_FACTOR_PIGAMMARHO_H__
#define __K_FACTOR_PIGAMMARHO_H__

#include "kfac_params.h"
#include "levi_civita.h"



/* build a factory of KFactor types */
/* many will be specific to the matrix elements */
/* the AvgFit passed in contains the data and the function, 
   so the quality can compute all kinds of stuff */

/*
  the kinematic factor constructor can depend on fixed parameters
  so an xml control is required
*/


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// base class of Kinematic Factor
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class KFactor{
 public:
  virtual ~KFactor(){};

  virtual vector<complex<double>> operator()( const KFacParams& params ) const = 0;
  /* return a double which is the Kinematic Factor */

  virtual string name() const = 0;
  /* return the name of the factor -- should match the string in the factory */
  
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// some common Kinematic Factors
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class KfacSVV : public KFactor{
 public:
  vector<complex<double>> operator()( const KFacParams& params ) const;
  string name() const { return "kin_fac_scalar_vector_vector"; } 
};


class KfacSSV : public KFactor{
 public:
  vector<complex<double>> operator()( const KFacParams& params ) const;
  string name() const { return "kin_fac_scalar_scalar_vector"; } 
};



#endif
