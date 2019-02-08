#include "k_factor_factory.h"

/* 
    Factory of functions used as Kinematic factor
*/

namespace{ 
  bool registered = false;
}

namespace{
  KFactor*  createKfacSSV( XMLReader& xml, const string& path ){ return new KfacSSV(); }  
  KFactor*  createKfacSVV( XMLReader& xml, const string& path ){ return new KfacSVV(); }   
  KFactor*  createKfacSSVwPhCorr( XMLReader& xml, const string& path ){ return new KfacSSVwPhCorr(); }  
  KFactor*  createKfacSVVwPhCorr( XMLReader& xml, const string& path ){ return new KfacSVVwPhCorr(); } 
}

namespace KFactorEnv
{ 
  
  bool registerAll(){
    bool success = true;
    
    if(!registered){
      
      success &= TheKFactorFactory::Instance().registerObject( "kin_fac_scalar_scalar_vector",
								  createKfacSSV);
      
      success &= TheKFactorFactory::Instance().registerObject( "kin_fac_scalar_vector_vector",
								  createKfacSVV);
      
      success &= TheKFactorFactory::Instance().registerObject( "kin_fac_scalar_scalar_vector_with_ph_corr",
								  createKfacSSVwPhCorr);      
      
      success &= TheKFactorFactory::Instance().registerObject( "kin_fac_scalar_vector_vector_with_ph_corr",
								  createKfacSVVwPhCorr);
      
      registered = true;
    }


    return success;
  };


}