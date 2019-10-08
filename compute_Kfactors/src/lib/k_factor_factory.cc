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
  KFactor*  createKfacSSS( XMLReader& xml, const string& path ){ return new KfacSSS(); }  
  KFactor*  createKfacSVS( XMLReader& xml, const string& path ){ return new KfacSVS(); }   

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
      
      success &= TheKFactorFactory::Instance().registerObject( "kin_fac_scalar_scalar_scalar",
								  createKfacSSS);
      
      success &= TheKFactorFactory::Instance().registerObject( "kin_fac_scalar_vector_scalar",
								  createKfacSVS);

      
      registered = true;
    }


    return success;
  };


}