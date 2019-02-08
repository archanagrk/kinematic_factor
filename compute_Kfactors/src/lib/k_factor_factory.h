#ifndef __KFACTOR_FACTORY_H__
#define __KFACTOR_FACTORY_H__

/* keep adding to this list of function containing header files */

#include "kfactor_pigammarho.h"

typedef
SingletonHolder< ObjectFactory<KFactor, 
  string, 
  TYPELIST_2(     XMLReader&, const string&),
  KFactor* (*)(XMLReader&, const string&),
  StringFactoryError> >TheKFactorFactory;

namespace KFactorEnv
{ 
  bool registerAll();
}



#endif