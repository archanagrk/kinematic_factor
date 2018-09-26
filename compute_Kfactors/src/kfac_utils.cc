/* Some utility codes */

#include "kfac_utils.h"


/* Truncate small numbers  */


namespace KfUt{
  double truncate(double num,int precision){
    if(abs(num) < pow(10,-precision)){ num = double(0);}
    return num;
  }}

namespace {  const double PI = (atan(double(1)) * double(4.0)); }

