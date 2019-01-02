/* Some utility codes */

#include "kfac_utils.h"


/* Truncate small numbers  */


namespace KfUt{
  double truncate(double num,int precision){
    if(abs(num) < pow(10,-precision)){ num = double(0);}
    return num;
  }
  
  XMLArray::Array<int> ToArray::toArray(Eigen::Vector3d input){
    XMLArray::Array<int> in_1(3); in_1[0] = input(0); in_1[1] = input(1); in_1[2] = input(2);
    return(in_1);
    
  }
  
  XMLArray::Array<int> ToArray::toArray(Array1dO<int> input){
    XMLArray::Array<int> in_1(3); in_1[0] = input[1]; in_1[1] = input[2]; in_1[2] = input[3];
    return(in_1);
  }
  
//  class ToArray{
//    
//  public:
//     toArray(Eigen::Vector3d input){
//      XMLArray::Array<int> in_1(3); in_1[0] = input(0); in_1[1] = input(1); in_1[2] = input(2);
//      return(in_1);
//    }
//      
//   static XMLArray::Array<int> toArray(Array1dO<int> input){
//        XMLArray::Array<int> in_1(3); in_1[0] = input[1]; in_1[1] = input[2]; in_1[2] = input[3];
//        return(in_1);
//    }
//
//
//  };

};


namespace {  const double PI = (atan(double(1)) * double(4.0)); }

