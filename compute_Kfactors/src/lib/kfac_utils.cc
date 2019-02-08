/* Some utility codes */

#include "kfac_utils.h"


  //**********************************************************************************************************************

  //**********************************************************************************************************************

namespace KfUt{

  /* Truncate small numbers  */

  double truncate(double num,int precision){
    if(abs(num) < pow(10,-precision)){ num = double(0);}
    return num;
  }

    //**********************************************************************************************************************
  
  /* Convert Eigen and Array1dO to array data type */

  XMLArray::Array<int> ToArray::toArray(Eigen::Vector3d input){
    XMLArray::Array<int> in_1(3); in_1[0] = input(0); in_1[1] = input(1); in_1[2] = input(2);
    return(in_1);
    
  }
  
  XMLArray::Array<int> ToArray::toArray(Array1dO<int> input){
    XMLArray::Array<int> in_1(3); in_1[0] = input[1]; in_1[1] = input[2]; in_1[2] = input[3];
    return(in_1);
  }

  //**********************************************************************************************************************

  /* Gmunu matrix */

  Eigen::MatrixXcd Gmunu(){
    
    typedef std::complex<double> cd;
    Eigen::MatrixXcd Gmunu = Eigen::MatrixXcd::Zero(4,4);
  
    Gmunu(0,0) = cd(1,0);
    Gmunu(1,1) = cd(-1,0);
    Gmunu(2,2) = cd(-1,0);
    Gmunu(3,3) = cd(-1,0);
    
    return Gmunu;
  }

//**********************************************************************************************************************

};


