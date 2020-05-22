/* Find Little Groups */

#include "little_group.h"

  //**********************************************************************************************************************

    /* function to check equality */

  //**********************************************************************************************************************

namespace{
  bool is_equal(double a, double b){
    return std::abs(a - b) < std::numeric_limits<float>::epsilon();
  }
}

  //**********************************************************************************************************************

    /* Find the little group based on the momenta */

  //**********************************************************************************************************************

namespace LittleGrp{

string generateLittleGroup(Eigen::Vector3d& mom_)
{ 
  XMLArray::Array<int> mom(3); mom[0] = mom_(0); mom[1] = mom_(1); mom[2] = mom_(2);
  XMLArray::Array<int> momCan = Hadron::canonicalOrder(mom);
  std::string littleGroup = "";
  if (momCan[2] == 0){
    if (momCan[1] == 0){
      if (momCan[0] == 0)   // 0 0 0
        littleGroup = "Oh"; 
      else                  // n 0 0
        littleGroup = "D4";
    }
    else{
      if (momCan[0] == momCan[1])   // n n 0
        littleGroup = "D2";
      else
        littleGroup = "C4nm0";         // n m 0
    }
  }
  else{
    if ( (momCan[0] == momCan[1]) && (momCan[0] == momCan[2]) )   // n n n
      littleGroup = "D3";
    else if ( (momCan[0] == momCan[1]) || (momCan[0] == momCan[2]) || (momCan[1] == momCan[2]) )   // m m n
      littleGroup = "C4nnm";
    else   // n m p
      littleGroup = "C2";
  }
  
  return littleGroup;
};



  //**********************************************************************************************************************

/* Angles that take pz to the given p */
/*  Get the reference angles for each LG - Appendix E Table VI - Helicity ops for Mesons - z-y-z  Jacob-Wick convention */

  //**********************************************************************************************************************

std::vector<double> refAngles(Eigen::Vector3d mom1){
    std::vector<double> ref;
    double R_phi, R_theta, R_psi;
    XMLArray::Array<int> mom(3);



    mom[0] = round(mom1(0));
    mom[1] = round(mom1(1));
    mom[2] = round(mom1(2));
    


    if(is_equal(mom[0],0.0) && is_equal(mom[1],0.0) && is_equal(mom[2],0.0) ){ref.push_back(0.0); ref.push_back(0.0); ref.push_back(0.0);}
    else{
      Hadron::CubicCanonicalRotation_t rot;

      rot = Hadron::cubicCanonicalRotation(mom);
      ref.push_back(rot.alpha);ref.push_back(rot.beta);ref.push_back(rot.gamma);
    }

    //cout << "refs" << ref[0] <<  "refs"  << ref [1] <<  "refs" << ref[2] << "refs" << "\n";
    
    return ref;
 }

};




