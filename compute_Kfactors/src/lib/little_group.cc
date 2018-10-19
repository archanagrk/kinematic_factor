/* Find Little Groups - Jo's code  */

#include "little_group.h"

namespace{
  bool is_equal(double a, double b){
    return std::abs(a - b) < std::numeric_limits<float>::epsilon();
  }
}



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
        littleGroup = "C4";         // n m 0
    }
  }
  else{
    if ( (momCan[0] == momCan[1]) && (momCan[0] == momCan[2]) )   // n n n
      littleGroup = "D3";
    else if ( (momCan[0] == momCan[1]) || (momCan[0] == momCan[2]) || (momCan[1] == momCan[2]) )   // m m n
      littleGroup = "C4";
    else   // n m p
      littleGroup = "C2";
  }
  return littleGroup;
};

/* Get the reference angles for each LG - copied from Jo - Appendix E Table VI - Helicity ops for Mesons - z-y-z  Jacob-Wick convention */


/*std::vector<double> refAngles(string little_group){
  std::vector<double> ref;
  double R1_phi, R1_theta, R1_psi; //reference rotation angles

  if(little_group == "Oh"){ R1_phi = 0.0; R1_theta = 0.0; R1_psi = 0.0;}
  else if(little_group == "D4"){ R1_phi = 0.0; R1_theta = 0.0; R1_psi = 0.0;} // (00n)
  else if(little_group == "D2"){ R1_phi = PI/2.0; R1_theta = PI/4.0; R1_psi = -PI/2.0;}// (0nn)
  else if(little_group == "D3"){ R1_phi = PI/4.0; R1_theta = acos(1.0/sqrt(3.0)); R1_psi = 0.0;}// (nnn)
  else if(little_group == "C4"){ cerr << "C4 not coded" << endl; exit(1); }// ???????
  else if(little_group == "C2"){ cerr << "C2 not coded" << endl; exit(1); }// ???????
  ref.push_back(R1_phi); ref.push_back(R1_theta); ref.push_back(R1_psi);

  return ref;
  }*/

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
    return ref;
 }

};




