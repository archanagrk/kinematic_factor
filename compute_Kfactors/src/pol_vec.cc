/* To compute the polarization given the momentum and helicity - Convention from Appendix B of Helicity Ops for Mesons  */

#include "pol_vec.h"

namespace PolVec {

 Eigen::MatrixXcd getPolarization(double& mom_sq, const int& two_helicity, double& mass_sq, double& phi, double& theta, double& psi){
   Eigen:: MatrixXcd pol = Eigen::MatrixXcd::Zero(4,1);
    Eigen::MatrixXcd pol_z = Eigen::MatrixXcd::Zero(4,1);
                                                                 // For spin 1 particles extra polarization degree of freedom
   if(two_helicity == 0){
     if(mass_sq == 0){pol_z << (0,0),(0,0),(0,0),(0,0);}            // Massless particles have only two physical polarizations

     else{
        if(mom_sq != 0){
          double energy = sqrt(mom_sq+mass_sq);
          pol_z << (sqrt(mom_sq/mass_sq),0),(0,0),(0,0),(energy/sqrt(mass_sq),0);   // k.pol = 0
          }
        else{pol_z << (0,0),(0,0),(0,0),(1,0);}
   }}

   else if(two_helicity == 2){pol_z <<  (0,0),(-sqrt(0.5),0),(0,-sqrt(0.5)),(0,0);}        // k.pol = 0
   else if(two_helicity == -2){pol_z << (0,0),(sqrt(0.5),0),(0,-sqrt(0.5)),(0,0);}

   else{cerr << "Not valid" << endl; exit(1);}
   pol = Rot::eulerRotMat(phi,theta,psi)*pol_z;                              // multiplies by the euler matrix to convert p_ref to p_canonical
   //cout << "\n" << "pol_s" << pol << "pol_e" << "\n" ;
   return pol;

 }

};
