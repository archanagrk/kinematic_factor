  //**********************************************************************************************************************

    /* To compute the polarization given the momentum and helicity - Convention from Appendix B of Helicity Ops for Mesons  */

  //**********************************************************************************************************************

#include "pol_vec.h"

namespace PolVec {


  //**********************************************************************************************************************

  /* Polarization  4 vector in the z direction */
  
  //**********************************************************************************************************************

  Eigen::MatrixXcd getPolz4(double& mom_sq, const int& two_helicity, double& mass_sq,bool& curr)

  {
  
    Eigen::MatrixXcd pol_z = Eigen::MatrixXcd::Zero(4,1);    // For spin 1 particles extra polarization degree of freedom
    
    typedef std::complex<double> cd;
    
    if(two_helicity == 0){
      
      if(mass_sq == 0 && !curr){pol_z << cd(0,0),cd(0,0),cd(0,0),cd(0,0);}            // Massless particles have only two physical polarizations
        
        else{
          if(mom_sq != 0){

            std::complex<double> E = std::sqrt(std::complex<double>(mom_sq+mass_sq,0));
          
            pol_z << std::sqrt(std::complex<double>(mom_sq/mass_sq,0)) ,cd(0,0),cd(0,0),E/std::sqrt(std::complex<double>(mass_sq,0));   // k.pol = 0
            
          }

          else{pol_z << cd(0,0),cd(0,0),cd(0,0),cd(1,0);}
        }
      }
        
        else if(two_helicity == 2){pol_z <<  cd(0,0),cd(-sqrt(0.5),0),cd(0,-sqrt(0.5)),cd(0,0);}        // k.pol = 0
        else if(two_helicity == -2){pol_z << cd(0,0),cd(sqrt(0.5),0),cd(0,-sqrt(0.5)),cd(0,0);}
      
        else{cerr << "Not valid" << endl; exit(1);}
        return pol_z;
  }
  

  //**********************************************************************************************************************

  /* Get the polarization by rotating from the z-axis */

  //**********************************************************************************************************************


 Eigen::MatrixXcd getPol4(double& mom_sq, const int& two_helicity, double& mass_sq, double& phi, double& theta, double& psi, bool curr)
 
 {
  Eigen:: MatrixXcd pol = Eigen::MatrixXcd::Zero(4,1);
  Eigen::MatrixXcd pol_z = Eigen::MatrixXcd::Zero(4,1);    // For spin 1 particles extra polarization degree of freedom

  if(curr){
    double zero = 0.0;

    // int J = 2;
    // complex<double> z_i(0.,1.);
    // complex<double> WignerD;

    // for(int m = -2; m <= 2; m = m+2){
    //   //pol_z = -z_i * (getPolz4(zero,m, mass_sq)).conjugate();
    //   pol_z = -1 * (getPolz4(zero,m, mass_sq));
    //   WignerD = conj(Hadron::Wigner_D(J, m, two_helicity, phi, theta, psi));
    //   pol += WignerD*pol_z;
    // }

    pol_z = getPolz4(zero, two_helicity, mass_sq,curr); //for the current the pol vector is just to rotation to circular basis
    pol = Rot::eulerRotMat(phi,theta,psi)*pol_z; 


    return pol;
  }
  else{
   
   pol_z = getPolz4(mom_sq,two_helicity,mass_sq,curr);

   pol = Rot::eulerRotMat(phi,theta,psi)*pol_z;             // multiplies by the euler matrix to convert p_ref to p_canonical
   
  }


  return pol;
 }
  

};
