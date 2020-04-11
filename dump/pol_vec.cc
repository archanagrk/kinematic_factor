  //**********************************************************************************************************************

    /* To compute the polarization given the momentum and helicity - Convention from Appendix B of Helicity Ops for Mesons  */

  //**********************************************************************************************************************

#include "pol_vec.h"

namespace PolVec {


  //**********************************************************************************************************************

  /* Polarization  4 vector */
  
  //**********************************************************************************************************************

  Eigen::MatrixXcd getPolz4(double& mom_sq, const int& two_helicity, double& mass_sq)

  {
  
    Eigen::MatrixXcd pol_z = Eigen::MatrixXcd::Zero(4,1);    // For spin 1 particles extra polarization degree of freedom
    
    typedef std::complex<double> cd;
    
    if(two_helicity == 0){
      
      if(mass_sq == 0){pol_z << cd(0,0),cd(0,0),cd(0,0),cd(0,0);}            // Massless particles have only two physical polarizations
        
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

  /* Polarization 3 vector */

  //**********************************************************************************************************************

  Eigen::MatrixXcd getPolz3(double& mom_sq, const int& two_helicity, double& mass_sq)

  {
    
    Eigen::MatrixXcd pol_z = Eigen::MatrixXcd::Zero(3,1);    // For spin 1 particles extra polarization degree of freedom
    
    typedef std::complex<double> cd;
    
    if(two_helicity == 0){
      if(mass_sq == 0){pol_z << cd(0,0),cd(0,0),cd(0,0);}            // Massless particles have only two physical polarizations
      
      else{
        if(mom_sq != 0){
        
          std::complex<double> E = std::sqrt(std::complex<double>(mom_sq+mass_sq,0));
            
          pol_z << cd(0,0),cd(0,0),E/std::sqrt(std::complex<double>(mass_sq,0));   // from helicity ops paper
                    
        }
        else{pol_z << cd(0,0),cd(0,0),cd(1,0);}
      }}
    
    else if(two_helicity == 2){pol_z <<  cd(-sqrt(0.5),0),cd(0,-sqrt(0.5)),cd(0,0);}        // from helicity ops paper
    else if(two_helicity == -2){pol_z << cd(sqrt(0.5),0),cd(0,-sqrt(0.5)),cd(0,0);}
    
    else{cerr << "Not valid" << endl; exit(1);}
    return pol_z;
  }
  
  //**********************************************************************************************************************

  /* Get the polarization by rotating from the z-axis */

  //**********************************************************************************************************************


 Eigen::MatrixXcd getPol4(double& mom_sq, const int& two_helicity, double& mass_sq, double& phi, double& theta, double& psi)
 
 {
   Eigen:: MatrixXcd pol = Eigen::MatrixXcd::Zero(4,1);
   
   Eigen::MatrixXcd pol_z = Eigen::MatrixXcd::Zero(4,1);    // For spin 1 particles extra polarization degree of freedom
   
   pol_z = getPolz4(mom_sq,two_helicity,mass_sq);

   pol = Rot::eulerRotMat(phi,theta,psi)*pol_z;                              // multiplies by the euler matrix to convert p_ref to p_canonical
   

   return pol;
 }

  //**********************************************************************************************************************

  /* Get the polarization by rotating from the z-axis */

  //**********************************************************************************************************************
  
  Eigen::MatrixXcd getPol3(double& mom_sq, const int& two_helicity, double& mass_sq, double& phi, double& theta, double& psi)
  
  {
    
    Eigen:: MatrixXcd pol = Eigen::MatrixXcd::Zero(3,1);
    
    Eigen::MatrixXcd pol_z3 = Eigen::MatrixXcd::Zero(3,1);
    
    for(int i = -2; i <= 2; i = i+2){
      pol_z3 = getPolz3(mom_sq,i,mass_sq);
      pol = pol + Hadron::Wigner_D(2,i,two_helicity,phi, theta, psi)*pol_z3;
    }
    // multiplies by the euler matrix to convert p_ref to p_canonical
    

    return pol;
    
  }
  

};
