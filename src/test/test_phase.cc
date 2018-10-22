
#include "phase.h"
#include "subduction.h"

complex <double> Ph::comp_Wigner_d(int twoJ, int twolam1,int twolam2, double a1, double b1, double c1, double a2, double b2, double c2, int n){
    complex <double> D = 0;
    switch(n)
    {
        case 1:
            for(int i = -twoJ; i <= twoJ; i = i+2){
                D += conj(Hadron::Wigner_D(twoJ, i, twolam1, a1, b1, c1))*
                     Hadron::Wigner_D(twoJ, i,  twolam2, a2, b2, c2);
            }
            break;
            
        case 2:
            for(int i = -twoJ; i <= twoJ; i = i+2){
                D += Hadron::Wigner_D(twoJ, twolam1, i, a1, b1, c1)*
                conj(Hadron::Wigner_D(twoJ, twolam2, i, a2, b2, c2));
            }
            break;
            
        default:
            for(int i = -twoJ; i <= twoJ; i = i+2){
                D += Hadron::Wigner_D(twoJ, twolam1, i, a1, b1, c1)*
                     Hadron::Wigner_D(twoJ, i, twolam2, a2, b2, c2);
            }
            break;
            
            
    }
    
    return D;
}


Eigen::MatrixXd Ph::phaseFactor(int twoJ1, int twoJ2,  int two_abs_lam1, int two_abs_lam2,
    Eigen::Vector3d mom1, Eigen::Vector3d mom2){


  std::string LG1 = LittleGrp::generateLittleGroup(mom1);
  std::string LG2 = LittleGrp::generateLittleGroup(mom2);
    
    XMLArray::Array<int> mom2_(3); mom2_[0] = mom2(0); mom2_[1] = mom2(1); mom2_[2] = mom2(2);
    XMLArray::Array<int> mom2_ = Hadron::canonicalOrder(mom2_);
    Eigen::Vector3d mom2_ref << mom2_[0],mom2_[1],mom2_[2];
    
    
	std::vector<double> r2 = LittleGrp::refAngles(mom2_ref);
    std::vector<double> r_mom2 = LittleGrp::refAngles(mom2);
    std::vector<double> r_mom1 = LittleGrp::refAngles(mom1);


  Eigen::MatrixXd R_mom2 = Rot::eulerRotMat(rot.alpha,rot.beta,rot.gamma);
  Eigen::MatrixXd R_ref_mom2 = Rot::eulerRotMat(r2[0],r2[1],r2[2]);
  Eigen::MatrixXd R = R_ref_mom2*R.inverse();



  Eigen::MatrixXd R_mom1 = Rot::eulerRotMat(rot.alpha,rot.beta,rot.gamma);

  mom1_f << 0, mom1(0), mom1(1), mom1(2), mom1(3);
  n_mom1_f = R*mom1_f;

    Eigen::Vector3d n_mom1(3);
    n_mom1 = n_mom1_f(1), n_mom1_f(2), n_mom1_f(3);
    
    std::vector<double> r_n_mom1 = LittleGrp::refAngles(n_mom1);



  Eigen::MatrixXd R_n_mom1 = Rot::eulerRotMat(rot.alpha,rot.beta,rot.gamma);
    
    complex<double> D2 = 0;
    
    for(int i = -twoJ2; i <= twoJ2; i = i+2){
        for(int j = -twoJ2; j <= twoJ2; j = j+2){
        
            D2 += conj(Hadron::Wigner_D(twoJ2,  i, two_abs_lam2,  r2[0], r2[1], r2[2])) * Ph::comp_Wigner_d(twoJ2, i, j, r2[0], r2[1], r2[2], r_mom2[0], r_mom2[1], r_mom2[2], 2) * Hadron::Wigner_D(twoJ2,  j, two_abs_lam2, r_mom2[0], r_mom2[1], r_mom2[2]);
        
        }}
    
    
    complex<double> D1 = 0;
    
    for(int i = -twoJ1; i <= twoJ1; i = i+2){
        for(int j = -twoJ1; j <= twoJ1; j = j+2){
            
            D1 += conj(Hadron::Wigner_D(twoJ1,  i, two_abs_lam1,  r_n_mom1[0], r_n_mom1[1], r_n_mom1[2])) *  Ph::comp_Wigner_d(twoJ1, i, j, r2[0], r2[1], r2[2], r_mom2[0], r_mom2[1], r_mom2[2], 2) * Hadron::Wigner_D(twoJ1,  j, two_abs_lam1, r_mom1[0], r_mom1[1], r_mom1[2]);
            
        }}
    
    return conj(D1)*D2*R.inverse();


}



//Hadron::Wigner_D(int twoJ, 2*mu, 2*abs_lam, R1_phi, R1_theta, R1_psi);



