
#include "phase.h"

std::complex <double> Ph::comp_Wigner_d(int twoJ, int twolam1, int twolam2,  double a1, double b1, double c1, double a2, double b2, double c2, int n){
  std::complex <double> D = 0;
    switch(n)
    {
        case 1:
            for(int i = -twoJ; i <= twoJ; i = i+2){
                D += conj(Hadron::Wigner_D(twoJ, i, twolam1, a1, b1, c1))*
                     Hadron::Wigner_D(twoJ, i, twolam2, a2, b2, c2);
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


Ph::phChars Ph::phaseFactor(int twoJ1, int twoJ2, int twoJCurr, Eigen::Vector3d mom1, Eigen::Vector3d mom2){
    
    Ph::phChars out;
    map < Ph::tripKey , complex<double>> phase_hel;
    

    XMLArray::Array<int> mom2_r(3); mom2_r[0] = mom2(0); mom2_r[1] = mom2(1); mom2_r[2] = mom2(2);
    mom2_r = Hadron::canonicalOrder(mom2_r);
    Eigen::Vector3d mom2_ref;
    mom2_ref  << mom2_r[0],mom2_r[1],mom2_r[2];
    
    
	std::vector<double> r2 = LittleGrp::refAngles(mom2_ref);
    std::vector<double> r_mom2 = LittleGrp::refAngles(mom2);
    std::vector<double> r_mom1 = LittleGrp::refAngles(mom1);


  Eigen::MatrixXd R_mom2 = Rot::eulerRotMat(r_mom2[0],r_mom2[1],r_mom2[2]);
  Eigen::MatrixXd R_ref_mom2 = Rot::eulerRotMat(r2[0],r2[1],r2[2]);
  Eigen::MatrixXd R = R_ref_mom2*R_mom2.inverse();



  //Eigen::MatrixXd R_mom1 = Rot::eulerRotMat(r_mom1[0],r_mom1[1],r_mom1[2]);

  Eigen::VectorXd  mom1_f(4); mom1_f << 0.0 , mom1(0), mom1(1), mom1(2);
  Eigen::VectorXd n_mom1_f = R*mom1_f;

    Eigen::Vector3d n_mom1;
    n_mom1 << round(n_mom1_f(1,0)), round(n_mom1_f(2,0)), round(n_mom1_f(3,0));


    std::vector<double> r_n_mom1 = LittleGrp::refAngles(n_mom1);


    Eigen::Vector3d n_mom_curr; n_mom_curr =  mom2_ref - n_mom1;
    Eigen::Vector3d mom_curr; mom_curr = mom2 - mom1;


    std::vector<double> r_mom_curr = LittleGrp::refAngles(mom_curr);
    std::vector<double> r_n_mom_curr = LittleGrp::refAngles(n_mom_curr);
 

   Eigen::MatrixXd R_mom_curr = Rot::eulerRotMat(r_mom_curr[0],r_mom_curr[1],r_mom_curr[2]);
   Eigen::MatrixXd R_n_mom_curr = Rot::eulerRotMat(r_n_mom_curr[0],r_n_mom_curr[1],r_n_mom_curr[2]);

    Eigen::MatrixXd R_n_mom1 = Rot::eulerRotMat(r_n_mom1[0],r_n_mom1[1],r_n_mom1[2]);
    
    cout << R_n_mom1 << "\n";
    cout << R_n_mom_curr << "\n";
    cout << R_ref_mom2 << "\n";

     Eigen::MatrixXcd ans;
     int app = 0;
     
     double sqw = n_mom_curr.squaredNorm();

     double ma = 0.0151; 

     ans = PolVec::getPolarization(sqw, app , ma , r_n_mom_curr[0],r_n_mom_curr[1],r_n_mom_curr[2]);
     cout << ans << "\n";

  Ph::tripKey lambd; 
    
  for(int two_abs_lam1 = -twoJ1; two_abs_lam1 <= twoJ1; two_abs_lam1 = two_abs_lam1 +2 ){
      for(int two_abs_lam2 = -twoJ2; two_abs_lam2 <= twoJ2; two_abs_lam2 = two_abs_lam2 +2 ){
        for(int two_abs_lamCurr = -twoJCurr; two_abs_lamCurr <= twoJCurr; two_abs_lamCurr = two_abs_lamCurr +2 ){
          
    complex<double> D2 = 0;
    complex<double> D1 = 0;
    complex<double> DCurr = 0;
      
      
    for(int i = -twoJ2; i <= twoJ2; i = i+2){
        for(int j = -twoJ2; j <= twoJ2; j = j+2){
        
            D2 += conj(Hadron::Wigner_D(twoJ2, i, two_abs_lam2, r2[0], r2[1], r2[2])) * Ph::comp_Wigner_d(twoJ2, i, j, r2[0], r2[1], r2[2], r_mom2[0], r_mom2[1], r_mom2[2], 2) * Hadron::Wigner_D(twoJ2,  j, two_abs_lam2, r_mom2[0], r_mom2[1], r_mom2[2]);
        
        }}
      
      
      for(int i = -twoJ1; i <= twoJ1; i = i+2){
          for(int j = -twoJ1; j <= twoJ1; j = j+2){
              
              D1 += conj(Hadron::Wigner_D(twoJ1, i, two_abs_lam1, r_n_mom1[0], r_n_mom1[1], r_n_mom1[2])) *  Ph::comp_Wigner_d(twoJ1, i, j, r2[0], r2[1], r2[2], r_mom2[0], r_mom2[1], r_mom2[2], 2) * Hadron::Wigner_D(twoJ1,  j, two_abs_lam1, r_mom1[0], r_mom1[1], r_mom1[2]);
              
          }}

      for(int i = -twoJCurr; i <= twoJCurr; i = i+2){
          for(int j = -twoJCurr; j <= twoJCurr; j = j+2){

            DCurr += conj(Hadron::Wigner_D(twoJCurr, i, two_abs_lamCurr, r_n_mom_curr[0], r_n_mom_curr[1], r_n_mom1[2])) *  Ph::comp_Wigner_d(twoJCurr,     i, j, r2[0], r2[1], r2[2], r_mom2[0], r_mom2[1], r_mom2[2], 2) * Hadron::Wigner_D(twoJCurr,  j, two_abs_lamCurr, r_mom_curr[0], r_mom_curr[1], r_mom_curr[2]);

          }}


      //cout << "lam_int" << two_abs_lam1 << two_abs_lam2 << "lam_int" << "\n";

      lambd = std::make_tuple(two_abs_lam1, two_abs_lam2 , two_abs_lamCurr);

      // phase_hel.insert( pair< pair<int,int> , Ph::phChars >(make_pair(two_abs_lam1, two_abs_lam2 ), conj(D1)*D2*R.inverse() ) );

      //Eigen::MatrixXcd ph =  conj(D1)*D2*R.inverse();
      //complex<double> ph =  conj(D1)*D2*conj(DCurr)*DCurr;
      complex<double> ph =  conj(D1)*D2*(DCurr);
      //complex<double> ph = 1;
      //cout << ph << "\n";
      //Eigen::MatrixXcd ph = MatrixXd::Identity(4,4);

      phase_hel.insert( std::make_pair( lambd, ph  ) );
      }}}
    
    out.mom1 = n_mom1;
    out.mom2 = mom2_ref;
    out.lam_phase = phase_hel;
    out.r = R;
    
    return out;


}

int main(){

  int twoJ1 = 0;
  int twoJ2 = 2;
  int twoJCurr = 2;
  Eigen::Vector3d mom1,mom2;

  mom1 << 0,-1,-1;
  mom2 << 0,0,-1;

  Ph::phChars chk;
  chk = Ph::phaseFactor(twoJ1,twoJ2,twoJCurr,mom1,mom2);

  for(map<Ph::tripKey , complex<double> >::iterator  it = chk.lam_phase.begin(); it != chk.lam_phase.end(); it++){
    std::tuple<int,int,int> a = it->first;

    std::cout << std::get<0>(a) << "1mom";
    std::cout << std::get<1>(a) << "2mom";
    std::cout << std::get<2>(a) << "curr";
    cout << it->second << "\n";}


}

//Hadron::Wigner_D(int twoJ, 2*mu, 2*abs_lam, R1_phi, R1_theta, R1_psi);


