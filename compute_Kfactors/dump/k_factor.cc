/* Loop over Lorentz indices to get the coefficient*/


#include "k_factor.h"

bool KFac::nearlyequal(complex<double>& a, complex<double>& b) {
  double diff = abs(a - b);
  double mag = (abs(a) + abs(b))/2;
  return diff <= (mag * DBL_EPSILON * (1ull << ULP_N));
}


Eigen::MatrixXcd KFac::subPhSum( map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2, Ph::phChars& phase ){
  
  Eigen::MatrixXcd sub_phase_sum = Eigen::MatrixXcd::Zero(4,4);
  Eigen::MatrixXcd S1 = Eigen::MatrixXcd::Zero(4,1);
  Eigen::MatrixXcd S2 = Eigen::MatrixXcd::Zero(4,1);
  Eigen::MatrixXcd SCurr = Eigen::MatrixXcd::Zero(4,1);
  Ph::tripKey lambd;
  
  
  
  for(map<  int, Eigen::MatrixXcd  >::iterator  it1 = Sub1.begin(); it1 != Sub1.end(); it1++){
    for(map<  int, Eigen::MatrixXcd >::iterator  it2 = Sub2.begin(); it2 != Sub2.end(); it2++){
      for(map<  int, Eigen::MatrixXcd >::iterator  it3 = SubCurr.begin(); it3 != SubCurr.end(); it3++){
        if((it1->first)+(it2->first)+(it3->first)){continue;}
        else{
          S1 = (it1->second);
          S2 = (it2->second);
          SCurr = (it3->second);
          lambd = std::make_tuple((it1->first), (it2->first), (it3->first));
          
          //cout << "lam" <<  (it3->first)  << "lam" << "\n";
          //cout << "phase" << phase.lam_phase[lambd] << "phase" << "\n";
          
          //sub_phase_sum(i,j) += phase.lam_phase[lambd]*S1(2,0)*S2(i,0)*SCurr(j,0);
          sub_phase_sum += phase.lam_phase[lambd] * S1(2,0) * S2 * (SCurr.transpose()).conjugate();
          //cout <<  phase.lam_phase[lambd] << "Phase" << "\n";
          //sub_phase_sum += S1(2,0) * S2 * SCurr.transpose();
          //sub_phase_sum(i,j) += S2(i,0)*SCurr(j,0);
          //cout << "sub" <<  sub_phase_sum.col(i) << "sub" << "\n";
        }
        
      }
      
    }
  }
  
  return sub_phase_sum;
}

Eigen::MatrixXcd KFac::subSum( map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2 ){
  
  Eigen::MatrixXcd sub_sum = Eigen::MatrixXcd::Zero(4,4);
  Eigen::MatrixXcd S1 = Eigen::MatrixXcd::Zero(4,1);
  Eigen::MatrixXcd S2 = Eigen::MatrixXcd::Zero(4,1);
  Eigen::MatrixXcd SCurr = Eigen::MatrixXcd::Zero(4,1);
  Ph::tripKey lambd;
  
  
  
  for(map<  int, Eigen::MatrixXcd  >::iterator  it1 = Sub1.begin(); it1 != Sub1.end(); it1++){
    for(map<  int, Eigen::MatrixXcd >::iterator  it2 = Sub2.begin(); it2 != Sub2.end(); it2++){
      for(map<  int, Eigen::MatrixXcd >::iterator  it3 = SubCurr.begin(); it3 != SubCurr.end(); it3++){
        if((it1->first)+(it2->first)+(it3->first)){continue;}
        else{
          S1 = (it1->second);
          S2 = (it2->second);
          SCurr = (it3->second);
          

          sub_sum += S1(2,0) * S2 * (SCurr.transpose()).conjugate();

        }

        
      }
      
    }
  }
  
  return sub_sum;
}


complex<double> KFac::KFacold( Eigen::VectorXd& qp, Eigen::VectorXd& qm, Eigen::MatrixXcd& Sub1 , Eigen::MatrixXcd& SubCurr , Eigen::MatrixXcd& Sub2 ){

   Eigen::MatrixXcd sum_sub =  Eigen::MatrixXcd::Zero(4,4);
   sum_sub = Sub1(2,0) * Sub2 * (SubCurr.transpose()).conjugate();
   //cout << sum_sub << "\n";
   complex<double> Coeff = 0;
   for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
       //cout << "il" << i << l << "il" << "\n";
       //cout << "sum" << (Sub1(2,0))*(SubCurr(i,0))*(Sub2(l,0)) << "sum" << "\n";
       for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
           int e[] = {i,j,k,l};
           //cout << "sum" << round(qp(j,0)) * round(qm(k,0))  << "sum" << "\n";
           Coeff += LevCiv::LeviCivita(e,4) * (qp(j,0)) * (qm(k,0)) * sum_sub(l,i);}}}}
           //Coeff += LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*(Sub1(2,0))*(SubCurr(i,0))*(Sub2(l,0));}}}}
           //Coeff += LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*(Sub2(i,0))*(SubCurr(l,0));}}}}
   Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
   return Coeff;
 };

complex<double> KFac::KFacSVV(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2 ){
  
  complex<double> Coeff = 0;
  Eigen::MatrixXcd sub_sum = KFac::subSum(Sub1, SubCurr, Sub2);
  //cout << "sub" <<  sub_phase_sum << "sub" << "\n";
  for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
    //cout << "il" << i << l << "il" << "\n";
    //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
      int e[] = {i,j,k,l};
      //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
      //cout << LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*sub_phase_sum(l,i) <<  "\n";
      Coeff += LevCiv::LeviCivita(e,4) * (qp(j,0)) * (qm(k,0)) * sub_sum(l,i);
    }}
    //cout << "cf" << Coeff << "cf" << "\n";
  }}
  Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
  return Coeff;
};





complex<double> KFac::KFacSSV( Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2 ){
    
    Eigen::MatrixXcd sum_sub =  Eigen::MatrixXcd::Zero(4,4);
    sum_sub = KFac::subSum(Sub1, SubCurr, Sub2);
    //cout << sum_sub << "\n";
    complex<double> Coeff = 0;
    int i =0;
    
    for(int l = 0; l < 4; l++ ){
        //cout << "il" << i << l << "il" << "\n";
        //cout << "sum" << (Sub1(2,0))*(SubCurr(i,0))*(Sub2(l,0)) << "sum" << "\n";
        for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
            int e[] = {i,j,k,l};
            //cout << "sum" << round(qp(j,0)) * round(qm(k,0))  << "sum" << "\n";
            Coeff += LevCiv::LeviCivita(e,4) * (qp(j,0)) * (qm(k,0)) * sum_sub(l,2);}}}
            //Coeff += LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*(Sub1(2,0))*(SubCurr(i,0))*(Sub2(l,0));}}}}
            //Coeff += LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*(Sub2(i,0))*(SubCurr(l,0));}}}}
    Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
    return Coeff;
};




complex<double> KFac::KinematicFactorwithPhase_j0(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2, Ph::phChars& phase ){
    
    complex<double> Coeff = 0;
    Eigen::MatrixXcd sub_phase_sum = KFac::subPhSum(Sub1, SubCurr, Sub2, phase);
    int i = 0;
    //cout << "sub" <<  sub_phase_sum << "sub" << "\n";
    for(int l = 0; l < 4; l++ ){
        //cout << "il" << i << l << "il" << "\n";
        //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
        for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
            int e[] = {i,j,k,l};
            //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
            //cout << LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*sub_phase_sum(l,i) <<  "\n";
            Coeff += LevCiv::LeviCivita(e,4) * (qp(j,0)) * (qm(k,0)) * sub_phase_sum(l,i);
        }}
        //cout << "cf" << Coeff << "cf" << "\n";
    }
    Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
    return Coeff;
};




complex<double> KFac::KinematicFactorwithPhase(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2, Ph::phChars& phase ){
    
    complex<double> Coeff = 0;
    Eigen::MatrixXcd sub_phase_sum = KFac::subPhSum(Sub1, SubCurr, Sub2, phase);
    //cout << "sub" <<  sub_phase_sum << "sub" << "\n";
    for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
      //cout << "il" << i << l << "il" << "\n";
      //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
      for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
        int e[] = {i,j,k,l};
        //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
        //cout << LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*sub_phase_sum(l,i) <<  "\n"; 
        Coeff += LevCiv::LeviCivita(e,4) * (qp(j,0)) * (qm(k,0)) * sub_phase_sum(l,i);
    }}
    //cout << "cf" << Coeff << "cf" << "\n";
    }}
    Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
    return Coeff;
};

complex<double> KFac::KFacwoS(Eigen::VectorXd& qp, Eigen::VectorXd& qm, std::vector<double>& r3, double& mom3_sq, const int& two_lam3, double& m3_sq, std::vector<double>& r_curr, double& mom_curr_sq, const int& two_lam_curr, double& m_curr_sq){
  Eigen::MatrixXcd pol;
  Eigen::MatrixXcd pol_curr;

 


  pol = KfUt::Gmunu()*PolVec::getPol4(mom3_sq, two_lam3, m3_sq, r3[0], r3[1], r3[2]);
  pol_curr = (KfUt::Gmunu()*PolVec::getPol4(mom_curr_sq, two_lam_curr, m_curr_sq, r_curr[0], r_curr[1], r_curr[2])).conjugate();
  

  complex<double> Coeff = 0;
  
  //cout << "pol" << pol.transpose() << "\n";
  //cout << "pol_c" << pol_curr.transpose() << "\n";
  //cout << qm.transpose() << "qm" << "\n";
  //cout << qp.transpose() << "qp" << "\n";

  for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
      int e[] = {i,j,k,l};
      //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
      //cout << LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*sub_phase_sum(l,i) <<  "\n";
    Coeff += LevCiv::LeviCivita(e,4) * (qp(j,0)) * (qm(k,0)) * pol(l,0) * pol_curr(i,0);
    //cout <<  mom_curr_sq << two_lam_curr << "\n";
    //cout << "check" << Coeff << "\n"; 

    }}}}
  
  Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
  return Coeff;
};

complex<double> KFac::KFacwoSwP(Eigen::VectorXd& qp, Eigen::VectorXd& qm, const int& two_lam1, std::vector<double>& r3, double& mom3_sq, const int& two_lam3, double& m3_sq, std::vector<double>& r_curr, double& mom_curr_sq, const int& two_lam_curr, double& m_curr_sq, Ph::phChars& phase){
  
  Eigen::MatrixXcd pol;
  Eigen::MatrixXcd pol_curr;
  

  
  pol =  KfUt::Gmunu()*PolVec::getPol4(mom3_sq, two_lam3, m3_sq, r3[0], r3[1], r3[2]);
  pol_curr = (KfUt::Gmunu()*PolVec::getPol4(mom_curr_sq, two_lam_curr, m_curr_sq, r_curr[0], r_curr[1], r_curr[2])).conjugate();
  //cout <<  pol.transpose() << "pol" << "\n";
  //cout << pol_curr.transpose() << "pol_c" << "\n";
  //cout << qm.transpose() << "qm" << "\n";
  //cout << qp.transpose() << "qp" << "\n";
  complex<double> Coeff = 0;
  
  Ph::tripKey lambd;
  lambd = std::make_tuple(two_lam1, two_lam_curr , two_lam3);
  //cout << "Phase" << phase.lam_phase[lambd] << "\n";
  //cout << "Abs" << abs(phase.lam_phase[lambd]) << "\n";
  
  
  
  for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
    int e[] = {i,j,k,l};
    //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
    //cout << LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*sub_phase_sum(l,i) <<  "\n";
    Coeff += phase.lam_phase[lambd] * LevCiv::LeviCivita(e,4) * (qp(j,0)) * (qm(k,0)) * pol(l,0) * pol_curr(i,0);
    //cout << "check" << Coeff << "\n";    
    
  }}}}
  //cout << "Phase is:" << phase.lam_phase[lambd] << "\n";
  
  //cout << "Phase is:" << std::abs(phase.lam_phase[lambd]) << "\n";
  Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
  return Coeff;
};
