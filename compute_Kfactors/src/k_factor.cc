/* Loop over Lorentz indices to get the coefficient*/


#include "k_factor.h"


complex<double> KFac::KinematicFactor( Eigen::VectorXd& qp, Eigen::VectorXd& qm, Eigen::MatrixXcd& Sub1 , Eigen::MatrixXcd& SubCurr , Eigen::MatrixXcd& Sub2 ){

    Eigen::MatrixXcd sum_sub =  Eigen::MatrixXcd::Zero(4,4);
    sum_sub = Sub1(2,0)*Sub2 * (SubCurr.transpose()).conjugate();
    //cout << sum_sub << "\n";
    complex<double> Coeff = 0;
    for(int i = 0; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
      //cout << "il" << i << l << "il" << "\n";
      //cout << "sum" << (Sub1(2,0))*(SubCurr(i,0))*(Sub2(l,0)) << "sum" << "\n";
      for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
      int e[] = {i,j,k,l};
      //cout << "sum" << (Sub1(2,0))*(SubCurr(i,0))*(Sub2(l,0)) << "sum" << "\n";
      Coeff += LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*sum_sub(l,i);}}}}
      //Coeff += LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*(Sub1(2,0))*(SubCurr(i,0))*(Sub2(l,0));}}}}
      //Coeff += LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*(Sub2(i,0))*(SubCurr(l,0));}}}}
      Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
      return Coeff;
};

Eigen::MatrixXcd KFac::subPhSum( map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2, Ph::phChars& phase ){
    
    Eigen::MatrixXcd sub_phase_sum = Eigen::MatrixXcd::Zero(4,4);
    Eigen::MatrixXcd S1 = Eigen::MatrixXcd::Zero(4,1);
    Eigen::MatrixXcd S2 = Eigen::MatrixXcd::Zero(4,1);
    Eigen::MatrixXcd SCurr = Eigen::MatrixXcd::Zero(4,1);
    Ph::tripKey lambd;
    


    for(map<  int, Eigen::MatrixXcd  >::iterator  it1 = Sub1.begin(); it1 != Sub1.end(); it1++){
        for(map<  int, Eigen::MatrixXcd >::iterator  it2 = Sub2.begin(); it2 != Sub2.end(); it2++){
          for(map<  int, Eigen::MatrixXcd >::iterator  it3 = SubCurr.begin(); it3 != SubCurr.end(); it3++){
            S1 = (it1->second);
            S2 = (it2->second);
            SCurr = (it3->second);
            //lambd = std::make_tuple((it1->first), (it2->first), (it3->first));

            //cout << "lam" <<  (it3->first) << (it2->first) << "lam" << "\n";
            //cout << "phase" << phase.lam_phase[make_pair((it1->first),(it2->first))] << "phase" << "\n";
            
            //sub_phase_sum(i,j) += phase.lam_phase[lambd]*S1(2,0)*S2(i,0)*SCurr(j,0);
            sub_phase_sum += S1(2,0)*S2*(SCurr.transpose()).conjugate();
            //sub_phase_sum(i,j) += S2(i,0)*SCurr(j,0);
            //cout << "i" << i << "\n";
            //cout << "sub" <<  sub_phase_sum.col(i) << "sub" << "\n";
            }
            
        }
    }

 return sub_phase_sum;   
 }



complex<double> KFac::KinematicFactorwithPhase(Eigen::VectorXd& qp, Eigen::VectorXd& qm, map< int, Eigen::MatrixXcd >& Sub1 , map< int, Eigen::MatrixXcd >& SubCurr , map< int, Eigen::MatrixXcd >& Sub2, Ph::phChars& phase ){
    
    complex<double> Coeff = 0;
    Eigen::MatrixXcd sub_phase_sum = KFac::subPhSum(Sub1, SubCurr, Sub2, phase);
    //cout << "sub" <<  sub_phase_sum << "sub" << "\n";
    for(int i = 0; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
      //cout << "il" << i << l << "il" << "\n";
      //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
      for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){
        int e[] = {i,j,k,l};
        //cout << "sub" <<  sub_phase_sum(i,l) << "sub" << "\n";
        Coeff += LevCiv::LeviCivita(e,4)*(qp(j,0))*(qm(k,0))*sub_phase_sum(l,i);
    }}
    //cout << "cf" << Coeff << "cf" << "\n";
    }}
    Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
    return Coeff;
};
