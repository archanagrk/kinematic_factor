#include "kfac_params.h"



  //**********************************************************************************************************************

    /* The kfactor params for different J's */

  //**********************************************************************************************************************

/* A fucntion that returns a three-dimensional matrix of ploarization matrix at each npt dotted with 
    thre respexctive subduction coeffs. A useful quantity to later compute the kinematic factors */

vector<MatrixXcd> KFacParams::subPhSum() const{

  /* Inialize variables */
  vector<MatrixXcd> sub_phase_sum;

  MatrixXcd S1 = MatrixXcd::Zero(4,1);
  MatrixXcd S2 = MatrixXcd::Zero(4,1);
  MatrixXcd SCurr = MatrixXcd::Zero(4,1);
  Ph::tripKey lambd;  

  for(int i = 0; i < 4;i++){
    MatrixXcd sub_ph_sum = MatrixXcd::Zero(4,4);

    for(map<  int, MatrixXcd  >::const_iterator  it1 = Sub1.begin(); it1 != Sub1.end(); it1++){
      for(map<  int, MatrixXcd >::const_iterator  it2 = Sub3.begin(); it2 != Sub3.end(); it2++){
        for(map<  int, MatrixXcd >::const_iterator  it3 = SubCurr.begin(); it3 != SubCurr.end(); it3++){

          S1 = (it1->second);
          S2 = (it2->second);
          SCurr = (it3->second); SCurr = KfUt::Gmunu() * SCurr;
          lambd = std::make_tuple((it1->first), (it2->first), (it3->first));
          sub_ph_sum += phase.lam_phase.find(lambd)->second * S1.conjugate()(i,0) * S2 * (((phase.r).inverse().transpose() * SCurr).transpose()); // dont conjugate as S is real but the pol is not transposed as it converts to cirecular basis

          
        }   
      }
    }

    sub_phase_sum.push_back(sub_ph_sum);
  }

  //cout << sub_phase_sum << endl;
  return sub_phase_sum;
};

//**********************************************************************************************************************

  /* A fucntion that returns the helicity at each npt */

Ph::tripKey KFacParams::two_abs_lam() const{

  Ph::tripKey twice_abs_lam; 
  
  for(map<  int, MatrixXcd  >::const_iterator  it1 = Sub1.begin(); it1 != Sub1.end(); it1++){
    for(map<  int, MatrixXcd >::const_iterator  it2 = Sub3.begin(); it2 != Sub3.end(); it2++){
      for(map<  int, MatrixXcd >::const_iterator  it3 = SubCurr.begin(); it3 != SubCurr.end(); it3++){
         
        twice_abs_lam = std::make_tuple(abs(it1->first),abs(it2->first),abs(it3->first));

      }
    }
  }

  return twice_abs_lam;
};


//**********************************************************************************************************************
//Constructor


KFacParams::KFacParams(map< int, Eigen::MatrixXcd > Sub1_, map< int, Eigen::MatrixXcd > SubCurr_, map< int, Eigen::MatrixXcd > Sub3_, Ph::phChars phase_, VectorXd  qp_, VectorXd  qm_)
{
  Sub1 = Sub1_;
  SubCurr = SubCurr_;
  Sub3 = Sub3_;
  phase = phase_;
  qp = qp_;
  qm = qm_;
};

//**********************************************************************************************************************




