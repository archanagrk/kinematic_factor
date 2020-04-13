#include "kfactor_pigammapi.h"

//**********************************************************************************************************************

  /* function to check equality */

//**********************************************************************************************************************

namespace{
  bool is_equal(double a, double b){
    return std::abs(a - b) < pow(10,-5);
  }
}

//*****************************************************************************************************************
vector<complex<double>> KfacSSS::operator()( const KFacParams& params ) const {

  vector<complex<double>> Coeff;
    
  MatrixXcd sum_sub =  MatrixXcd::Zero(4,4);
  sum_sub = params.subPhSum().at(2);
  complex<double> factor;

  double Q_sq = - params.qm.squaredNorm();
  VectorXd q_f = params.qp + params.qm;
  VectorXd q_i = params.qp - params.qm;
  double m_f_sq  = (pow(q_f[0],2) - pow(q_f[1],2) - pow(q_f[2],2) - pow(q_f[3],2) )/4.0;
  double m_i_sq  = (pow(q_i[0],2) - pow(q_i[1],2) - pow(q_i[2],2) - pow(q_i[3],2) )/4.0;


  if(is_equal(m_i_sq,m_f_sq) ){factor = {0.0,0.0};}
  else{factor = {((m_f_sq - m_i_sq)/Q_sq),0.0}; }

  complex<double> tmp = 0;
  
  tmp += ( (params.qp(0,0)) + (  factor * (params.qm(0,0)) ) ) * sum_sub(2,0);

  tmp = {KfUt::truncate(tmp.real(),5),KfUt::truncate(tmp.imag(),5)};

  Coeff.push_back(tmp);
  return Coeff;
};
//*****************************************************************************************************************

//*****************************************************************************************************************
vector<complex<double>> KfacSVS::operator()( const KFacParams& params ) const {

  vector<complex<double>> Coeff;

  complex<double> tmp = 0;
  MatrixXcd sub_sum = params.subPhSum().at(2);
  complex<double> factor;

  complex<double> one(-1.,0.); // to convert four vector upstairs to downstairs

  double Q_sq = - params.qm.squaredNorm();
  VectorXd q_f = params.qp + params.qm;
  VectorXd q_i = params.qp - params.qm;
  double m_f_sq  = (pow(q_f[0],2) - pow(q_f[1],2) - pow(q_f[2],2) - pow(q_f[3],2) )/4.0;
  double m_i_sq  = (pow(q_i[0],2) - pow(q_i[1],2) - pow(q_i[2],2) - pow(q_i[3],2) )/4.0;


  if(is_equal(m_i_sq,m_f_sq) ){factor = {0.0,0.0};}
  else{factor = {((m_f_sq - m_i_sq)/Q_sq),0.0};}

  for(int i = 1; i < 4; i++ ){
    tmp += one * (( (params.qp(i,0)) + (  factor * (params.qm(i,0)) ) ) * sub_sum(2,i));
  }

  tmp = {KfUt::truncate(tmp.real(),5),KfUt::truncate(tmp.imag(),5)};

  Coeff.push_back(tmp);
  return Coeff;
      
};

//*****************************************************************************************************************
