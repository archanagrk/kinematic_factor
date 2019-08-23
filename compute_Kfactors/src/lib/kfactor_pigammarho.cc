#include "kfactor_pigammarho.h"

//*****************************************************************************************************************
complex<double> KfacSSV::operator()( const KFacParams& params ) const {
    
  MatrixXcd sum_sub =  MatrixXcd::Zero(4,4);
  sum_sub = params.subPhSum();

  VectorXd q_f = params.qp + params.qm;
  double m_f_sq  = (pow(q_f[0],2) - pow(q_f[1],2) - pow(q_f[2],2) - pow(q_f[3],2) )/4.0;

  complex<double> norm = 1/(sqrt(m_f_sq)); //much more sensible norm

  complex<double> Coeff = 0;
  int i =0;
  
  for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      Coeff += LevCiv::LeviCivita(e,4) * (params.qp(k,0)) * (params.qm(l,0)) * sum_sub(j,0) ;
    }}
  }

  Coeff *= norm; 

  Coeff = {KfUt::truncate(Coeff.real(),5),KfUt::truncate(Coeff.imag(),5)};
  return Coeff;
};
//*****************************************************************************************************************

//*****************************************************************************************************************
complex<double> KfacSVV::operator()( const KFacParams& params ) const {

  complex<double> Coeff = 0;
  MatrixXcd sub_sum = params.subPhSum();

  VectorXd q_f = params.qp + params.qm;
  double m_f_sq  = (pow(q_f[0],2) - pow(q_f[1],2) - pow(q_f[2],2) - pow(q_f[3],2) )/4.0;

  complex<double> norm = 1/(sqrt(m_f_sq)); //much more sensible norm

  for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      Coeff += LevCiv::LeviCivita(e,4) * (params.qp(k,0)) * (params.qm(l,0)) * sub_sum(j,i);

    }}
  }}

  Coeff *= norm; 

  Coeff = {KfUt::truncate(Coeff.real(),5),KfUt::truncate(Coeff.imag(),5)};
  return Coeff;
      
};

//*****************************************************************************************************************
