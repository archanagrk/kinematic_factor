#include "kfactor_pigammarho.h"
 //**********************************************************************************************************************
 //**********************************************************************************************************************
/*
  the form of kinematic factors for the transition PS -> gamma + V

  There is one decomposition for the temporal part of photon current (KfacVSS / KfacSSV) /mu = 0
  and one for the spatial (KfacVVS / KfacSVV) /mu = 1,2,3
*/

    /* 
    coeff = i * 1/m_PS * \epsilon(q,\lamda_gamma)_\mu * [\levichivita^{\mu\nu\rho\sigma}\epsilon^{*}(p_f,\lambda_f)_\nu (p_i + p_f)_\rho + (p_f - p_i)_\sigma]*Subductions]
    */

//*****************************************************************************************************************
 //**********************************************************************************************************************

vector<complex<double>> KfacSSV::operator()( const KFacParams& params ) const {

  vector<complex<double>> Coeff;
    
  MatrixXcd sum_sub =  MatrixXcd::Zero(4,4);
  sum_sub = params.subPhSum().at(2);

  VectorXd q_f = params.qp + params.qm;
  double m_f_sq  = (pow(q_f[0],2) - pow(q_f[1],2) - pow(q_f[2],2) - pow(q_f[3],2) )/4.0;

  complex<double> norm = 1/(sqrt(m_f_sq)); //much more sensible norm

  complex<double> tmp = 0;
  int i =0;
  
  for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      tmp += LevCiv::LeviCivita(e,4) * (params.qp(k,0)) * (params.qm(l,0)) * sum_sub(j,0) ;
    }}
  }

  tmp *= norm; 

  tmp = {KfUt::truncate(tmp.real(),5),KfUt::truncate(tmp.imag(),5)};

  Coeff.push_back(tmp);
  return Coeff;
};
//*****************************************************************************************************************

//*****************************************************************************************************************
vector<complex<double>> KfacSVV::operator()( const KFacParams& params ) const {

  vector<complex<double>> Coeff;

  complex<double> tmp = 0;
  MatrixXcd sub_sum = params.subPhSum().at(2);

  complex<double> z_i(0.,1.);
    
  VectorXd q_f = params.qp + params.qm;
  double m_f_sq  = (pow(q_f[0],2) - pow(q_f[1],2) - pow(q_f[2],2) - pow(q_f[3],2) )/4.0;

  complex<double> norm = 1/(sqrt(m_f_sq)); //much more sensible norm

  for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      tmp += z_i * LevCiv::LeviCivita(e,4) * (params.qp(k,0)) * (params.qm(l,0)) * sub_sum(j,i);

    }}
  }}

  tmp *= norm; 

  tmp = {KfUt::truncate(tmp.real(),5),KfUt::truncate(tmp.imag(),5)};

  Coeff.push_back(tmp);
  return Coeff;
      
};

//*****************************************************************************************************************
 //**********************************************************************************************************************
/*
  the form of kinematic factors for the transition V + gamma -> PS

  There is one decomposition for the temporal part of photon current (KfacVSS) /mu = 0
  and one for the spatial (KfacVVS) /mu = 1,2,3
*/

    /* 
    coeff = i * 1/m_PS * \epsilon(q,\lamda_gamma)_\mu * [\levichivita^{\mu\nu\rho\sigma}\epsilon(p_i,\lambda_i)_\nu (p_i + p_f)_\rho + (p_f - p_i)_\sigma * Subductions]
    */

//*****************************************************************************************************************
//*****************************************************************************************************************
vector<complex<double>> KfacVSS::operator()( const KFacParams& params ) const {

  vector<complex<double>> Coeff;


  vector<MatrixXcd> sum_sub = params.subPhSum();

  VectorXd q_i = params.qp - params.qm;
  double m_i_sq  = (pow(q_i[0],2) - pow(q_i[1],2) - pow(q_i[2],2) - pow(q_i[3],2) )/4.0;

  complex<double> norm = 1/(sqrt(m_i_sq)); //much more sensible norm

  complex<double> tmp = 0;
  int i =0;
  
  for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      tmp += LevCiv::LeviCivita(e,4) * (params.qp(k,0)) * (params.qm(l,0)) * sum_sub.at(j)(2,0) ;
    }}
  }

  tmp *= norm; 

  tmp = {KfUt::truncate(tmp.real(),5),KfUt::truncate(tmp.imag(),5)};

  Coeff.push_back(tmp);
  return Coeff;
};
//*****************************************************************************************************************

//*****************************************************************************************************************
vector<complex<double>> KfacVVS::operator()( const KFacParams& params ) const {

  vector<complex<double>> Coeff;

  complex<double> tmp = 0;
  vector<MatrixXcd> sub_sum = params.subPhSum();

  complex<double> z_i(0.,1.);
    
  VectorXd q_i = params.qp - params.qm;
  double m_i_sq  = (pow(q_i[0],2) - pow(q_i[1],2) - pow(q_i[2],2) - pow(q_i[3],2) )/4.0;

  complex<double> norm = 1/(sqrt(m_i_sq)); //much more sensible norm

  for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      tmp += z_i * LevCiv::LeviCivita(e,4) * (params.qp(k,0)) * (params.qm(l,0)) * sub_sum.at(j)(2,i);

    }}
  }}

  tmp *= norm; 

  tmp = {KfUt::truncate(tmp.real(),5),KfUt::truncate(tmp.imag(),5)};

  Coeff.push_back(tmp);
  return Coeff;
      
};

//*****************************************************************************************************************
