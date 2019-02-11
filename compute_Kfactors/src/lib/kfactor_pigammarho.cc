#include "kfactor_pigammarho.h"

//*****************************************************************************************************************
complex<double> KfacSSV::operator()( const KFacParams& params ) const {
    
  MatrixXcd sum_sub =  MatrixXcd::Zero(4,4);
  sum_sub = params.subPhSum();

  complex<double> Coeff = 0;
  int i =0;
  
  for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      Coeff += LevCiv::LeviCivita(e,4) * (params.qp(j,0)) * (params.qm(k,0)) * sum_sub(l,0);
    }}
  }

  Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
  return Coeff;
};
//*****************************************************************************************************************

//*****************************************************************************************************************
complex<double> KfacSVV::operator()( const KFacParams& params ) const {

  complex<double> Coeff = 0;
  MatrixXcd sub_sum = params.subPhSum();

  for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      Coeff += LevCiv::LeviCivita(e,4) * (params.qp(j,0)) * (params.qm(k,0)) * sub_sum(l,i);

    }}
  }}

  Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
  return Coeff;
      
};
//*****************************************************************************************************************

//*****************************************************************************************************************
complex<double> KfacSSVwPhCorr::operator()( const KFacParams& params ) const {
  
  complex<double> Coeff = 0;
  Eigen::MatrixXcd sub_phase_sum = params.subPhSum();
  int i = 0;

  for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      Coeff += LevCiv::LeviCivita(e,4) * (params.qp(j,0)) * (params.qm(k,0)) * sub_phase_sum(l,i);

    }}
  }

  Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
  return Coeff;

};
//*****************************************************************************************************************

//*****************************************************************************************************************
complex<double> KfacSVVwPhCorr::operator()( const KFacParams& params ) const {

  complex<double> Coeff = 0;
  Eigen::MatrixXcd sub_phase_sum = params.subPhSum();


  for(int i = 1; i < 4; i++ ){for(int l = 0; l < 4; l++ ){
    for(int k = 0; k < 4; k++ ){for(int j = 0; j< 4; j++ ){

      int e[] = {i,j,k,l};
      Coeff += LevCiv::LeviCivita(e,4) * (params.qp(j,0)) * (params.qm(k,0)) * sub_phase_sum(l,i);

    }}
  }}

  Coeff = {KfUt::truncate(Coeff.real(),10),KfUt::truncate(Coeff.imag(),10)};
  return Coeff;

};
//*****************************************************************************************************************
