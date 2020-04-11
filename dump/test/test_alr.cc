//std lib
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include "math.h"
#include <stdio.h>
#include <Eigen/Dense>
#include <vector>
#include <map>

using namespace std;

double truncate(double num,int precision){
  if(abs(num) < pow(10,-precision)){ num = double(0);}
  return num;
}



bool is_equal(double a, double b)
{
    return std::abs(a - b) < std::numeric_limits<float>::epsilon();
}



map< double,double >  retAng(Eigen::Vector3d ax,double cos_tht,double sin_tht){
  
  ax = ax.cwiseAbs();
  map< double,double >  csa;

  std::sort(ax.data(),ax.data()+ax.size());

  const double PI = (atan(double(1)) * double(4.0));

  if(is_equal(ax(0),0.0)){

    if(is_equal(ax(1),0.0)){
      csa.insert( pair< double , double >(cos_tht, sin_tht));
      return csa;
        }

    else{
         csa.insert( pair< double , double >(cos(PI), sin(PI)));
         csa.insert( pair< double , double >(cos(2*PI), sin(2*PI)));
         return csa;}

  }

  else{
          csa.insert( pair< double , double >(cos(2*PI/3), sin(2*PI/3)));
          csa.insert( pair< double , double >(cos(4*PI/3), sin(4*PI/3)));
          csa.insert( pair< double , double >(cos(2*PI), sin(2*PI)));
          return csa;}

};



/* Rotation matrix between two vectors using Rodrigues formula */


Eigen::MatrixXd rotMat(Eigen::Vector3d ax, double cos_tht, double sin_tht){

  Eigen::MatrixXd rot_mat = Eigen::MatrixXd::Zero(4,4);
  rot_mat(0,0) = 1;
  rot_mat(1,1) = cos_tht+pow(ax(0),2)*(1- cos_tht); rot_mat(1,2) = ax(0)*ax(1)*(1-cos_tht) - ax(2) * sin_tht;
                                                                  rot_mat(1,3) = ax(0)*ax(2)*(1-cos_tht) + ax(1) * sin_tht;
  rot_mat(2,1) = rot_mat(1,2) + 2.0 * ax(2) * sin_tht;  rot_mat(2,2) = cos_tht + pow(ax(1),2)*(1- cos_tht);
                                                                  rot_mat(2,3) = ax(1)*ax(2)*(1-cos_tht) - ax(0) * sin_tht;
  rot_mat(3,1) = rot_mat(1,3) - 2.0 * ax(1) * sin_tht;  rot_mat(3,2) = rot_mat(2,3) + 2.0 * ax(0) * sin_tht;
                                                                  rot_mat(3,3) = cos_tht + pow(ax(2),2)*(1- cos_tht);
  return rot_mat;
}




Eigen::MatrixXd getRotMat(Eigen::Vector3d mom2, Eigen::Vector3d mom1)

{
  double cos_t;

  if(mom1.norm() == 0.0 || mom2.norm() == 0.0 ){cos_t = 0;}
  else{cos_t  = double(mom1.dot(mom2))/double((mom1.norm())*(mom2.norm()));}

  double sin_t;

  if(std::abs(cos_t) < 1 ) {sin_t = sqrt(1- pow(cos_t,2));}
  else{sin_t = sqrt(1- round(pow(cos_t,2)));}

  Eigen::Vector3d perp = mom1.cross(mom2);

  Eigen::Vector3d ax;
  if(!is_equal(perp.norm(),0.0)){ax = perp/double(perp.norm());}
  else{ax << 0,0,0;}


  map< double,double >  csa = retAng(ax,cos_t,sin_t);

  for(map<  double, double >::iterator  it = csa.begin(); it != csa.end(); it++){

  double cos_tht = it->first;
  double sin_tht = it->second;


  Eigen::MatrixXd rot_mat = rotMat(ax, cos_tht, sin_tht);
  Eigen::MatrixXd mom2_f(4,1);
  mom2_f << 0.0, mom2(0), mom2(1), mom2(2);
  Eigen::MatrixXd ans = rot_mat*mom2_f;

  if(is_equal(mom1(0),ans(1)) && is_equal(mom1(1),ans(2)) && is_equal(mom1(2),ans(3))){return rot_mat;}
  else if(is_equal(mom1(0),-ans(1)) && is_equal(mom1(1),-ans(2)) && is_equal(mom1(2),-ans(3))){
    Eigen::MatrixXd  R = Eigen::MatrixXd::Identity(4,4);
    return -R*rot_mat;
  }

  }
  
  Eigen::MatrixXd rot_mat = rotMat(ax, cos_t, sin_t);
  cout << "nope" << "\n";
    //cout << rot_mat << "\n";
    cout << "per" << perp << "\n";
    cout << "cos" <<cos_t << "\n";
    cout << "mom1" << mom1 << "\n";
    cout << "mom2" << mom2 << "\n";
  return rot_mat;
}


int main(){
  Eigen::Vector3d mom1,mom2;
  Eigen::MatrixXd mo1(4,1),mo2(4,1);
    /*for(int i= 0; i <= 2; i++){
      for(int j = 0; j <= 2; j++){                                 // Looping over all mom_1
        for(int k = 0; k <= 2; k++){
            mom1 << i,j,k;
            for(int l= 0; l <= 2; l++){
              for(int m = 0; m <= 2; m++){                                 // Looping over all mom_1
                for(int n = 0; n <= 2; n++){
                  mom2 << l,m,n;
                  //mom1 = mom2;
                  //mom1 = mom1.cwiseAbs();
                  //std::sort(mom1.data(),mom1.data()+mom1.size());

                  mo2 << 0,l,m,n;
  if(mom1.squaredNorm() == mom2.squaredNorm()){
  Eigen::MatrixXd rot(4,4);
  rot = getRotMat(mom1,mom2);
  //cout << rot << "\n";
  //cout << mom1 << "\n";
  //cout << mom2 << "\n";
  //cout << "n" << rot*mo1;
  Eigen::MatrixXd a = rot*mo2;
  //if(!(std::abs(a(1)- mo1(1)) < std::numeric_limits<float>::epsilon()  && std::abs(a(2)- mo1(2)) < std::numeric_limits<float>::epsilon()   && std::abs(a(3)- mo1(3)) < std::numeric_limits<float>::epsilon() )){
  //cout << mo2 << "\n";
  //cout << a << "\n";
  //cout << rot << "\n";
  //cout << mo1 << "\n";
  //cout << "******************" << "\n";
  //}
   //cout << "******************" << "\n";
  }
                }}}}}}*/
    mom1 << 1,1,1;
    mom2 << 1,1,-1;
    mo1 << 0,0,1,1;
    mo2 << 0,1,0,-1;
  Eigen::MatrixXd rot(4,4);
    rot =  getRotMat(mom1,mom2);
    //cout << mo2 << "\n";
    Eigen::MatrixXd a = rot*mo2;
    //cout << a << "\n";
    //cout << rot << "\n";
    //cout << mo1 << "\n";



}


/* Get euler angles from the rotation matrix */

std::vector<double> euAng(Eigen::MatrixXd R){
  std::vector<double> Eangle;

  double R_alpha,R_gamma;
  double R_beta = acos(R(3,3));

  if (sin(R_beta) != 0.00){
    R_alpha = asin((R(2,3)/sin(R_beta)));
    R_gamma = asin((R(2,3)/sin(R_beta)));}

  else {
    R_gamma = asin(R(2,1)) - asin(R(1,2)) ;
    R_alpha = asin(R(2,1)) - R_gamma ;}

  Eangle.push_back(R_alpha), Eangle.push_back(R_beta), Eangle.push_back(R_gamma);

  return Eangle;
 }
