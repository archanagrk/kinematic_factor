int main(){
   const double PI = (atan(double(1)) * double(4.0));
  double alpha = PI/4;
  double beta = acos(1/sqrt(3));
  double gamma = 0;

  Eigen::MatrixXd tes = Rot::eulerRotMat(alpha, beta, gamma);
  cout << tes;
}
