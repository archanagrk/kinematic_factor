
int main(){

  int twoJ1 = 0;
  int twoJ2 = 2;
  Eigen::Vector3d mom1,mom2;

  mom1 << 1,0,1;
  mom2 << 1,1,0;

  Ph::phaseFactor(twoJ1,twoJ2,mom1,mom2);

}

//Hadron::Wigner_D(int twoJ, 2*mu, 2*abs_lam, R1_phi, R1_theta, R1_psi);


