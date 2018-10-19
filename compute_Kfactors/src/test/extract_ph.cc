/* main code to return the irreps to compute */

#include "kfactors.h"

int main(int argc, char** argv){
  if( argc != 12 ){
    cerr << "get_irreps <twoJ1> <P1> <m1sq> <twoJ_curr> <P_curr>  <twoJ3> <P3> <m3sq> <max_mom1> <max_mom2> <max_mom3> \n ";
    exit(1); }

  int two_J1;  {istringstream a(argv[1]); a >> two_J1;};
  int P1;   {istringstream a(argv[2]); a >> P1;};
  double m1_sq;   {istringstream a(argv[3]); a >> m1_sq;};
  int two_J2;  {istringstream a(argv[4]); a >> two_J2;};
  int P2;   {istringstream a(argv[5]); a >> P2;};
  int two_J3;  {istringstream a(argv[6]); a >> two_J3;};
  int P3;   {istringstream a(argv[7]); a >> P3;};
  double m3_sq;   {istringstream a(argv[8]); a >> m3_sq;};
  int max_mom1;   {istringstream a(argv[9]); a >> max_mom1;};
  int max_mom2;   {istringstream a(argv[10]); a >> max_mom2;};
  int max_mom3;   {istringstream a(argv[11]); a >> max_mom3;};

  complex<double> Coeff;

  if((two_J1==0) && (P1==-1) && (two_J2==2) && (P2==-1) && (two_J3==2) && (P3==-1)){

    Vector3d mom1(3,1);    Vector3d mom3(3,1);    Vector3d mom_curr(3,1);

    for(int i=-int(sqrt(max_mom1)); i <= int(sqrt(max_mom1)); i++){
      for(int j = -int(sqrt(max_mom1)); j <= int(sqrt(max_mom1)); j++){                                 // Looping over all mom_1
        for(int k = -int(sqrt(max_mom1)); k <= int(sqrt(max_mom1)); k++){


          double mom1_sq = (pow(i,2)+pow(j,2)+pow(k,2));

          if(max_mom1 >= mom1_sq){

            for(int l = -int(sqrt(max_mom3)); l <= int(sqrt(max_mom3)); l++){
               for(int m=-int(sqrt(max_mom3)); m <= int(sqrt(max_mom3)) ; m++){                         // Looping over all mom_3
                 for(int n=-int(sqrt(max_mom3)); n<=int(sqrt(max_mom3)); n++){

                  double mom3_sq = (pow(l,2)+pow(m,2)+pow(n,2));

                  if(max_mom3 >= mom3_sq){

                   mom1 << i,j,k;  mom3 << l,m,n;  mom_curr << (l-i),(m-j),(n-k);


                   double mom_curr_sq = (pow((i-l),2)+pow((j-m),2)+pow((k-n),2));

                   if(max_mom2 >= mom_curr_sq){

                    VectorXd  qp(4,1);  VectorXd  qm(4,1);  // q+ = (p1-p2) q- = (p1+p2)
                       
                    Ph::phChars phase;
                    
                    phase = Ph::phaseFactor(two_J1, two_J3, two_J2, mom1, mom3);
                    
                    double n_mom1_sq = phase.mom1.squaredNorm();
                    double n_mom3_sq = phase.mom2.squaredNorm();
                       
                    Eigen::Vector3d n_mom_curr;
                    n_mom_curr << (phase.mom2(0)-phase.mom1(0)),(phase.mom2(1)-phase.mom1(1)),(phase.mom2(2)-phase.mom1(2));

                    double n_mom_curr_sq = n_mom_curr.squaredNorm(); 

                    qp  << (sqrt(m1_sq+n_mom1_sq)+sqrt(m3_sq+n_mom3_sq)),-(phase.mom1(0)+phase.mom2(0)),-(phase.mom1(1)+phase.mom2(1)),-(phase.mom1(2)+phase.mom2(2));
                    qm  << (sqrt(m3_sq+n_mom3_sq)-sqrt(m1_sq+n_mom1_sq)),-(phase.mom2(0)-phase.mom1(0)),-(phase.mom2(1)-phase.mom1(1)),-(phase.mom2(2)-phase.mom1(2));

                    //qp  << (sqrt(m1_sq+mom1_sq)+sqrt(m3_sq+mom3_sq)),(mom1(0)+mom3(0)),(mom1(1)+mom3(1)),(mom1(2)+mom3(2));
                    //qm  << (sqrt(m1_sq+mom1_sq)-sqrt(m3_sq+mom3_sq)),(mom3(0)-mom1(0)),(mom3(1)-mom1(1)),(mom3(2)-mom1(2));

                    //cout << qp.transpose() << "qp" << "\n";
                    //cout << qm.transpose() << "qm" << "\n";

                    double n_m_curr_sq =  pow(sqrt(n_mom3_sq + m3_sq)-sqrt(n_mom1_sq + m1_sq) ,2) - n_mom_curr_sq;
                    double m_curr_sq =  pow(sqrt(mom3_sq + m3_sq)-sqrt(mom1_sq + m1_sq) ,2) - mom_curr_sq;

                    for(int two_lam1 = -two_J1; two_lam1 <= two_J1; two_lam1 = two_lam1 +2){
                      for(int two_lam2 = -two_J3; two_lam2 <= two_J3; two_lam2 = two_lam2 +2){
                        for(int two_lamCurr = -two_J2; two_lamCurr <= two_J2; two_lamCurr = two_lamCurr +2){

                          Ph::tripKey lambd;
                          lambd  = std::make_tuple(two_lam1, two_lam2 , two_lamCurr);
                          
                          cout << mom1.transpose() << "mom_i" << "\n";
                          cout << mom3.transpose() << "mom_f" << "\n";
                          
                          complex<double> ph = {KfUt::truncate(phase.lam_phase[lambd].real(),10),KfUt::truncate(phase.lam_phase[lambd].imag(),10)};

                          cout << ph  << "phase" << "\n";
                    

                   }}}
                  }
                 }
             }}}
           }
         }}}
  }
};


