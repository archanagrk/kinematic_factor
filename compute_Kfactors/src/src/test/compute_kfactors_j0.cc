/* main code to return the irreps to compute */

#include "lib/kfactors.h"

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

  if((two_J1==0) && (P1==-1) && (two_J2==0) && (P2==-1) && (two_J3==2) && (P3==-1)){

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

                    qp  << (sqrt(m1_sq+mom1_sq)+sqrt(m3_sq+mom3_sq)),-(i+l),-(j+m),-(k+n);
                    qm  << (sqrt(m3_sq+mom3_sq)-sqrt(m1_sq+mom1_sq)),-(l-i),-(m-j),-(n-k);

                    double m_curr_sq =  pow(sqrt(mom3_sq + m3_sq)-sqrt(mom1_sq + m1_sq) ,2) - mom_curr_sq;

                    string LG1 = generateLittleGroup(mom1);
                    string LG3 = generateLittleGroup(mom3);
                    string LG_curr = generateLittleGroup(mom_curr);



		                std::vector<double> r1 = refAngles(mom1);
		                std::vector<double> r_curr = refAngles(mom_curr);
		                std::vector<double> r3 = refAngles(mom3);                       
                       

                    std::vector<std::string> irrep1 = getIrrep(two_J1,P1,LG1);
                    std::vector<std::string> irrep_curr = getIrrep(two_J2,P2,LG_curr);
                    std::vector<std::string> irrep3 = getIrrep(two_J3,P3,LG3);

                    irrep_label rep1; irrep_label rep_curr; irrep_label rep3;


                    rep1.twoJ = two_J1; rep3.twoJ = two_J3; rep_curr.twoJ = two_J2;
                    rep1.P = P1; rep3.P = P3; rep_curr.P = P2;


                    for(auto p = irrep1.begin(); p != irrep1.end(); p++){
                        for(auto q = irrep3.begin(); q != irrep3.end(); q++){         // Looping over all irreps at the source, sink & isertion
                           for(auto r = irrep_curr.begin(); r != irrep_curr.end(); r++){

                                rep1.irrep = *p; rep3.irrep = *q; rep_curr.irrep = *r;


                                rep1.n = find_n_subduced_embeddings(LG1, rep1.irrep, rep1.twoJ, (rep1.P*pow(-1,(two_J1/2))));
                                rep_curr.n = find_n_subduced_embeddings(LG_curr, rep_curr.irrep, rep_curr.twoJ, (rep_curr.P*pow(-1,(two_J2/2))));
                                rep3.n = find_n_subduced_embeddings(LG3, rep3.irrep, rep3.twoJ, (rep3.P*pow(-1,(two_J3/2))));


                                for(int row1 = 1; row1 <= irrepRows(*p); row1++){     // Looping over all irrep rows at source, sink & insertion
                                  for(int row_curr = 1; row_curr <= irrepRows(*r); row_curr++){
                                    for(int row3 = 1; row3 <= irrepRows(*q); row3++){

                                        rep1.row = row1; rep3.row = row3; rep_curr.row = row_curr;


                                        sub_hel Sub1 = Subduce_all(mom1_sq, m1_sq, two_J1 , rep1, LG1, r1[0], r1[1], r1[2]);
                                        sub_hel Sub3 = Subduce_all(mom3_sq, m3_sq, two_J3 , rep3, LG3, r3[0], r3[1], r3[2]);
                                        sub_hel SubCurr = Subduce_all(mom_curr_sq, m_curr_sq, two_J2 , rep_curr, LG_curr, r_curr[0], r_curr[1], r_curr[2]);




                                        Coeff = KFacSSV(qp,qm,Sub1.sum,SubCurr.sum,Sub3.sum);
                                        if(std::real(Coeff) || std::imag(Coeff)){
                                          cout << mom1.transpose() << rep1.irrep << "["<< rep1.row <<"]" << "\n";
                                          cout << mom_curr.transpose()  << rep_curr.irrep << "["<< rep_curr.row <<"]"<< "\n";
                                          cout << mom3.transpose()  << rep3.irrep << "["<< rep3.row <<"]"<< "\n";
                                        //cout << "The abs_factor is:" << pow(std::real(Coeff),2)+pow(std::imag(Coeff),2) << "\n";
                                          cout << "The factor is:" << Coeff << "\n";
                                        }
                      }}}
                   }}}
                  }
                 }
             }}}
           }
         }}}
  }
};

