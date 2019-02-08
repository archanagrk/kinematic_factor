/* main code to return the irreps to compute */

#include "lib/kfactors.h"



int main(int argc, char** argv){

  //=================
  KFactorEnv::registerAll();
  KFactorEnv::registerAll();
  //=================

  if( argc != 3 ){
    cerr << "get_non_zero_prefactors <hadron_xml> <output_xml>\n ";
    exit(1); }
    
  std::string in;      {istringstream a(argv[1]); a >> in;};
  std::string out;     {istringstream a(argv[2]); a >> out;};

    

    
  int npt;
  std::vector<hadron> had;
  int ui;
  hadron had_tmp;
  string matrix_type; 
  string compute_phase;


  //==============================
  //READ THE INPUT XML hadron_xml
  //==============================
  
  XMLReader xml_in(in);
  read(xml_in, "/hadron/npt", npt);
  for(int j =1;j<=npt;j++){
  read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/twoJ", had_tmp.twoJ);
  read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/P", had_tmp.P);
  read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/msq", had_tmp.msq);
  read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/max_mom", had_tmp.max_mom);
  had.push_back(had_tmp);
  }
  read(xml_in, "/hadron/matrix_type", matrix_type);
  read(xml_in, "/hadron/phase", compute_phase);
  
  //=================
  //LOAD DATA
  //=================

  int two_J1 = had[0].twoJ;
  int P1  = had[0].P;
  double m1_sq = had[0].msq;
  int two_J2 = had[1].twoJ;
  int P2 = had[1].P;
  int two_J3 = had[2].twoJ;
  int P3 = had[2].P;
  double m3_sq = had[2].msq;
  int max_mom1 = had[0].max_mom;
  int max_mom2 = had[1].max_mom;
  int max_mom3 = had[2].max_mom;

  complex<double> Coeff;
  KFactor* kfac;
  Ph::phChars phase;


  //output xml
  XMLFileWriter xml_out(out);
  push(xml_out, "kfac");
      
  
  //====================================================================================================
  //LOOP OVER MOMS AT SOURCE AND SINK,THE IRREPS AND IRREP ROWS AT THE 3 PTS TO GET NON ZERO PREFACTORS
  //====================================================================================================


  int count =0;

  Vector3d mom1(3,1);    Vector3d mom3(3,1);    Vector3d mom_curr(3,1);    Vector3d n_mom_curr(3,1);

  // Looping over all mom_1
  for(int i=-int(sqrt(max_mom1)); i <= int(sqrt(max_mom1)); i++){
    for(int j = -int(sqrt(max_mom1)); j <= int(sqrt(max_mom1)); j++){
      for(int k = -int(sqrt(max_mom1)); k <= int(sqrt(max_mom1)); k++){


        double mom1_sq = (pow(i,2)+pow(j,2)+pow(k,2));
        
        // Cut-off on mom1
        if(max_mom1 >= mom1_sq){
          
          // Looping over all mom_3
          for(int l = -int(sqrt(max_mom3)); l <= int(sqrt(max_mom3)); l++){
            for(int m = -int(sqrt(max_mom3)); m <= int(sqrt(max_mom3)); m++){
              for(int n = -int(sqrt(max_mom3)); n<= int(sqrt(max_mom3)); n++){

                double mom3_sq = (pow(l,2)+pow(m,2)+pow(n,2));
                
                // Cut-off on mom2
                if(max_mom3 >= mom3_sq){

                  mom1 << i,j,k;  mom3 << l,m,n;  mom_curr << (l-i),(m-j),(n-k);

                  double mom_curr_sq = (pow((i-l),2)+pow((j-m),2)+pow((k-n),2));
                  
                  // Cut-off on mom_curr
                  if(max_mom2 >= mom_curr_sq){


                    phase = Ph::phaseFactor(two_J1, two_J3, two_J2, mom1, mom3, compute_phase=="true"?true:false);

                    double n_mom1_sq = phase.mom1.squaredNorm();
                    double n_mom3_sq = phase.mom2.squaredNorm();

                    n_mom_curr = phase.mom2 - phase.mom1;



                    VectorXd  qp(4,1);  VectorXd  qm(4,1);    // q+ = (p1-p2) q- = (p1+p2)

                    
                    qp << (sqrt(m1_sq+n_mom1_sq)+sqrt(m3_sq+n_mom3_sq)),-(phase.mom1(0) + phase.mom2(0)),-(phase.mom1(1) + phase.mom2(1)),-(phase.mom1(2) + phase.mom2(2));
                    qm  << (sqrt(m3_sq+n_mom3_sq)-sqrt(m1_sq+n_mom1_sq)),-(phase.mom2(0)-phase.mom1(0)),-(phase.mom2(1)-phase.mom1(1)),-(phase.mom2(2)-phase.mom1(2));


                    double m_curr_sq =  pow(sqrt(n_mom3_sq + m3_sq)-sqrt(n_mom1_sq + m1_sq) ,2) - n_mom_curr.squaredNorm();


                    string LG1 = generateLittleGroup(phase.mom1);
                    string LG3 = generateLittleGroup(phase.mom2);
                    string LG_curr = generateLittleGroup(n_mom_curr);


                    //the ref angles for each mom from adat
		                std::vector<double> r1 = refAngles(phase.mom1);
		                std::vector<double> r_curr = refAngles(n_mom_curr);
		                std::vector<double> r3 = refAngles(phase.mom2);                       
                       

                    std::vector<std::string> irrep1 = getIrrep(two_J1,P1,LG1);
                    std::vector<std::string> irrep_curr = getIrrep(two_J2,P2,LG_curr);
                    std::vector<std::string> irrep3 = getIrrep(two_J3,P3,LG3);

                    irrep_label rep1; irrep_label rep_curr; irrep_label rep3;


                    rep1.twoJ = two_J1; rep3.twoJ = two_J3; rep_curr.twoJ = two_J2;
                    rep1.P = P1; rep3.P = P3; rep_curr.P = P2;

                    
                    // Looping over all irreps at the source, sink & isertion
                    for(auto p = irrep1.begin(); p != irrep1.end(); p++){
                      for(auto q = irrep3.begin(); q != irrep3.end(); q++){
                        for(auto r = irrep_curr.begin(); r != irrep_curr.end(); r++){

                          rep1.irrep = *p; rep3.irrep = *q; rep_curr.irrep = *r;

                          //find embeddings, always 1
                          rep1.n = find_n_subduced_embeddings(LG1, rep1.irrep, rep1.twoJ, (rep1.P*pow(-1,(two_J1/2))));
                          rep_curr.n = find_n_subduced_embeddings(LG_curr, rep_curr.irrep, rep_curr.twoJ, (rep_curr.P*pow(-1,(two_J2/2))));
                          rep3.n = find_n_subduced_embeddings(LG3, rep3.irrep, rep3.twoJ, (rep3.P*pow(-1,(two_J3/2))));


                          // Looping over all irrep rows at source, sink & insertion
                          for(int row1 = 1; row1 <= irrepRows(*p); row1++){
                            for(int row_curr = 1; row_curr <= irrepRows(*r); row_curr++){
                              for(int row3 = 1; row3 <= irrepRows(*q); row3++){

                                rep1.row = row1; rep3.row = row3; rep_curr.row = row_curr;


                                map< int, Eigen::MatrixXcd > Sub1    = Subduce_with_pol(mom1_sq, m1_sq, two_J1 , rep1, LG1, r1[0], r1[1], r1[2]);
                                map< int, Eigen::MatrixXcd > Sub3    = Subduce_with_pol(mom3_sq, m3_sq, two_J3 , rep3, LG3, r3[0], r3[1], r3[2]);
                                map< int, Eigen::MatrixXcd >SubCurr  = Subduce_with_pol(mom_curr_sq, m_curr_sq, two_J2 , rep_curr, LG_curr, r_curr[0], r_curr[1], r_curr[2]);

                                //int helicity1 = Sub1.two_hel/2;
                                //int helicity_curr = SubCurr.two_hel/2;
                                //int helicity3 = Sub3.two_hel/2;

                                KFacParams* kfac_params = new KFacParams(Sub1,SubCurr,Sub3,phase,qp,qm);
                                        
                                // for scalar vector with vector insertion

                                XMLReader xml;
                                kfac = TheKFactorFactory::Instance().createObject(matrix_type, xml, "/Stuff");

                                Coeff = (*kfac)(*kfac_params);

                                Ph::tripKey two_abs_lam = (*kfac_params).two_abs_lam();
                                  

                                if(std::real(Coeff) || std::imag(Coeff)){
                                  
                                  cout << mom1.transpose() << rep1.irrep << "["<< rep1.row <<"]" << "\n";
                                  cout << mom_curr.transpose()  << rep_curr.irrep << "["<< rep_curr.row <<"]"<< "\n";
                                  cout << mom3.transpose()  << rep3.irrep << "["<< rep3.row <<"]"<< "\n";  
                                  cout << "The factor is:" << Coeff << endl;                              



                                  //=================
                                  //WRITE THE OUTPUT
                                  //=================

                                  //write in xml                                
                                  push(xml_out, "elem");
                                  push(xml_out, "elem");
                                  if(mom1_sq){write(xml_out, "irrep", LG1+rep1.irrep);}
                                  else{write(xml_out, "irrep", rep1.irrep);}
                                  write(xml_out, "row", rep1.row);
                                  write_ei(xml_out, "mom", mom1);
                                  write(xml_out, "psq", mom1_sq);
                                  write(xml_out, "abs_lam", get<0>(two_abs_lam));
                                  pop(xml_out);
                                  
                                  push(xml_out, "elem");
                                  if(mom_curr_sq){write(xml_out, "irrep", LG_curr+rep_curr.irrep);}
                                  else{write(xml_out, "irrep", rep_curr.irrep);}
                                  write(xml_out, "row", rep_curr.row);
                                  write_ei(xml_out, "mom", mom_curr);
                                  write(xml_out, "psq", mom_curr_sq);
                                  write(xml_out, "abs_lam", get<2>(two_abs_lam));
                                  pop(xml_out);
                                  
                                  push(xml_out, "elem");
                                  if(mom3_sq){write(xml_out, "irrep", LG3+rep3.irrep);}
                                  else{write(xml_out, "irrep", rep3.irrep);}
                                  write(xml_out, "row", rep3.row);
                                  write_ei(xml_out, "mom", mom3);
                                  write(xml_out, "psq", mom3_sq);
                                  write(xml_out, "abs_lam", get<1>(two_abs_lam));
                                  pop(xml_out);
                                  pop(xml_out);
                                  count++;


                                  delete kfac_params;
                                }
                                        
                      }}}
                   }}}
                  }
                 }
             }}}
           }
         }}}
  
  
  write(xml_out, "pts", count);
  pop(xml_out);
  xml_out.close();
  
  cout << count << "\n" ;

  delete kfac;
  
};

