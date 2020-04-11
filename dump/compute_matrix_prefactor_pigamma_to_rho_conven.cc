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

    

    
  int npt, ui, L;
  double Xi, XiE;
  std::vector<hadron> had;
  hadron had_tmp;
  string matrix_type; 
  string compute_phase;


  //==============================
  //READ THE INPUT XML hadron_xml
  //==============================
  
  XMLReader xml_in(in);
   try{
    read(xml_in, "/hadron/npt", npt);
    read(xml_in, "/hadron/L", L);
    read(xml_in, "/hadron/Xi", Xi);
    read(xml_in, "/hadron/XiError", XiE);

    for(int j =1;j<=npt;j++){
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/twoJ", had_tmp.twoJ);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/P", had_tmp.P);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/ell", had_tmp.ell);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/msq", had_tmp.msq);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/max_mom", had_tmp.max_mom);

    had.push_back(had_tmp);
    }
    read(xml_in, "/hadron/matrix_type", matrix_type);
    read(xml_in, "/hadron/phase", compute_phase);
   }
    catch( const string& error ){
    cerr << "Error reading hadron file : " << error << endl;
    }
  
  //=================
  //LOAD DATA
  //=================

  int two_J1 = had[0].twoJ;
  int P1  = had[0].P*pow(-1,had[0].ell);
  double m1_sq = had[0].msq;
  int two_J2 = had[1].twoJ;
  int P2 = had[1].P*pow(-1,had[1].ell);
  int two_J3 = had[2].twoJ;
  int P3 = had[2].P*pow(-1,had[2].ell);
  double m3_sq = had[2].msq;
  int max_mom1 = had[0].max_mom;
  int max_mom2 = had[1].max_mom;
  int max_mom3 = had[2].max_mom;

  complex<double> Coeff;
  KFactor* kfac;
  Ph::phChars phase;

  complex<double> r_Coeff;
  KFactor* r_kfac;
  Ph::phChars r_phase;

  complex<double> l_Coeff;
  KFactor* l_kfac;
  Ph::phChars l_phase;



  //output xml
  XMLFileWriter xml_out(out);
  push(xml_out, "kfac");
  write(xml_out, "L", L);
  write(xml_out, "Xi", Xi);
  write(xml_out, "XiE", XiE);
  write(xml_out, "m1Sq", m1_sq);
  write(xml_out, "m3Sq", m3_sq);
      
  
  //====================================================================================================
  //LOOP OVER MOMS AT SOURCE AND SINK,THE IRREPS AND IRREP ROWS AT THE 3 PTS TO GET NON ZERO PREFACTORS
  //====================================================================================================


  int count =0;

  Vector3d mom1(3,1);    Vector3d mom3(3,1);    Vector3d mom_curr(3,1);    Vector3d r_mom_curr(3,1);    Vector3d l_mom_curr(3,1);

  double mom_coeff = (2*PI)/(Xi*L);
  double mom_coeff_sq = pow(mom_coeff,2);

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

                    XMLArray::Array<int> mom_tmp(3); mom_tmp[0] = i; mom_tmp[1] = j; mom_tmp[2] = k;
                    XMLArray::Array<int> canon_mom_1 = Hadron::canonicalOrder(mom_tmp);

                    mom_tmp[0] = l; mom_tmp[1] = m; mom_tmp[2] = n;
                    XMLArray::Array<int> canon_mom_2 = Hadron::canonicalOrder(mom_tmp);

                    mom_tmp[0] = l-i; mom_tmp[1] = m-j; mom_tmp[2] = n-k;
                    XMLArray::Array<int> mom_curr_can = Hadron::canonicalOrder(mom_tmp);

                    if((mom_tmp == mom_curr_can) && (mom_curr_sq != 0))  
                    {

                    r_phase = Ph::phaseFactor(two_J1, two_J3, two_J2, mom1, mom3, compute_phase=="true"?true:false);
                    phase   = Ph::phaseFactor(two_J1, two_J3, two_J2, mom1, mom3, false);


                    r_mom_curr = r_phase.mom2 - r_phase.mom1;

                    double mom_in_sq   = mom_coeff_sq*mom1_sq;
                    double mom_out_sq  = mom_coeff_sq*mom3_sq;
                    double mom_c_sq    = mom_coeff_sq*mom_curr_sq;

                    VectorXd  qp(4,1);  VectorXd  qm(4,1);    // q+ = (p1-p2) q- = (p1+p2)
                    VectorXd  r_qp(4,1);  VectorXd  r_qm(4,1);    // q+ = R(p1-p2) q- = R(p1+p2)

                    
                    r_qp << (sqrt(m1_sq+ mom_in_sq)+sqrt(m3_sq+mom_out_sq)),-mom_coeff*(r_phase.mom1(0) + r_phase.mom2(0)),
                                                                -mom_coeff*(r_phase.mom1(1) + r_phase.mom2(1)),-mom_coeff*(r_phase.mom1(2) + r_phase.mom2(2));

                    r_qm  << (sqrt(m3_sq+mom_out_sq)-sqrt(m1_sq+mom_in_sq)),-mom_coeff*(r_phase.mom2(0)-r_phase.mom1(0)),
                                                                -mom_coeff*(r_phase.mom2(1)-r_phase.mom1(1)),-mom_coeff*(r_phase.mom2(2)-r_phase.mom1(2));


                    qp << (sqrt(m1_sq+mom_in_sq)+sqrt(m3_sq+mom_out_sq)),-mom_coeff*(mom1(0) + mom3(0)),-mom_coeff*(mom1(1) + mom3(1)),-mom_coeff*(mom1(2) + mom3(2));
                    qm  << (sqrt(m3_sq+mom_out_sq)-sqrt(m1_sq+mom_in_sq)),-mom_coeff*(mom3(0)-mom1(0)),-mom_coeff*(mom3(1)-mom1(1)),-mom_coeff*(mom3(2)-mom1(2));


                    double m_curr_sq =  pow(sqrt( mom_out_sq + m3_sq) - sqrt( mom_in_sq + m1_sq) ,2) - (mom_coeff_sq*mom_curr.squaredNorm());


                    string LG1     = generateLittleGroup(mom1);
                    string LG3     = generateLittleGroup(mom3);
                    string LG_curr = generateLittleGroup(mom_curr);


                    //the ref angles for each mom from adat
		                std::vector<double> r_r1     = refAngles(r_phase.mom1);
		                std::vector<double> r_r_curr = refAngles(r_mom_curr);
		                std::vector<double> r_r3     = refAngles(r_phase.mom2);     

		                std::vector<double> r1     = refAngles(mom1);
		                std::vector<double> r_curr = refAngles(mom_curr);
		                std::vector<double> r3     = refAngles(mom3);                     
                       

                    std::vector<std::string> irrep1     = getIrrep(two_J1,P1,LG1);
                    std::vector<std::string> irrep_curr = getIrrep(two_J2,P2,LG_curr);
                    std::vector<std::string> irrep3     = getIrrep(two_J3,P3,LG3);

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


                                map< int, Eigen::MatrixXcd > Sub1    = Subduce_with_pol(mom_in_sq, m1_sq, two_J1 , rep1, LG1, r1[0], r1[1], r1[2], false);
                                map< int, Eigen::MatrixXcd > Sub3    = Subduce_with_pol(mom_out_sq, m3_sq, two_J3 , rep3, LG3, r3[0], r3[1], r3[2], false);
                                map< int, Eigen::MatrixXcd >SubCurr  = Subduce_with_pol(mom_c_sq, m_curr_sq, two_J2 , rep_curr, LG_curr, r_curr[0], r_curr[1], r_curr[2], true);


                                map< int, Eigen::MatrixXcd > r_Sub1    = Subduce_with_pol(mom_in_sq, m1_sq, two_J1 , rep1, LG1, r_r1[0], r_r1[1], r_r1[2], false);
                                map< int, Eigen::MatrixXcd > r_Sub3    = Subduce_with_pol(mom_out_sq, m3_sq, two_J3 , rep3, LG3, r_r3[0], r_r3[1], r_r3[2], false);
                                map< int, Eigen::MatrixXcd > r_SubCurr = Subduce_with_pol(mom_c_sq, m_curr_sq, two_J2 , rep_curr, LG_curr, r_curr[0], r_curr[1], r_curr[2], true);


                                KFacParams* kfac_params   = new KFacParams(Sub1,SubCurr,Sub3,phase,qp,qm);
                                KFacParams* r_kfac_params = new KFacParams(r_Sub1,SubCurr,r_Sub3,r_phase,r_qp,r_qm);


                                // for scalar vector with vector insertion

                                XMLReader xml;
                                kfac   = TheKFactorFactory::Instance().createObject(matrix_type, xml, "/Stuff");
                                r_kfac = TheKFactorFactory::Instance().createObject(matrix_type, xml, "/Stuff");

                                Coeff   = (*kfac)(*kfac_params);
                                r_Coeff = (*r_kfac)(*r_kfac_params);                              

                                Ph::tripKey two_abs_lam = (*kfac_params).two_abs_lam();
                                  

                                //if( !std::real(Coeff) || !std::imag(Coeff) || !std::imag(r_Coeff) || !std::real(r_Coeff) ){
                                if( !std::real(Coeff) && !std::imag(Coeff) ){

                                  // cout << mom1.transpose() << rep1.irrep << "["<< rep1.row <<"]" << "\n";
                                  // cout << mom_curr.transpose() << rep_curr.irrep << "["<< rep_curr.row <<"]"<< "\n";
                                  // cout << mom3.transpose() << rep3.irrep << "["<< rep3.row <<"]"<< "\n";  

                                  // cout << "The rotated factor is:" << r_Coeff << endl;
                                  // cout << "The factor is:" << Coeff << endl;

                                  // if( std::abs(abs(Coeff) - abs(r_Coeff)) > std::numeric_limits<float>::epsilon() ){cout << "Discrepency!" << endl; }
                                  // if( std::abs(Coeff - r_Coeff) < std::numeric_limits<float>::epsilon() ){cout << "No Phase at all" << endl; }

                                  //=================
                                  //WRITE THE CORR FN
                                  //=================
                                  
                                  string corrfn = "";
                                  string tmpfn;

                                  for(int cz = 0; cz < 3; cz++){
                                    if(cz == 0){
                                      string canon = std::to_string(canon_mom_1[0]) + std::to_string(canon_mom_1[1]) + std::to_string(canon_mom_1[2]);
                                      string mst = std::to_string(-i) + std::to_string(-j) + std::to_string(-k);
                                      if(mom1_sq != 0){
                                        tmpfn = "t32,fI2Y0i2,r"+std::to_string(rep1.row)+","+ mst +",pion_proj0_p"+ canon + "_H"+ std::to_string(get<0>(two_abs_lam)/2) + LG1 + rep1.irrep + "__"+ canon;
                                      }
                                      else{
                                        tmpfn = "t32,fI2Y0i2,r"+std::to_string(rep1.row)+","+ mst +",pion_proj0_p"+ canon + rep1.irrep + "__"+ canon;
                                      }
                                    }
                                    if(cz == 1){
                                      string canon = std::to_string(mom_curr_can[0]) + std::to_string(mom_curr_can[1]) + std::to_string(mom_curr_can[2]);
                                      if(mom_curr_sq != 0){
                                        tmpfn = ".tm3,fI0Y0i0,r"+std::to_string(rep_curr.row)+","+ canon +",omegal_rhoxD0_J0__J1_H"+ std::to_string(get<1>(two_abs_lam)/2) + LG_curr + rep_curr.irrep + "__"+ canon;
                                      }
                                      else{
                                        tmpfn = ".tm3,fI0Y0i0,r"+std::to_string(rep_curr.row)+","+ canon +",omegal_rhoxD0_J0__J1_" + rep_curr.irrep + "__"+ canon;
                                      }

                                    }
                                    if(cz == 2){
                                      string canon = std::to_string(canon_mom_2[0]) + std::to_string(canon_mom_2[1]) + std::to_string(canon_mom_2[2]);
                                      string mst = std::to_string(-l) + std::to_string(-m) + std::to_string(-n);
                                      if(mom3_sq != 0){
                                        tmpfn = ".t0,fI2Y0i2,r"+std::to_string(rep3.row)+","+ mst +",rho_proj0_p"+ canon + "_H"+ std::to_string(get<2>(two_abs_lam)/2) + LG3 + rep3.irrep + "__"+ canon + ".dat";
                                      }
                                      else{
                                        tmpfn = ".t0,fI2Y0i2,r"+std::to_string(rep3.row)+","+ mst +",rho_proj0_p"+ canon + rep3.irrep + "__"+ canon + ".dat";
                                      }
                                    }

                                    corrfn += tmpfn;
                                  }

                                  cout << corrfn << endl;

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
                                  write(xml_out, "absLam", get<0>(two_abs_lam));
                                  pop(xml_out);
                                  
                                  push(xml_out, "elem");
                                  if(mom_curr_sq){write(xml_out, "irrep", LG_curr+rep_curr.irrep);}
                                  else{write(xml_out, "irrep", rep_curr.irrep);}
                                  write(xml_out, "row", rep_curr.row);
                                  write_ei(xml_out, "mom", mom_curr);
                                  write(xml_out, "psq", mom_curr_sq);
                                  write(xml_out, "absLam", get<2>(two_abs_lam));
                                  pop(xml_out);
                                  
                                  push(xml_out, "elem");
                                  if(mom3_sq){write(xml_out, "irrep", LG3+rep3.irrep);}
                                  else{write(xml_out, "irrep", rep3.irrep);}
                                  write(xml_out, "row", rep3.row);
                                  write_ei(xml_out, "mom", mom3);
                                  write(xml_out, "psq", mom3_sq);
                                  write(xml_out, "absLam", get<1>(two_abs_lam));
                                  pop(xml_out);
                                  write(xml_out, "cReal", std::real(Coeff));
                                  write(xml_out, "cImag", std::imag(Coeff));
                                  pop(xml_out);
                                  count++;


                                  delete kfac_params;
                                }
                                        
                      }}}
                   }}}
                  }}
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

