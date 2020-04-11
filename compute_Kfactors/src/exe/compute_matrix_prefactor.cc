/* main code to return the irreps to compute */

#include "lib/kfactors.h"

// semble
#include"/Users/archanar/LQCDSoftware/three_pt_analysis/semble_install/include/semble/semble_file_management.h"
#include "/Users/archanar/LQCDSoftware/three_pt_analysis/semble_install/include/semble/semble_meta.h"



int main(int argc, char** argv){

  //=================
  KFactorEnv::registerAll();
  KFactorEnv::registerAll();
  //=================

  if( argc != 3 ){
    cerr << "get_prefactors <hadron_xml> <output_xml>\n ";
    exit(1); }
    
  std::string in;      {istringstream a(argv[1]); a >> in;};
  std::string out;     {istringstream a(argv[2]); a >> out;};

    

    
  int npt, ui, L;
  double Xi, XiE;
  std::vector<hadron> had;
  hadron had_tmp;
  string matrix_type; 
  int num_matrix_elem;
  Array1dO<string> matrix_name;
  string compute_phase;
  string print_zero;


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
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/name", had_tmp.name);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/twoJ", had_tmp.twoJ);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/P", had_tmp.P);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/ell", had_tmp.ell);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/levels", had_tmp.levels);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/Elab", had_tmp.elab);
    read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/max_mom", had_tmp.max_mom);

    had.push_back(had_tmp);
    }
    read(xml_in, "/hadron/MatrixType", matrix_type);
    read(xml_in, "/hadron/numMatrixElem", num_matrix_elem);
    matrix_name.resize(num_matrix_elem);
    read(xml_in, "/hadron/MatrixName", matrix_name);
    read(xml_in, "/hadron/PrintZero", print_zero);

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
  int two_J2 = had[1].twoJ;
  int P2 = had[1].P*pow(-1,had[1].ell);
  int two_J3 = had[2].twoJ;
  int P3 = had[2].P*pow(-1,had[2].ell);
  int max_mom1 = had[0].max_mom;
  int max_mom2 = had[1].max_mom;
  int max_mom3 = had[2].max_mom;

  Array1dO<string> levs1 = had[0].levels;
  Array1dO<string> levs3 = had[2].levels;


  vector<complex<double>> Coeff;
  KFactor* kfac;
  Ph::phChars phase;

  vector<complex<double>> r_Coeff;
  KFactor* r_kfac;
  Ph::phChars r_phase;


  //output xml
  XMLFileWriter xml_out(out);
  push(xml_out, "kfac");

  write(xml_out, "L", L);
  write(xml_out, "Xi", Xi);
  write(xml_out, "XiError", XiE);
  write(xml_out, "ElabFilesIn", had[2].elab);
  write(xml_out, "ElabFilesOut", had[0].elab);
  write(xml_out, "kfacFile", matrix_name);
      
  
  //====================================================================================================
  //LOOP OVER MOMS AT SOURCE AND SINK,THE IRREPS AND IRREP ROWS AT THE 3 PTS TO GET NON ZERO PREFACTORS
  //====================================================================================================


  int count =0;

  Vector3d mom1(3,1);    Vector3d mom3(3,1);    Vector3d mom_curr(3,1);    Vector3d r_mom_curr(3,1);    Vector3d l_mom_curr(3,1);

  double mom_coeff = (2*PI)/(Xi*L);
  double mom_coeff_sq = pow(mom_coeff,2);

  // Looping over all mom_1
  vector<Vector3d> mom1_iter = iter::itermom(max_mom1);
  vector<Vector3d> mom3_iter = iter::itermom(max_mom3);  

  for(auto it1 = mom1_iter.begin(); it1 != mom1_iter.end(); it1++ ){
    mom1 = *it1;
    double mom1_sq = mom1.squaredNorm();

    for(auto it2 = mom3_iter.begin(); it2 != mom3_iter.end(); it2++ ){  
      mom3 = *it2;   
        
      double mom3_sq = mom3.squaredNorm();
      mom_curr = mom1 - mom3;

      double mom_curr_sq = mom_curr.squaredNorm();
                  
      // Cut-off on mom_curr
      if(max_mom2 >= mom_curr_sq){

      XMLArray::Array<int> mom_tmp(3); mom_tmp[0] = mom1(0); mom_tmp[1] = mom1(1); mom_tmp[2] = mom1(2);
      XMLArray::Array<int> canon_mom_1 = Hadron::canonicalOrder(mom_tmp);

      mom_tmp[0] = mom3(0); mom_tmp[1] = mom3(1); mom_tmp[2] = mom3(2);
      XMLArray::Array<int> canon_mom_2 = Hadron::canonicalOrder(mom_tmp);

      mom_tmp[0] = mom_curr(0); mom_tmp[1] = mom_curr(1); mom_tmp[2] = mom_curr(2);
      XMLArray::Array<int> mom_curr_can = Hadron::canonicalOrder(mom_tmp);

      if((mom_tmp == mom_curr_can))  //&& (mom_curr_sq != 0)
      {

      r_phase = Ph::phaseFactor(two_J1, two_J3, two_J2, mom1, mom3, compute_phase=="true"?true:false);
      phase   = Ph::phaseFactor(two_J1, two_J3, two_J2, mom1, mom3, false);


      r_mom_curr = r_phase.mom1 - r_phase.mom2;

      double mom_in_sq   = mom_coeff_sq*mom3_sq;
      double mom_out_sq  = mom_coeff_sq*mom1_sq;
      double mom_c_sq    = mom_coeff_sq*mom_curr_sq;



      string LG1     = generateLittleGroup(mom1);
      string LG3     = generateLittleGroup(mom3);
      string LG_curr = generateLittleGroup(mom_curr);


      //the ref angles for each mom from adat using the convention in 10.1103/PhysRevD.85.014507
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
      for(auto p = irrep1.begin(); p != irrep1.end(); p++){for(auto q = irrep3.begin(); q != irrep3.end(); q++){for(auto r = irrep_curr.begin(); r != irrep_curr.end(); r++){

        rep1.irrep = *p; rep3.irrep = *q; rep_curr.irrep = *r;

        //find embeddings, always 1
        rep1.n = find_n_subduced_embeddings(LG1, rep1.irrep, rep1.twoJ, (rep1.P*pow(-1,(two_J1/2))));
        rep_curr.n = find_n_subduced_embeddings(LG_curr, rep_curr.irrep, rep_curr.twoJ, (rep_curr.P*pow(-1,(two_J2/2))));
        rep3.n = find_n_subduced_embeddings(LG3, rep3.irrep, rep3.twoJ, (rep3.P*pow(-1,(two_J3/2))));


        // Looping over all irrep rows at source, sink & insertion
        for(int row1 = 1; row1 <= irrepRows(*p); row1++){for(int row_curr = 1; row_curr <= irrepRows(*r); row_curr++){for(int row3 = 1; row3 <= irrepRows(*q); row3++){

          for(int lev1 = 1; lev1 <= levs1.size(); lev1++){
            for(int lev3 = 1; lev3 <= levs3.size(); lev3++){

              string level_1 = levs1[lev1];
              string level_3 = levs3[lev3];
              

              rep1.row = row1; rep3.row = row3; rep_curr.row = row_curr;


              // Find the ELab jackknife file for the source and sink

              string E1str, E3str, p1str, p3str;
              string E1_name, E3_name;

              p1str = std::string("p") + to_string(canon_mom_1[0]) + to_string(canon_mom_1[1]) + to_string(canon_mom_1[2]);
              p3str = std::string("p") + to_string(canon_mom_2[0]) + to_string(canon_mom_2[1]) + to_string(canon_mom_2[2]);

              

              if(LG1 == "Oh"){E1str =  rep1.irrep;}
              else{E1str = LG1 + rep1.irrep;}


              if(LG3 == "Oh"){E3str = rep3.irrep;}
              else{E3str = LG3 + rep3.irrep;}    
                            

              cout << "Ei: " << E3str << endl; cout << "Ef: " << E1str << endl;

              bool found = false;

              for(int elab_count = 1; elab_count <= had[0].elab.size(); elab_count++){
                if((had[0].elab[elab_count].find(E1str) != std::string::npos) && (had[0].elab[elab_count].find(p1str) != std::string::npos) 
                          && (had[0].elab[elab_count].find(level_1) != std::string::npos) ){
                  found = true;
                  E1_name = had[0].elab[elab_count];
                  cout << "Reading Ef: " << had[0].elab[elab_count] << " " << level_1 << endl;
                  break;
                }
                else{continue;}
              }

              if(!found){break;}
              found = false;

              for(int elab_count = 1; elab_count <= had[2].elab.size(); elab_count++){
                if( (had[2].elab[elab_count].find(E3str) != std::string::npos ) && (had[2].elab[elab_count].find(p3str) != std::string::npos )
                            && (had[2].elab[elab_count].find(level_3) != std::string::npos) ){
                  found = true;    
                  E3_name = had[2].elab[elab_count];
                  cout << "Reading Ei: " << had[2].elab[elab_count] << " " << level_3 << endl;
                  break;
                    }

                else{continue;}
              }
              if(!found){break;}

              //Read the Elab jack files

              EnsemReal Elab1; read(E1_name, Elab1); Elab1 = rescaleEnsemDown(Elab1);
              EnsemReal Elab3; read(E3_name, Elab3); Elab3 = rescaleEnsemDown(Elab3);
              vector<EnsemComplex> prefactor(num_matrix_elem); 
              vector<EnsemComplex> r_prefactor(num_matrix_elem);
      
              if(Elab1.size() != Elab3.size()){cerr << "Ensemble sizes: " << E1_name << " and " << E3_name << "don't match"; exit(1); }
              
              int num = 0;
              while(num < num_matrix_elem){ prefactor[num].resize(Elab1.size());  r_prefactor[num].resize(Elab1.size()); num++;}

              double check_zero = 0.0;
              Ph::tripKey two_abs_lam;
              double  q2 = 0.0;

              for(int bin = 0; bin < Elab1.size(); bin++){

                double E1, E3;

                E1 = SEMBLE::toScalar(Elab1.elem(bin));
                E3 = SEMBLE::toScalar(Elab3.elem(bin));

                double m1_sq, m3_sq;
                m1_sq = pow(E1,2) - (mom_coeff_sq*mom1_sq);
                m3_sq = pow(E3,2) - (mom_coeff_sq*mom3_sq);

                VectorXd  qp(4,1);  VectorXd  qm(4,1);    // q+ = (p1-p2) q- = (p1+p2)
                VectorXd  r_qp(4,1);  VectorXd  r_qm(4,1);    // q+ = R(p1-p2) q- = R(p1+p2)

                
                r_qp << (E3+E1),-mom_coeff*(r_phase.mom1(0) + r_phase.mom2(0)),
                                                            -mom_coeff*(r_phase.mom1(1) + r_phase.mom2(1)),-mom_coeff*(r_phase.mom1(2) + r_phase.mom2(2));

                r_qm  << (E1 - E3),-mom_coeff*(-r_phase.mom2(0)+r_phase.mom1(0)),
                                                            -mom_coeff*(-r_phase.mom2(1)+r_phase.mom1(1)),-mom_coeff*(-r_phase.mom2(2)+r_phase.mom1(2));


                qp << (E3 + E1),-mom_coeff*(mom1(0) + mom3(0)),-mom_coeff*(mom1(1) + mom3(1)),-mom_coeff*(mom1(2) + mom3(2));
                qm  << (E1 - E3),-mom_coeff*(-mom3(0)+mom1(0)),-mom_coeff*(-mom3(1)+mom1(1)),-mom_coeff*(-mom3(2)+mom1(2));

                q2 += (- pow(qm(0),2) + pow(qm(1),2) + pow(qm(2),2) + pow(qm(3),2))/Elab1.size();

                double m_curr_sq =  pow(E1 - E3 ,2) - (mom_coeff_sq*mom_curr.squaredNorm());                           


                map< int, Eigen::MatrixXcd > Sub1    = Subduce_with_pol(mom_out_sq, m1_sq, two_J1 , rep1, LG1, r1[0], r1[1], r1[2], false);
                map< int, Eigen::MatrixXcd > Sub3    = Subduce_with_pol(mom_in_sq, m3_sq, two_J3 , rep3, LG3, r3[0], r3[1], r3[2], false);
                map< int, Eigen::MatrixXcd >SubCurr  = Subduce_with_pol(mom_c_sq, m_curr_sq, two_J2 , rep_curr, LG_curr, r_curr[0], r_curr[1], r_curr[2], true);


                map< int, Eigen::MatrixXcd > r_Sub1    = Subduce_with_pol(mom_out_sq, m1_sq, two_J1 , rep1, LG1, r_r1[0], r_r1[1], r_r1[2], false);
                map< int, Eigen::MatrixXcd > r_Sub3    = Subduce_with_pol(mom_in_sq, m3_sq, two_J3 , rep3, LG3, r_r3[0], r_r3[1], r_r3[2], false);
                map< int, Eigen::MatrixXcd > r_SubCurr = Subduce_with_pol(mom_c_sq, m_curr_sq, two_J2 , rep_curr, LG_curr, r_curr[0], r_curr[1], r_curr[2], true);

                KFacParams* kfac_params   = new KFacParams(Sub1,SubCurr,Sub3,phase,qp,qm);
                KFacParams* r_kfac_params = new KFacParams(r_Sub1,SubCurr,r_Sub3,r_phase,r_qp,r_qm);


                // for scalar vector with vector insertion

                XMLReader xml;
                kfac   = TheKFactorFactory::Instance().createObject(matrix_type, xml, "/Stuff");
                r_kfac = TheKFactorFactory::Instance().createObject(matrix_type, xml, "/Stuff");

                
                Coeff   = (*kfac)(*kfac_params);
                r_Coeff = (*r_kfac)(*r_kfac_params); 


                for(int num = 0; num < num_matrix_elem; num++){

                  Complex cc = SEMBLE::toScalar(Coeff[num]); Complex r_cc = SEMBLE::toScalar(r_Coeff[num]);
                  prefactor[num].elem(bin) = cc.elem();  r_prefactor[num].elem(bin) = r_cc.elem(); 

                }

                two_abs_lam = (*kfac_params).two_abs_lam();

                delete kfac_params; delete r_kfac_params;

              }   


              //============================
              //======= DIR NAMING ======== 

              /* Creates a subdirectory based on the irrep names. The kfac jackknife file for the specific irrep is stored 
              in the subdirectory and when the fitting code is run, it uses the same directory naming to access the kfac files */

                string name = naming::name(npt,two_abs_lam,r_phase.mom1,r_mom_curr,r_phase.mom2,rep1,rep_curr,rep3,LG1,LG_curr,LG3,level_1,level_3);
                string name_irrep = naming::name(npt,two_abs_lam,mom1,mom_curr,mom3,rep1,rep_curr,rep3,LG1,LG_curr,LG3,level_1,level_3);

                std::stringstream ss;
                ss << "Q2_";
                ss << std::fixed << std::setprecision(6) << q2; 
                ss <<  "/" << name << "/" << name_irrep;
                std::string path = SEMBLE::SEMBLEIO::getPath() += ss.str();
                SEMBLE::SEMBLEIO::makeDirectoryPath(path);
                path += std::string("/");

                cout << path << endl;


              for(int num = 0; num < num_matrix_elem; num++){

                prefactor[num]   = rescaleEnsemUp(prefactor[num]);
                r_prefactor[num] = rescaleEnsemUp(r_prefactor[num]);                      
                  

                cout << mom1.transpose() << rep1.irrep << "["<< rep1.row <<"]" << "\n";
                cout << mom_curr.transpose() << rep_curr.irrep << "["<< rep_curr.row <<"]"<< "\n";
                cout << mom3.transpose() << rep3.irrep << "["<< rep3.row <<"]"<< "\n";  

                cout << "The rotated " << num << " factor is:" << ENSEM::mean(r_prefactor[num]) << endl;
                cout << "The factor " << num <<  " is:" << ENSEM::mean(prefactor[num]) << endl;

                if( std::abs(abs(SEMBLE::toScalar(ENSEM::mean(prefactor[num]))) - abs(SEMBLE::toScalar(ENSEM::mean(r_prefactor[num])))) > std::numeric_limits<float>::epsilon() ){cout << "Discrepency!" << endl; }
                if( std::abs(abs(SEMBLE::toScalar(ENSEM::mean(prefactor[num]))) - abs(SEMBLE::toScalar(ENSEM::mean(r_prefactor[num])))) < std::numeric_limits<float>::epsilon() ){cout << "No Phase at all" << endl; }

                {
                  ostringstream outfile; outfile << path << matrix_name[num+1]; 
                  write(outfile.str(), prefactor[num]);
                }

                {
                  ostringstream outfile; outfile << path << "Rot" << matrix_name[num+1]; 
                  write(outfile.str(), r_prefactor[num]);
                }

                check_zero += abs(SEMBLE::toScalar(ENSEM::mean(prefactor[num])));
              }


              //probably a stupid way to write the file but works

              //=================
              //WRITE THE OUTPUT
              //=================

              //write in xml if non-zero   

              if(check_zero){

                ostringstream outfile; outfile << path << "../non_zero_correlators.txt";
                std::ofstream file;
                file.open(outfile.str(), std::ios_base::app);
                file << name_irrep << endl; 
              }
              else
              {
                ostringstream outfile; outfile << path << "../zeroed_correlators.txt";
                std::ofstream file;
                file.open(outfile.str(), std::ios_base::app);
                file << name_irrep << endl;                                     
              }

              if(check_zero || (print_zero=="true"?true:false )){

                push(xml_out, "elem");
                push(xml_out, "elem");
                if(mom1_sq){write(xml_out, "irrep", LG1+rep1.irrep);}
                else{write(xml_out, "irrep", rep1.irrep);}
                write(xml_out, "row", rep1.row);
                write_ei(xml_out, "mom", mom1);
                write(xml_out, "psq", mom1_sq);
                write(xml_out, "level", levs1[lev1]);
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
                write(xml_out, "level", levs3[lev3]);
                write(xml_out, "absLam", get<1>(two_abs_lam));
                pop(xml_out);
                pop(xml_out);
                Array1dO<Complex> tmp(num_matrix_elem);

                for(int num = 0; num < num_matrix_elem; num++){tmp[num+1] = cmplx( Real(std::real(SEMBLE::toScalar(ENSEM::mean(prefactor[num])))) , Real(std::imag(SEMBLE::toScalar(ENSEM::mean(prefactor[num])))) );}

                write(xml_out, "FF", tmp );
                count++; 

              }

          }}
                              
          }}}
        }}}
      }}
      }
}

  
  
  write(xml_out, "pts", count);

  push(xml_out, "hadron");
  write(xml_out, "npt", npt);
  write(xml_out, "L", L);
  write(xml_out, "Xi", Xi);
  write(xml_out, "XiError", XiE);

  push(xml_out, "had");

  for(int j =0;j<npt;j++){
    push(xml_out, "elem");
    write(xml_out, "name", had[j].name);
    write(xml_out, "twoJ", had[j].twoJ);
    write(xml_out, "P", had[j].P);
    write(xml_out, "ell", had[j].ell);
    write(xml_out, "Elab", had[j].elab);
    write(xml_out, "max_mom", had[j].max_mom);
    pop(xml_out);
  }

  pop(xml_out);
  write(xml_out, "MatrixType", matrix_type);
  write(xml_out, "MatrixName", matrix_name);
  write(xml_out, "numMatrixElem", num_matrix_elem);
  write(xml_out, "PrintZero", print_zero);  

  write(xml_out, "phase", compute_phase);

  pop(xml_out);


  pop(xml_out);


  xml_out.close();
  
  cout << count << "\n" ;

  delete kfac; delete r_kfac;
  
};

