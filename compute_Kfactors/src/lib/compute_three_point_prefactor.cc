/* main code to return the irreps to compute */
#include "compute_three_point_prefactor.h"

void compute_three_point_prefactor(string& in, string& out){
    
  int npt, L;
  double Xi, XiE;
  std::vector<NPtCorr_t> had;
  string matrix_type; 
  int num_matrix_elem;
  Array1dO<string> matrix_name;
  string compute_phase;
  string print_zero, elab_dir;
  bool make_redstar_xml;
  string redstar_xml;
  vector<NPtIrrepLam_t> irreps;


  //==============================
  //READ THE INPUT XML hadron_xml
  //==============================

  XMLReader xml_in(in);

  read_xml_ini(xml_in, npt, L, Xi, XiE, had, num_matrix_elem, matrix_type, matrix_name, print_zero,
                        compute_phase, elab_dir, make_redstar_xml, redstar_xml );

  
  //=================
  //LOAD DATA
  //=================

  int P1  = had[0].P*pow(-1,had[0].ell);
  int P2 = had[1].P*pow(-1,had[1].ell);
  int P3 = had[2].P*pow(-1,had[2].ell);


  Array1dO<string> levs1 = had[0].levels;
  Array1dO<string> levs2 = had[1].levels;
  Array1dO<string> levs3 = had[2].levels;


  vector<complex<double>> Coeff;
  KFactor* kfac;
  Ph::phChars phase;

  vector<complex<double>> r_Coeff;
  KFactor* r_kfac;
  Ph::phChars r_phase;

      
  
  //====================================================================================================
  //LOOP OVER MOMS AT SOURCE AND SINK,THE IRREPS AND IRREP ROWS AT THE 3 PTS TO GET NON ZERO PREFACTORS
  //====================================================================================================


  int count =0;

  Vector3d mom1(3,1);    Vector3d mom3(3,1);    Vector3d mom_curr(3,1);    Vector3d r_mom_curr(3,1);    Vector3d l_mom_curr(3,1);

  double mom_coeff = (2*PI)/(Xi*L);
  double mom_coeff_sq = pow(mom_coeff,2);


  // Looping over all mom_1
  vector<Vector3d> mom1_iter = iter::itermom(had[0].max_mom, had[0].min_mom, had[0].omit_mom, had[0].canonical);
  vector<Vector3d> mom2_iter = iter::itermom(had[1].max_mom, had[1].min_mom, had[1].omit_mom, had[1].canonical);
  vector<Vector3d> mom3_iter = iter::itermom(had[2].max_mom, had[2].min_mom, had[2].omit_mom, had[2].canonical);
  
  
  vector<pair<EnsemReal,EnsemReal>> ecm1_qsq, ecm3_qsq;

  for(auto it1 = mom1_iter.begin(); it1 != mom1_iter.end(); it1++ ){
    mom1 = *it1;
    double mom1_sq = mom1.squaredNorm();
  
    for(auto it2 = mom3_iter.begin(); it2 != mom3_iter.end(); it2++ ){  
      mom3 = *it2;   
      
      double mom3_sq = mom3.squaredNorm();
      mom_curr = mom1 - mom3; 
      double mom_curr_sq = mom_curr.squaredNorm();
   
      bool ispresent = false;
      // Cut-off on mom_curr
      if(mom2_iter.size() > 0){
        for(auto iter = mom2_iter.begin(); iter != mom2_iter.end(); iter++){
          if(*iter == mom_curr) ispresent = true;
        }
      }


      if(ispresent)
      {
        XMLArray::Array<int> mom_tmp(3); mom_tmp[0] = mom1(0); mom_tmp[1] = mom1(1); mom_tmp[2] = mom1(2);
        XMLArray::Array<int> canon_mom_1 = Hadron::canonicalOrder(mom_tmp);

        mom_tmp[0] = mom3(0); mom_tmp[1] = mom3(1); mom_tmp[2] = mom3(2);
        XMLArray::Array<int> canon_mom_2 = Hadron::canonicalOrder(mom_tmp);


        r_phase = Ph::phaseFactor(had[0].twoJ, had[2].twoJ, had[1].twoJ, mom1, mom3, compute_phase=="true"?true:false);
        phase   = Ph::phaseFactor(had[0].twoJ, had[2].twoJ, had[1].twoJ, mom1, mom3, false);

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
            

        std::vector<std::string> irrep1     = getIrrep(had[0].twoJ,P1,LG1);
        std::vector<std::string> irrep_curr = getIrrep(had[1].twoJ,P2,LG_curr);
        std::vector<std::string> irrep3     = getIrrep(had[2].twoJ,P3,LG3);

        irrep_label rep1; irrep_label rep_curr; irrep_label rep3;


        rep1.twoJ = had[0].twoJ; rep3.twoJ = had[2].twoJ; rep_curr.twoJ = had[1].twoJ;
        rep1.P = P1; rep3.P = P3; rep_curr.P = P2;


        // Looping over all irreps at the source, sink & isertion
        for(auto p = irrep1.begin(); p != irrep1.end(); p++){for(auto q = irrep3.begin(); q != irrep3.end(); q++){for(auto r = irrep_curr.begin(); r != irrep_curr.end(); r++){

          rep1.irrep = *p; rep3.irrep = *q; rep_curr.irrep = *r;

          //find embeddings, always 1
          rep1.n = find_n_subduced_embeddings(LG1, rep1.irrep, rep1.twoJ, (rep1.P*pow(-1,(had[0].twoJ/2))));
          rep_curr.n = find_n_subduced_embeddings(LG_curr, rep_curr.irrep, rep_curr.twoJ, (rep_curr.P*pow(-1,(had[1].twoJ/2))));
          rep3.n = find_n_subduced_embeddings(LG3, rep3.irrep, rep3.twoJ, (rep3.P*pow(-1,(had[2].twoJ/2))));


          // Looping over all irrep rows at source, sink & insertion
          for(int row1 = 1; row1 <= irrepRows(*p); row1++){for(int row_curr = 1; row_curr <= irrepRows(*r); row_curr++){for(int row3 = 1; row3 <= irrepRows(*q); row3++){

            for(int lev1 = 1; lev1 <= levs1.size(); lev1++){
              for(int lev3 = 1; lev3 <= levs3.size(); lev3++){
                // for(int lev2 = 1; lev2 <= levs2.size(); lev2++){

                string level_1 = levs1[lev1];
                string level_3 = levs3[lev3];
                // string level_2 = levs2[lev2];


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

                EnsemReal Elab1; read(elab_dir+"/"+E1_name, Elab1); Elab1 = rescaleEnsemDown(Elab1);
                EnsemReal Elab3; read(elab_dir+"/"+E3_name, Elab3); Elab3 = rescaleEnsemDown(Elab3);
                vector<EnsemComplex> prefactor(num_matrix_elem); 
                vector<EnsemComplex> r_prefactor(num_matrix_elem);
                EnsemReal qsq; qsq.resize(Elab1.size());
                EnsemReal Ecm1; Ecm1.resize(Elab1.size());
                EnsemReal Ecm3; Ecm3.resize(Elab1.size());

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

                  // get the polarization dot subdction coefficient of each npt make a map of helicity, pol*Subd

                  map< int, Eigen::MatrixXcd > Sub1    = Subduce_with_pol(mom_out_sq, m1_sq, had[0].twoJ , rep1, LG1, r1[0], r1[1], r1[2], false);
                  map< int, Eigen::MatrixXcd > Sub3    = Subduce_with_pol(mom_in_sq, m3_sq, had[2].twoJ , rep3, LG3, r3[0], r3[1], r3[2], false);
                  map< int, Eigen::MatrixXcd >SubCurr  = Subduce_with_pol(mom_c_sq, m_curr_sq, had[1].twoJ , rep_curr, LG_curr, r_curr[0], r_curr[1], r_curr[2], true);


                  map< int, Eigen::MatrixXcd > r_Sub1    = Subduce_with_pol(mom_out_sq, m1_sq, had[0].twoJ , rep1, LG1, r_r1[0], r_r1[1], r_r1[2], false);
                  map< int, Eigen::MatrixXcd > r_Sub3    = Subduce_with_pol(mom_in_sq, m3_sq, had[2].twoJ , rep3, LG3, r_r3[0], r_r3[1], r_r3[2], false);
                  map< int, Eigen::MatrixXcd > r_SubCurr = Subduce_with_pol(mom_c_sq, m_curr_sq, had[1].twoJ , rep_curr, LG_curr, r_r_curr[0], r_r_curr[1], r_r_curr[2], true);

                  // get the kfactors with the form given in the input xml
                  KFacParams* kfac_params   = new KFacParams(Sub1,SubCurr,Sub3,phase,qp,qm);
                  KFacParams* r_kfac_params = new KFacParams(r_Sub1,SubCurr,r_Sub3,r_phase,r_qp,r_qm);

                  XMLReader xml;
                  kfac   = TheKFactorFactory::Instance().createObject(matrix_type, xml, "/Stuff");
                  r_kfac = TheKFactorFactory::Instance().createObject(matrix_type, xml, "/Stuff");

                  // get the coefficients
                  Coeff   = (*kfac)(*kfac_params);
                  r_Coeff = (*r_kfac)(*r_kfac_params); 

                  qsq.elem(bin) = - pow(qm(0),2) + pow(qm(1),2) + pow(qm(2),2) + pow(qm(3),2);
                  Ecm1.elem(bin) = sqrt(m1_sq);
                  Ecm3.elem(bin) = sqrt(m3_sq);


                  for(int num = 0; num < num_matrix_elem; num++){

                    Complex cc = SEMBLE::toScalar(Coeff[num]); Complex r_cc = SEMBLE::toScalar(r_Coeff[num]);
                    prefactor[num].elem(bin) = cc.elem();  r_prefactor[num].elem(bin) = r_cc.elem(); 

                  }

                  two_abs_lam = (*kfac_params).two_abs_lam();

                  delete kfac_params; delete r_kfac_params;

                }  

                //rescale the ensembles // 
                qsq   = rescaleEnsemUp(qsq);
                Ecm1 = rescaleEnsemUp(Ecm1);
                Ecm3 = rescaleEnsemUp(Ecm3);


                //============================
                //======= DIR NAMING ======== 

                /* Creates a subdirectory based on the irrep names. The kfac jackknife file for the specific irrep is stored 
                in the subdirectory and when the fitting code is run, it uses the same directory naming to access the kfac files */

                  string name = naming::name(npt,two_abs_lam,r_phase.mom1,r_mom_curr,r_phase.mom2,rep1,rep_curr,rep3,LG1,LG_curr,LG3,level_1,level_3);
                  string name_irrep = naming::name(npt,two_abs_lam,mom1,mom_curr,mom3,rep1,rep_curr,rep3,LG1,LG_curr,LG3,level_1,level_3);

                  std::stringstream ss;
                  ss << "Q2_";
                  ss << std::fixed << std::setprecision(6) << SEMBLE::toScalar(ENSEM::mean(qsq)); 
                  ss <<  "/" << name << "/" << name_irrep;
                  std::string path = SEMBLE::SEMBLEIO::getPath() += ss.str();
                  SEMBLE::SEMBLEIO::makeDirectoryPath(path);
                  path += std::string("/");

                  cout << path << endl;

                //======================
                //WRITE THE OUTPUT FILES
                //======================

                // number of matrix elems 
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
                //write in xml if non-zero   

                if(check_zero){

                  ostringstream outfile; outfile << path << "../non_zero_correlators.txt";
                  std::ofstream file;
                  file.open(outfile.str(), std::ios_base::app);
                  file << name_irrep << endl; 

                  ecm1_qsq.push_back(make_pair(Ecm1,qsq));
                  ecm3_qsq.push_back(make_pair(Ecm3,qsq));
                }

                else
                {
                  ostringstream outfile; outfile << path << "../zeroed_correlators.txt";
                  std::ofstream file;
                  file.open(outfile.str(), std::ios_base::app);
                  file << name_irrep << endl;                                     
                }

                if(check_zero || (print_zero=="true"?true:false )){

                  NPtIrrepLam_t irrep_tmp;
                  IrrepLam_t irrep_lam_tmp;

                  if(mom1_sq) irrep_lam_tmp.irrep = LG1+rep1.irrep;
                  else irrep_lam_tmp.irrep  = rep1.irrep;
                  irrep_lam_tmp.row         = rep1.row;
                  irrep_lam_tmp.mom         = mom1;
                  irrep_lam_tmp.mom_sq      = mom1_sq;
                  irrep_lam_tmp.lev         = level_1;
                  irrep_lam_tmp.two_lam     = get<0>(two_abs_lam);

                  irrep_tmp.Npt.push_back(irrep_lam_tmp);

                  if(mom_curr_sq) irrep_lam_tmp.irrep = LG_curr+rep_curr.irrep;
                  else irrep_lam_tmp.irrep  = rep_curr.irrep;
                  irrep_lam_tmp.row         = rep_curr.row;
                  irrep_lam_tmp.mom         = mom_curr;
                  irrep_lam_tmp.mom_sq      = mom_curr_sq;
                  irrep_lam_tmp.lev         = "";
                  irrep_lam_tmp.two_lam     = get<2>(two_abs_lam);

                  irrep_tmp.Npt.push_back(irrep_lam_tmp);

                  if(mom3_sq) irrep_lam_tmp.irrep = LG3+rep3.irrep;
                  else irrep_lam_tmp.irrep  = rep3.irrep;
                  irrep_lam_tmp.row         = rep3.row;
                  irrep_lam_tmp.mom         = mom3;
                  irrep_lam_tmp.mom_sq      = mom3_sq;
                  irrep_lam_tmp.lev         = level_3;
                  irrep_lam_tmp.two_lam     = get<1>(two_abs_lam);

                  irrep_tmp.Npt.push_back(irrep_lam_tmp);


                  Array1dO<Complex> tmp(num_matrix_elem);

                  for(int num = 0; num < num_matrix_elem; num++){tmp[num+1] = cmplx( Real(std::real(SEMBLE::toScalar(ENSEM::mean(prefactor[num])))) , Real(std::imag(SEMBLE::toScalar(ENSEM::mean(prefactor[num])))) );}
                    irrep_tmp.kfac = tmp;
                  irreps.push_back(irrep_tmp);
                  count++; 


                }

            }}
            //}
                                
            }}}
          }}}
      }


      else{continue;}
  }
}

    // plot to see the phase space explored  E vs Q2 plot
  {  
    ostringstream outfile; outfile << "phase_space.plot";
    std::ofstream file;
    file.open(outfile.str(), std::ios_base::app);
    file << "# ecm || ecm_err || qsq || qsq_err" << endl << endl << endl;
    
    file << "# npt_one" << endl;
    for(auto it = ecm1_qsq.begin(); it != ecm1_qsq.end();it++){
      file << SEMBLE::toScalar(ENSEM::mean(it->first)) << " " <<  sqrt(SEMBLE::toScalar(ENSEM::variance(it->first))) <<  " " << SEMBLE::toScalar(ENSEM::mean(it->second)) << 
              " " << sqrt(SEMBLE::toScalar(ENSEM::variance(it->second))) << endl;  
    }
    
    file << endl << endl << "# npt_three" << endl;
    for(auto it = ecm3_qsq.begin(); it != ecm3_qsq.end();it++){
      file << SEMBLE::toScalar(ENSEM::mean(it->first)) << " " <<  sqrt(SEMBLE::toScalar(ENSEM::variance(it->first))) <<  " " << SEMBLE::toScalar(ENSEM::mean(it->second)) << 
              " " << sqrt(SEMBLE::toScalar(ENSEM::variance(it->second))) << endl;  
    }
  }  



  //output xml
  //==============================
  //WRITE THE OUTPUT XML kfac_xml
  //==============================
  XMLFileWriter xml_out(out);
  write_xml_out(xml_out, npt, L, Xi, XiE, had, num_matrix_elem, matrix_type, matrix_name, print_zero,
                        compute_phase, elab_dir, make_redstar_xml, irreps , count );
  xml_out.close();
  
  cout << count << "\n" ;

  //==============================
  //WRITE THE OUTPUT REDSTAR XML
  //==============================

  if(make_redstar_xml){
    XMLFileWriter red_xml(redstar_xml);
    gen_redstar_xml(had, irreps, red_xml);
    red_xml.close();
  }

  
};