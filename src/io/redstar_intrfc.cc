
/* main code to return the irreps to compute */

#include "../lib/kfactors.h"
#include "hadron_3pt.xml.h"

int main(int argc, char** argv){
  if( argc != 13 ){
    cerr << "get_irreps <twoJ1> <P1> <m1sq> <twoJ_curr> <P_curr>  <twoJ3> <P3> <m3sq> <max_mom1> <max_mom2> <max_mom3> <xml_in>\n ";
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
  std::string in;  {istringstream a(argv[12]); a >> in;};

  complex<double> Coeff;

  if((two_J1==0) && (P1==-1) && (two_J2==2) && (P2==-1) && (two_J3==2) && (P3==-1)){

    std::list<Hadron::KeyHadronSUNNPartNPtCorr_t> corrs;
    had_npt_layout label;
    
    Vector3d mom1(3,1);    Vector3d mom3(3,1);    Vector3d mom_curr(3,1);
    
      int npt;
      Array1dO<string> creation_op;
      Array1dO<string> smearedP;
      Array1dO<int> twoI;
      Array1dO<int> threeY;
      Array1dO<int> twoI_z;
      Array1dO<int> t_slice;
      
      
      XMLReader xml_in(in);
      
      
    
      
      read(xml_in, "/doc/config/Nt_corr", label.Nt_corr);
      read(xml_in, "/doc/config/t_origin", label.t_origin);
      read(xml_in, "/doc/config/bc_spec", label.bc_spec);
      read(xml_in, "/doc/config/convertUDtoL", label.convertUDtoL);
      read(xml_in, "/doc/config/convertUDtoS", label.convertUDtoS);
      read(xml_in, "/doc/config/average_1pt_diagrams", label.average_1pt_diagrams);
      read(xml_in, "/doc/config/zeroUnsmearedGraphsP", label.zeroUnsmearedGraphsP);
      read(xml_in, "/doc/config/ensemble", label.ensemble);
      read(xml_in, "/doc/config/decayDir", label.decayDir);
      read(xml_in, "/doc/config/latticeSize", label.latticeSize);
      
      read(xml_in, "/doc/Npts", npt);
      read(xml_in, "/doc/creation_op", creation_op);
      read(xml_in, "/doc/smearedP", smearedP);
      read(xml_in, "/doc/t_slice", t_slice);
      
      read(xml_in, "/doc/fl/twoI", twoI);
      read(xml_in, "/doc/fl/threeY", threeY);
      read(xml_in, "/doc/fl/twoI_z", twoI_z);
      
      
      //xml_in.close();
      
      XMLFileWriter xml_out("kfac.xml");
      write_had_layout(xml_out, "RedstarNPt", label);

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

                    qp  << (sqrt(m1_sq+mom1_sq)+sqrt(m3_sq+mom3_sq)),(i+l),(j+m),(k+n);
                    qm  << (sqrt(m1_sq+mom1_sq)-sqrt(m3_sq+mom3_sq)),(l-i),(m-j),(n-k);

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




                                        Coeff = KinematicFactor(qp,qm,Sub1.sum,SubCurr.sum,Sub3.sum);
                                        if(std::real(Coeff) || std::imag(Coeff)){
                                          cout << mom1.transpose() << rep1.irrep << "["<< rep1.row <<"]" << "\n";
                                          cout << mom_curr.transpose()  << rep_curr.irrep << "["<< rep_curr.row <<"]"<< "\n";
                                          cout << mom3.transpose()  << rep3.irrep << "["<< rep3.row <<"]"<< "\n";
                                        //cout << "The abs_factor is:" << pow(std::real(Coeff),2)+pow(std::imag(Coeff),2) << "\n";
                                          cout << "The factor is:" << Coeff << "\n";

					                                Hadron::KeyHadronSUNNPartNPtCorr_t key;
					                                key.npoint.resize(npt);

					                                key.npoint[1].t_slice = t_slice[1];
					                                key.npoint[1].irrep.creation_op = false;
					                                key.npoint[1].irrep.smearedP = true;
					                                key.npoint[1].irrep.flavor.twoI = 1;
					                                key.npoint[1].irrep.flavor.threeY = 3;
					                                key.npoint[1].irrep.flavor.twoI_z = 1;
					                                key.npoint[1].irrep.irrep_mom.row = rep1.row;
					                                key.npoint[1].irrep.irrep_mom.mom = KfUt::toArray(mom1);

					                                KeyHadronSUNNPartIrrepOp_t op1;
					                                op1.CGs.resize(0);
					                                op1.ops.resize(1);

                                                    XMLArray::Array<int> canon_mom1 = Hadron::canonicalOrder(KfUt::toArray(mom1));
					                                op1.ops[1]= KeyParticleOp_t("Kneg_proj0_p" + std::to_string(canon_mom1[0]) + std::to_string(canon_mom1[1]) + std::to_string(canon_mom1[2]) + "_H0" + rep1.irrep, "", canon_mom1);
					                                key.npoint[1].irrep.op = op1; 

					                                key.npoint[2].t_slice = t_slice[2];
                                                    key.npoint[2].irrep.creation_op = false;
                                                    key.npoint[2].irrep.smearedP = false;
                                                    key.npoint[2].irrep.flavor.twoI = 2;
                                                    key.npoint[2].irrep.flavor.threeY = 0;
                                                    key.npoint[2].irrep.flavor.twoI = 0;
                                                    key.npoint[2].irrep.irrep_mom.row = rep_curr.row;
                                                    key.npoint[2].irrep.irrep_mom.mom = KfUt::toArray(mom_curr);
                                            
                                            
					                                KeyHadronSUNNPartIrrepOp_t op2;
                                                    op2.CGs.resize(0);
                                                    op2.ops.resize(1);
                                            
                                                    XMLArray::Array<int> canon_mom2 = Hadron::canonicalOrder(KfUt::toArray(mom_curr));
					                                op2.ops[1] = KeyParticleOp_t("omegal_rhoxD0_J0__J1_H0" + rep_curr.irrep , "", canon_mom2);
					                                key.npoint[2].irrep.op = op2; 

					                                key.npoint[3].t_slice = t_slice[3];
                                                    key.npoint[3].irrep.creation_op = true;
                                                    key.npoint[3].irrep.smearedP = true;
                                                    key.npoint[3].irrep.flavor.twoI = 2;
                                                    key.npoint[3].irrep.flavor.threeY = 3;
                                                    key.npoint[3].irrep.flavor.twoI = 2;
                                                    key.npoint[3].irrep.irrep_mom.row = rep3.row;
                                                    key.npoint[3].irrep.irrep_mom.mom = KfUt::toArray(mom3);
                                
					                                KeyHadronSUNNPartIrrepOp_t op3;
					                                op3.CGs.resize(0);
					                                op3.ops.resize(1);
                                                    int helicity = Sub3.two_hel/2;
					                                //irrep3.ops[1]= KeyParticleOp_t("Kneg_proj0_p" + canon_mom1[0] + canon_mom1[1] + canon_mom2[2] + "_H0" + rep1, "", canon_mom3);
                                                    XMLArray::Array<int> canon_mom3 = Hadron::canonicalOrder(KfUt::toArray(mom3));
                                                    op3.ops[1] = KeyParticleOp_t("Kneg_rhoxD0_J0__J1_H"+ std::to_string(helicity) + rep3.irrep, "", canon_mom3);
					                                key.npoint[3].irrep.op = op3;

					                                corrs.push_back(key);
                                        }
                      }}}
                   }}}
                  }
                 }
             }}}
           }
         }}}

    write(xml_out, "NPointList", corrs);
    pop(xml_out);
    
    db db_lab;
      
      read(xml_in, "/doc/DBFiles/proj_op_xmls", db_lab.proj_op_xmls);
      read(xml_in, "/doc/DBFiles/corr_graph_db", db_lab.corr_graph_db);
      read(xml_in, "/doc/DBFiles/noneval_graph_xml", db_lab.noneval_graph_xml);
      read(xml_in, "/doc/DBFiles/smeared_hadron_node_xml", db_lab.smeared_hadron_node_xml);
      read(xml_in, "/doc/DBFiles/unsmeared_hadron_node_xml", db_lab.unsmeared_hadron_node_xml);
      read(xml_in, "/doc/DBFiles/hadron_npt_graph_db", db_lab.hadron_npt_graph_db);
      read(xml_in, "/doc/DBFiles/hadron_node_dbs", db_lab.hadron_node_dbs);
      read(xml_in, "/doc/DBFiles/output_db", db_lab.output_db);

      
      xml_in.close();
      
    
      write_db_keys(xml_out, "DBFiles", db_lab);
      
      
    xml_out.close();
    
  }
};
