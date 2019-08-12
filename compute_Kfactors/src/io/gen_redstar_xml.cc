
/* main code to return the irreps to compute */

#include "gen_redstar_xml.h"

int main(int argc, char** argv){
  if( argc != 4 ){
    cerr << "create redstar xml <xml_in> <xml_kf_in> <out_xml>\n ";
    exit(1); }


  std::string in;  {istringstream a(argv[1]); a >> in;};
  std::string kf_in;  {istringstream a(argv[2]); a >> kf_in;};
  std::string out;  {istringstream a(argv[3]); a >> out;};



  std::list<Hadron::KeyHadronSUNNPartNPtCorr_t> corrs;
  had_npt_layout label;
  
    
  int npt;
  Array1dO<int> creation_op;
  Array1dO<int> smearedP;
  Array1dO<int> twoI;
  Array1dO<int> threeY;
  Array1dO<int> twoI_z;
  Array1dO<int> t_slice;
  std::vector<Array1dO<int>> mom;
  //std::vector<int> row;
  //std::vector<int> lam;
  //std::vector<string> irrep;
  Array1dO<int> mom_tmp;
  int row_tmp;
  int lam_tmp;
  string irrep_tmp;
  XMLArray::Array<int> canon_mom_tmp;
  double psq_tmp;
  //std::vector<XMLArray::Array<int>> canon_mom;

  
  int pts;
  
      
      
  XMLReader xml_in(in);
      
      
    
   try{    
    read(xml_in, "/doc/config/Nt_corr", label.Nt_corr);
    read(xml_in, "/doc/config/t_origin", label.t_origin);
    read(xml_in, "/doc/config/bc_spec", label.bc_spec);
    read(xml_in, "/doc/config/convertUDtoL", label.convertUDtoL);
    read(xml_in, "/doc/config/convertUDtoS", label.convertUDtoS);
    read(xml_in, "/doc/config/average_1pt_diagrams", label.average_1pt_diagrams);
    read(xml_in, "/doc/config/zeroUnsmearedGraphsP", label.zeroUnsmearedGraphsP);
    read(xml_in, "/doc/config/ensemble", label.ensemble);
    read(xml_in, "/doc/config/decayDir", label.decayDir);
    read(xml_in, "/doc/config/lattSize", label.lattSize);
        
        
    read(xml_in, "/doc/Npts", npt);
    read(xml_in, "/doc/creation_op", creation_op);
    read(xml_in, "/doc/smearedP", smearedP);
    read(xml_in, "/doc/t_slice", t_slice);
        
    read(xml_in, "/doc/fl/twoI", twoI);
    read(xml_in, "/doc/fl/threeY", threeY);
    read(xml_in, "/doc/fl/twoI_z", twoI_z);
   }

   catch( const string& error ){
    cerr << "Error reading input file : " << error << endl;
    }
    
  //xml_in.close();
  
  
  Hadron::KeyHadronSUNNPartNPtCorr_t key;
  key.npoint.resize(npt);
  KeyHadronSUNNPartIrrepOp_t op_tmp;




  XMLFileWriter xml_out(out);
  write_had_layout(xml_out, "RedstarNPt", label);


  XMLReader xml_kf_in(kf_in);
  
  
  try{read(xml_kf_in, "/kfac/pts", pts);}
  catch( const string& error ){cerr << "Error reading kfac input file : " << error << endl;}
  
  
  for(int k =1;k <= pts;k++){
    for(int i = 1; i <= npt; i++ ){
      
      
      
     try{
      read(xml_kf_in, "/kfac/elem["+std::to_string(k)+"]/elem["+std::to_string(i)+"]/irrep", irrep_tmp);
      read(xml_kf_in, "/kfac/elem["+std::to_string(k)+"]/elem["+std::to_string(i)+"]/row", row_tmp);
      read(xml_kf_in, "/kfac/elem["+std::to_string(k)+"]/elem["+std::to_string(i)+"]/mom", mom_tmp );
      read(xml_kf_in, "/kfac/elem["+std::to_string(k)+"]/elem["+std::to_string(i)+"]/absLam", lam_tmp);
      read(xml_kf_in, "/kfac/elem["+std::to_string(k)+"]/elem["+std::to_string(i)+"]/psq", psq_tmp);
     }
     catch( const string& error ){
      cerr << "Error reading kfac input file : " << error << endl;
     }
      

      key.npoint[i].t_slice = t_slice[i];
      key.npoint[i].irrep.creation_op = creation_op[i]==1?true:false;
      key.npoint[i].irrep.smearedP = smearedP[i]==1?true:false;

      key.npoint[i].irrep.flavor.twoI = twoI[i];
      key.npoint[i].irrep.flavor.threeY = threeY[i];
      key.npoint[i].irrep.flavor.twoI_z = twoI_z[i];
      

      key.npoint[i].irrep.irrep_mom.row = row_tmp;
      // if(i != 2){mom_tmp[1] = -mom_tmp[1]; mom_tmp[2] = -mom_tmp[2]; mom_tmp[3] = -mom_tmp[3]; cout << "curr";}

      key.npoint[i].irrep.irrep_mom.mom = KfUt::ToArray::toArray(mom_tmp);
      
      canon_mom_tmp = Hadron::canonicalOrder(KfUt::ToArray::toArray(mom_tmp));
      
      op_tmp.CGs.resize(0);
      op_tmp.ops.resize(1);

      
      if(i == 1){
        if(psq_tmp != 0){op_tmp.ops[1]= KeyParticleOp_t("pion_proj0_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_H"+ std::to_string(lam_tmp/2) + irrep_tmp, "", canon_mom_tmp);}
        else{op_tmp.ops[1]= KeyParticleOp_t("pion_proj0_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_"+ irrep_tmp, "", canon_mom_tmp);}
      }
      else if(i == 2){
        if(psq_tmp != 0){op_tmp.ops[1]= KeyParticleOp_t("omegal_rhoxD0_J0__J1_H"+ std::to_string(lam_tmp/2) + irrep_tmp, "", canon_mom_tmp);}
        else{op_tmp.ops[1]= KeyParticleOp_t("omegal_rhoxD0_J0__J1_" + irrep_tmp, "", canon_mom_tmp);}        
      }
      else if(i == 3){
        if(psq_tmp != 0){op_tmp.ops[1]= KeyParticleOp_t("rho_proj0_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_H"+ std::to_string(lam_tmp/2) + irrep_tmp, "", canon_mom_tmp);}
        else{op_tmp.ops[1]= KeyParticleOp_t("rho_proj0_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_"+ irrep_tmp, "", canon_mom_tmp);}        
      }

      key.npoint[i].irrep.op = op_tmp;

      }


    //canon_mom_tmp = Hadron::canonicalOrder(mom_tmp);
    //irrep.push_back(irrep_tmp); row.push_back(row_tmp); mom.push_back(mom_tmp); lam.push_back(lam_tmp); canon_mom.push_back(lam_tmp);


    corrs.push_back(key);
    
    

  }
  






 /* if(mom_curr_sq != 0){op2.ops[1] = KeyParticleOp_t("omegal_rhoxD0_J0__J1_H"+ std::to_string(helicity_curr) + LG_curr + rep_curr.irrep , "", canon_mom2);}
  else{op2.ops[1] = KeyParticleOp_t("omegal_rhoxD0_J0__J1_"  + rep_curr.irrep , "", canon_mom2);}
  key.npoint[2].irrep.op = op2;



  if(mom3_sq != 0){op3.ops[1] = KeyParticleOp_t("Kneg_rhoxD0_J0__J1_H"+ std::to_string(helicity3) + LG3+ rep3.irrep, "", canon_mom3);}
  else{op3.ops[1] = KeyParticleOp_t("Kneg_rhoxD0_J0__J1_" + rep3.irrep, "", canon_mom3);}
  key.npoint[3].irrep.op = op3; */



  xml_kf_in.close();


  
  write(xml_out, "NPointList", corrs);
  pop(xml_out);

  db db_lab;

 try{
  read(xml_in, "/doc/DBFiles/proj_op_xmls", db_lab.proj_op_xmls);
  read(xml_in, "/doc/DBFiles/corr_graph_db", db_lab.corr_graph_db);
  read(xml_in, "/doc/DBFiles/noneval_graph_xml", db_lab.noneval_graph_xml);
  read(xml_in, "/doc/DBFiles/smeared_hadron_node_xml", db_lab.smeared_hadron_node_xml);
  read(xml_in, "/doc/DBFiles/unsmeared_hadron_node_xml", db_lab.unsmeared_hadron_node_xml);
  read(xml_in, "/doc/DBFiles/hadron_npt_graph_db", db_lab.hadron_npt_graph_db);
  read(xml_in, "/doc/DBFiles/hadron_node_dbs", db_lab.hadron_node_dbs);
  read(xml_in, "/doc/DBFiles/output_db", db_lab.output_db);
 }

 catch( const string& error ){
  cerr << "Error reading input file : " << error << endl;
  }


  xml_in.close();
  

  write_db_keys(xml_out, "DBFiles", db_lab);

  pop(xml_out);


  xml_out.close();

  
};
