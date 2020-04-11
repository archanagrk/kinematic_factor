
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
    
  int npt;
  Array1dO<int> creation_op;
  Array1dO<int> smearedP;
  Array1dO<int> twoI;
  Array1dO<int> threeY;
  Array1dO<int> twoI_z;
  Array1dO<int> t_slice;
  std::vector<Array1dO<int>> mom;
  Array1dO<int> mom_tmp;
  int row_tmp;
  int lam_tmp;
  string irrep_tmp, name_tmp;
  XMLArray::Array<int> canon_mom_tmp;
  double psq_tmp;
  string lev_tmp;

  
  int pts;
  
      
      
  XMLReader xml_in(in);
       
   try{            
        
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
      read(xml_kf_in, "/kfac/hadron/had/elem["+std::to_string(i)+"]/name", name_tmp);
      if(i != 2){read(xml_kf_in, "/kfac/elem["+std::to_string(k)+"]/elem["+std::to_string(i)+"]/level", lev_tmp);}
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
        if(psq_tmp != 0){op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_" + lev_tmp + "_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_H"+ std::to_string(lam_tmp/2) + irrep_tmp, "", canon_mom_tmp);}
        else{op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_"+ irrep_tmp, "", canon_mom_tmp);}
      }
      else if(i == 2){
        if(psq_tmp != 0){op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_H"+ std::to_string(lam_tmp/2) + irrep_tmp, "", canon_mom_tmp);}
        else{op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_" + irrep_tmp, "", canon_mom_tmp);}        
      }
      else if(i == 3){
        if(psq_tmp != 0){op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_" + lev_tmp + "_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_H"+ std::to_string(lam_tmp/2) + irrep_tmp, "", canon_mom_tmp);}
        else{op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_"+ irrep_tmp, "", canon_mom_tmp);}        
      }

      key.npoint[i].irrep.op = op_tmp;

      }


    corrs.push_back(key);
    
    

  }
  


  xml_kf_in.close();


  
  write(xml_out, "NPointList", corrs);
  xml_out.close();

  
};
