
/* fucntion to write out the xml that can be read by redstar */

#include "gen_redstar_xml.h"

void gen_redstar_xml(vector<NPtCorr_t>& had, vector<NPtIrrepLam_t>& irreps, XMLWriter& red_xml){

  std::list<Hadron::KeyHadronSUNNPartNPtCorr_t> corrs;
  int pts = irreps.size();
  int npt = had.size();

  std::vector<Array1dO<int>> mom;


  int lam_tmp;
  string irrep_tmp, name_tmp;
  XMLArray::Array<int> canon_mom_tmp;
  double psq_tmp;
  string lev_tmp;
  bool projected;


\
  
  Hadron::KeyHadronSUNNPartNPtCorr_t key;
  key.npoint.resize(npt);
  KeyHadronSUNNPartIrrepOp_t op_tmp;

  
  
  for(int k =1;k <= pts;k++){
    for(int i = 1; i <= npt; i++ ){
      
      key.npoint[i].t_slice           = had[i-1].t_slice;
      key.npoint[i].irrep.creation_op = had[i-1].creation_op;
      key.npoint[i].irrep.smearedP    = had[i-1].smearedP;

      key.npoint[i].irrep.flavor.twoI   = had[i-1].flavor.twoI;
      key.npoint[i].irrep.flavor.threeY = had[i-1].flavor.threeY;
      key.npoint[i].irrep.flavor.twoI_z = had[i-1].flavor.twoIz;
      

      key.npoint[i].irrep.irrep_mom.row = irreps[k-1].Npt[i-1].row;
      key.npoint[i].irrep.irrep_mom.mom = KfUt::ToArray::toArray(irreps[k-1].Npt[i-1].mom);
      
      name_tmp      = had[i-1].name;
      canon_mom_tmp = Hadron::canonicalOrder(KfUt::ToArray::toArray(irreps[k-1].Npt[i-1].mom));
      irrep_tmp     = irreps[k-1].Npt[i-1].irrep;
      lam_tmp       = irreps[k-1].Npt[i-1].two_lam;
      psq_tmp       = irreps[k-1].Npt[i-1].mom_sq;
      lev_tmp       = irreps[k-1].Npt[i-1].lev;
      projected     = had[i-1].projected;


      
      op_tmp.CGs.resize(0);
      op_tmp.ops.resize(1);

      
      if(projected){
        if(psq_tmp != 0){op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_" + lev_tmp + "_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_H"+ std::to_string(lam_tmp/2) + irrep_tmp, "", canon_mom_tmp);}
        else{op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_" + lev_tmp + "_p" + std::to_string(canon_mom_tmp[0]) + std::to_string(canon_mom_tmp[1]) + std::to_string(canon_mom_tmp[2]) + "_"+ irrep_tmp, "", canon_mom_tmp);}
      }
      else{
        if(psq_tmp != 0){op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_H"+ std::to_string(lam_tmp/2) + irrep_tmp, "", canon_mom_tmp);}
        else{op_tmp.ops[1]= KeyParticleOp_t(name_tmp + "_" + irrep_tmp, "", canon_mom_tmp);}        
      }

      key.npoint[i].irrep.op = op_tmp;

      }


    corrs.push_back(key);
    
    

  }
  

  write(red_xml, "NPointList", corrs);

  
};
