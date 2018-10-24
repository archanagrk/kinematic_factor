#include "hadron_3pt.xml.h"


int main(int argc, char** argv){
    if( argc != 4 ){
        cerr << "xml_out <xml_in> <Ckfac.xml> <filename_xml_out> \n ";
        exit(1); }
    
    std::string in;  {istringstream a(argv[1]); a >> in;};
    std::string kfac;  {istringstream a(argv[2]); a >> kfac;};
    std::string out;  {istringstream a(argv[3]); a >> out;};
    
    int npt;
    Array1dO<string> creation_op;
    Array1dO<string> smearedP;
    Array1dO<int> twoI;
    Array1dO<int> threeY;
    Array1dO<int> twoI_z;
    Array1dO<int> t_slice;
    
    
    
    had_npt_layout label;
    
    
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
    
    
    xml_in.close();
    
    XMLFileWriter xml_out(out);
    write_had_layout(xml_out, "RedstarNPt", label);
    
   
    
    
   XMLReader xml_kfac(kfac);

    for(int i = 1; i <= 1666; i++){
        
        KeyHadronSUNNPartNPtCorr_t npt_irrep;
        
        npt_irrep.npoint.resize(npt);
        
        for(int j = 1; j <= npt; j++){
            
            npt_irrep.npoint[j].t_slice = t_slice[j];
            
            istringstream(creation_op[j]) >> std::boolalpha >> npt_irrep.npoint[j].irrep.creation_op;
            istringstream(smearedP[j]) >> std::boolalpha >> npt_irrep.npoint[j].irrep.smearedP;
            
            
            
            //npt_irrep.npoint[j].irrep.creation_op = creation_op[j];
            //npt_irrep.npoint[j].irrep.smearedP = smearedP[j];
            npt_irrep.npoint[j].irrep.flavor.twoI = twoI[j];
            npt_irrep.npoint[j].irrep.flavor.threeY = threeY[j];
            npt_irrep.npoint[j].irrep.flavor.twoI_z = twoI_z[j];
            
        
            read(xml_kfac, "/kfac/elem["+std::to_string(i)+"]/row"+std::to_string(j), npt_irrep.npoint[j].irrep.irrep_mom.row );
            read(xml_kfac, "/kfac/elem["+std::to_string(i)+"]/mom"+std::to_string(j), npt_irrep.npoint[j].irrep.irrep_mom.mom );
            
            npt_irrep.npoint[j].irrep.op.CGs.resize(0);
            npt_irrep.npoint[j].irrep.op.ops.resize(1);
            

            read(xml_kfac, "/kfac/elem["+std::to_string(i)+"]/irrep"+std::to_string(j), npt_irrep.npoint[j].irrep.op.ops[0].name);
            
            read(xml_kfac, "/kfac/elem["+std::to_string(i)+"]/mom"+std::to_string(j), npt_irrep.npoint[j].irrep.op.ops[0].mom_type);
            
            xml_kfac.close();
            
            
            
            
            //label.npointL[i].npoint[j].irrep.op.CGs[0]
            //label.npointL[i].npoint[j].irrep.op.ops[0].smear
            //mom type is wrong do it later
            
        }
        
        write(xml_out, "RedstarNPt", npt_irrep);
    }
    

    
/*

    

    for(int i = 1; i <= 1666; i++){
        for(int j = 1; j <= npt; j++){
            label.n.npoint[j].t_slice = t_slice[j];
            
            istringstream(creation_op[j]) >> std::boolalpha >> label.n.npoint[j].irrep.creation_op;
            istringstream(smearedP[j]) >> std::boolalpha >> label.n.npoint[j].irrep.smearedP;
            
            //label.npoint[j].irrep.creation_op = creation_op[j];
            //label.npoint[j].irrep.smearedP = smearedP[j];
            label.n.npoint[j].irrep.flavor.twoI = twoI[j];
            label.n.npoint[j].irrep.flavor.threeY = threeY[j];
            label.n.npoint[j].irrep.flavor.twoI_z = twoI_z[j];
            
            
            read(xml_kfac, "/kfac/elem["+std::to_string(i)+"]/row"+std::to_string(j), label.n.npoint[j].irrep.irrep_mom.row );
            read(xml_kfac, "/kfac/elem["+std::to_string(i)+"]/mom"+std::to_string(j), label.n.npoint[j].irrep.irrep_mom.mom );
            
            label.n.npoint[j].irrep.op.CGs.resize(0);
            label.n.npoint[j].irrep.op.ops.resize(1);
            
            
            read(xml_kfac, "/kfac/elem["+std::to_string(i)+"]/irrep"+std::to_string(j), label.n.npoint[j].irrep.op.ops[1].name);
            
            read(xml_kfac, "/kfac/elem["+std::to_string(i)+"]/mom"+std::to_string(j), label.n.npoint[j].irrep.op.ops[1].mom_type);
            
            xml_kfac.close();
            
            
            
            
            //label.npoint[j].irrep.op.CGs[0]
            //label.npoint[j].irrep.op.ops[0].smear
            //mom type is wrong do it later
            
        }

					  Hadron::KeyHadronSUNNPartNPtCorr_t key;
                      label.npointL.resize(1);
					  label.npointL.npoint.resize(1);
    


					  key.npoint[1].t_slice = ;
					  key.npoint[1].irrep.creation_op = false;
					  key.npoint[1].irrep.smearedP = true;
					  key.npoint[1].irrep.flavor.twoI = 1;
					  key.npoint[1].irrep.flavor.threeY = 3;
					  key.npoint[1].irrep.flavor.twoI = 1;
					  key.npoint[1].irrep.irrep_mom.row = rep1.row;
					  key.npoint[1].irrep.irrep_mom.mom = mom1;

					  KeyHadronSUNNPartIrrepOp_t irrep1;
					  irrep1.CGs.resize(0);
					  irrep1.ops.resize(1);
					  irrep1.ops[1]= KeyParticleOp_t("Kneg_proj0_p" + canon_mom1[0] + canon_mom1[1] + canon_mom2[2] + "_H0" + rep1, "", canon_mom1);
					  key.npoint[1].irrep.op = irrep1; 

					  key.npoint[2].t_slice = -3;
					  KeyHadronSUNNPartIrrepOp_t irrep2;
					  irrep2.ops[1] = KeyParticleOp_t("omegal_rhoxD0_J0__J1_" + "H0" + rep2, "", canon_mom2);
					  key.npoint[2].irrep.op = irrep2; 

					  key.npoint[3].t_slice = 0;
					  KeyHadronSUNNPartIrrepOp_t irrep3;
					  irrep3.CGs.resize(0);
					  irrep3.ops.resize(1);
					  //irrep3.ops[1]= KeyParticleOp_t("Kneg_proj0_p" + canon_mom1[0] + canon_mom1[1] + canon_mom2[2] + "_H0" + rep1, "", canon_mom3);
					  irrep3.ops[1] = KeyParticleOp_t("Kneg_rhoxD0_J0__J1_" + "H0" + rep3, "", canon_mom3);
					  key.npoint[3].irrep.op = irrep3; 

					  corrs.push_back(key); */
    
    
    
    
    
    
};

