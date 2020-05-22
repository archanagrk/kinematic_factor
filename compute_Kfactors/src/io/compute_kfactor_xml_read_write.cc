#include "compute_kfactor_xml_read_write.h"

//==============================
//READ THE INPUT XML
//==============================
  //**********************************************************************************************************************
void read_xml_ini( XMLReader&  xml_in, int& npt, int& L, double& Xi, double& XiE, std::vector<NPtCorr_t>& had, int& num_matrix_elem, string& matrix_type, Array1dO<string>& matrix_name,
                          string& print_zero, string& compute_phase, string& elab_dir,bool& make_redstar_xml, string& redstar_xml ){
  
    NPtCorr_t had_tmp;

    try
    {
      read(xml_in, "/hadron/npt", npt);
      read(xml_in, "/hadron/L", L);
      read(xml_in, "/hadron/Xi", Xi);
      read(xml_in, "/hadron/XiError", XiE);


      for(int j =1;j<=npt;j++){
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/name", had_tmp.name);
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/twoJ", had_tmp.twoJ);
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/P", had_tmp.P);
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/ell", had_tmp.ell);
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/projected", had_tmp.projected); 
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/levels", had_tmp.levels);
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/Elab", had_tmp.elab);
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/max_mom", had_tmp.max_mom);
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/min_mom", had_tmp.min_mom);
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/canonical", had_tmp.canonical);    
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/omit_mom", had_tmp.omit_mom); 
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/flavor/threeY", had_tmp.flavor.threeY); 
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/flavor/twoI", had_tmp.flavor.twoI); 
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/flavor/twoIz", had_tmp.flavor.twoIz);         
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/smearedP", had_tmp.smearedP);  
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/creation_op", had_tmp.creation_op);  
      read(xml_in, "/hadron/had/elem["+std::to_string(j)+"]/t_slice", had_tmp.t_slice);      
      had.push_back(had_tmp);
      }
      read(xml_in, "/hadron/ElabDir", elab_dir);
      read(xml_in, "/hadron/MatrixType", matrix_type);
      read(xml_in, "/hadron/numMatrixElem", num_matrix_elem);
      matrix_name.resize(num_matrix_elem);
      read(xml_in, "/hadron/MatrixName", matrix_name);
      read(xml_in, "/hadron/PrintZero", print_zero);

      read(xml_in, "/hadron/phase", compute_phase);
      read(xml_in, "/hadron/makeRedstarXml", make_redstar_xml);
      if(make_redstar_xml)read(xml_in, "/hadron/redstarXml", redstar_xml);

    }
      catch( const string& error ){
      cerr << "Error reading hadron file : " << error << endl;
      }
};


//==============================
//WRITE THE OUTPUT XML
//==============================

  //**********************************************************************************************************************  
void write_ei( XMLWriter& xml, const std::string& path, const Eigen::Vector3d& input){
    XMLArray::Array<int> in_1(3); in_1[0] = input(0); in_1[1] = input(1); in_1[2] = input(2);
    write(xml , path , in_1);
    };

  //**********************************************************************************************************************

void write_irrep(XMLWriter&  xml_out, IrrepLam_t& irrep_lam){
    push(xml_out, "elem");
    write(xml_out, "irrep", irrep_lam.irrep);
    write(xml_out, "row", irrep_lam.row);
    write_ei(xml_out, "mom", irrep_lam.mom);
    write(xml_out, "psq", irrep_lam.mom_sq);
    write(xml_out, "level", irrep_lam.lev);
    write(xml_out, "absLam", irrep_lam.two_lam);  
    pop(xml_out);
}

  //**********************************************************************************************************************
  
void write_xml_out( XMLWriter&  xml_out, int& npt, int& L, double& Xi, double& XiE, std::vector<NPtCorr_t>& had, int& num_matrix_elem, string& matrix_type, Array1dO<string>& matrix_name,
                          string& print_zero, string& compute_phase, string& elab_dir,bool& make_redstar_xml, vector<NPtIrrepLam_t>& irreps, int& count ){

  push(xml_out, "kfac");

  write(xml_out, "L", L);
  write(xml_out, "Xi", Xi);
  write(xml_out, "XiError", XiE);
  write(xml_out, "ElabFilesIn", had[2].elab);
  write(xml_out, "ElabFilesOut", had[0].elab);
  write(xml_out, "kfacFile", matrix_name);
      
  for(int i =0;i<count;i++){
    push(xml_out, "elem");

    for(int k=0; k<npt; k++){
      write_irrep(xml_out, irreps[i].Npt[k]);
    }
    pop(xml_out);

    write(xml_out, "FF", irreps[i].kfac );
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



 };


  //**********************************************************************************************************************  

