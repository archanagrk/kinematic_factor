#include "gen_redstar_xml.h"


void write_had_layout( XMLWriter& xml, const std::string& path, const had_npt_layout& label){
    
    push(xml, path);
    push(xml, "Param");
    write(xml, "version", 10);
    write(xml, "diagnostic_level", 0);
    write(xml, "Nt_corr", label.Nt_corr);
    write(xml, "t_origin", label.t_origin);
    write(xml, "bc_spec", label.bc_spec);
    write(xml, "convertUDtoL", label.convertUDtoL);
    write(xml, "convertUDtoS", label.convertUDtoS);
    write(xml, "average_1pt_diagrams", label.average_1pt_diagrams);
    write(xml, "zeroUnsmearedGraphsP", label.zeroUnsmearedGraphsP);
    write(xml, "ensemble", label.ensemble);
    push(xml, "Layout");
    write(xml, "decayDir", label.decayDir);
    write(xml, "lattSize", label.lattSize);
    pop(xml);

};

/*
void write_had_nptL( XMLWriter& xml, const std::string& path, const had_npt_Llabel& label){
    
    
    push(xml, path);
    push(xml, "Param");
    write(xml, "version", 10);
    write(xml, "diagnostic_level", 2);
    write(xml, "Nt_corr", label.Nt_corr);
    write(xml, "t_origin", label.t_origin);
    write(xml, "bc_spec", label.bc_spec);
    write(xml, "convertUDtoL", label.convertUDtoL);
    write(xml, "convertUDtoS", label.convertUDtoS);
    write(xml, "average_1pt_diagrams", label.average_1pt_diagrams);
    write(xml, "zeroUnsmearedGraphsP", label.zeroUnsmearedGraphsP);
    write(xml, "ensemble", label.ensemble);
    push(xml, "Layout");
    write(xml, "decayDir", label.decayDir);
    write(xml, "latticeSize", label.latticeSize);
    pop(xml);
    //XMLArrayWriter xml_ar(xml, 1);  //label.npointL.size()
    //Hadron::write(xml_ar, "NPoint", label.n); //npointL
    //push(xml_ar, "NpointList");
    //for(int i=0; i < 1; ++i)
    //{
      //  push(xml_ar);
      //  Hadron::write(xml_ar, "NPoint", label.npointL); //[i] missing
       // pop(xml_ar);
   // }
    
    //pop(xml);
    
    write(xml, "NPoint", label.npointL);
    pop(xml);
    pop(xml);
    
}; */

void write_ei( XMLWriter& xml, const std::string& path, const Eigen::Vector3d& input){
    XMLArray::Array<int> in_1(3); in_1[0] = input(0); in_1[1] = input(1); in_1[2] = input(2);
    write(xml , path , in_1);
    
};

void write_db_keys (XMLWriter& xml, const std::string& path, const db& label){
    push(xml, path);
  write(xml, "proj_op_xmls", label.proj_op_xmls);
  write(xml, "corr_graph_db", label.corr_graph_db);
  write(xml, "noneval_graph_xml", label.noneval_graph_xml);
  write(xml, "smeared_hadron_node_xml", label.smeared_hadron_node_xml);
  write(xml, "unsmeared_hadron_node_xml", label.unsmeared_hadron_node_xml);
  write(xml, "hadron_npt_graph_db", label.hadron_npt_graph_db);
  write(xml, "hadron_node_dbs", label.hadron_node_dbs);
  write(xml, "output_db",label.output_db);
    pop(xml);


};
