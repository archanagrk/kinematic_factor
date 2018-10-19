#include "hadron_3pt.xml.h"


int main(int argc, char** argv){
    if( argc != 3 ){
        cerr << "xml_out <xml_in_file> <filename_xml_out> \n ";
        exit(1); }
    
    std::string in;  {istringstream a(argv[1]); a >> in;};
    std::string out;  {istringstream a(argv[2]); a >> out;};
    
    had_npt_Llabel label;
    
    
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
    
    
    
    
    
    XMLFileWriter xml_out(out);
    
    write_had_nptL(xml_out, "RedstarNPt", label);
    
    
    
};

