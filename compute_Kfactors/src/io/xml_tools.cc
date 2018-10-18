#include "hadron_3pt.xml.h"

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
    XMLArrayWriter xml_ar(xml, label.npointL.size());
    push(xml_ar, "NpointList");
    for(int i=0; i < xml_ar.size(); ++i)
    {
        push(xml_ar);
        write(xml_ar, "NPoint", label.npointL[i]);
        pop(xml_ar);
    }
    
    pop(xml);
    pop(xml);
    
};


