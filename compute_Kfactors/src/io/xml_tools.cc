#include "gen_redstar_xml.h"


void write_ei( XMLWriter& xml, const std::string& path, const Eigen::Vector3d& input){
    XMLArray::Array<int> in_1(3); in_1[0] = input(0); in_1[1] = input(1); in_1[2] = input(2);
    write(xml , path , in_1);
    
};

