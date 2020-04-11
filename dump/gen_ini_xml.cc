#include "gen_redstar_xml.h"

int main(int argc, char** argv){
  if( argc != 13 ){
    cerr << "get_irreps <twoJ1> <P1> <m1sq> <twoJ_curr> <P_curr>  <twoJ3> <P3> <m3sq> <max_mom1> <max_mom2> <max_mom3> <out_xml>\n ";
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
  string out_xml; {istringstream a(argv[12]); a >> out_xml;}
    
  XMLFileWriter xml_out(out_xml);
  write(xml_out, "RedstarNPt", hadron);

};
