/* main code to return the irreps to compute */

#include "lib/compute_three_point_prefactor.h"


int main(int argc, char** argv){

  //=================
  KFactorEnv::registerAll();
  KFactorEnv::registerAll();
  //=================

  if( argc != 4 ){
    cerr << "get_prefactors <hadron_xml> <output_xml> <npt>\n ";
    exit(1); }
    
  std::string in;      {istringstream a(argv[1]); a >> in;};
  std::string out;     {istringstream a(argv[2]); a >> out;};
  int pts;             {istringstream a(argv[3]); a >> pts;};  

    
  if(pts == 3) compute_three_point_prefactor(in,out);
  else cout << "No support for " << pts << " yet" << endl;

  
};

