#ifndef __SUBDUCTION__
#define __SUBDUCTION__



//std lib
#include <vector>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <complex>
#include <utility>
#include <algorithm>
#include <cmath>
#include "math.h"
#include <stdio.h>
#include <Eigen/Dense>


//adat
#include "io/adat_io.h"
#include "io/adat_xml_group_reader.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include <adat/handle.h>
#include "hadron/clebsch.h"
#include "hadron/subduce_tables_factory.h"
#include "hadron/subduce_tables_oct_factory.h"
#include "hadron/subduce_tables_lg_factory.h"
#include "hadron/cgc_irrep_mom_auto.h"
#include "hadron/irrep_util.h" //adat_devel
#include "ensem/ensem.h"


using namespace std;
using namespace ADAT;
using namespace ADATXML;
using namespace Util;
using namespace ENSEM;
using namespace Hadron;
using namespace Eigen;

//dependency
#include "pol_vec.h"

//structures

struct irrep_label{
  string irrep;   /* irrep name */
  int row;        /* row of the irrep - 1 based */

  int twoJ;       /* total ang mom subduced from (*2)           */
  int n;          /* embedding of spin-J into irrep "irrep" - 1 based    */

  int P;       /* product of the intrinsic parities of the two particles */

  bool operator<(const irrep_label &rhs) const; /*usual map label problems */
  //bool operator!=(const irrep_label &rhs) const; /*usual map label problems */
};


//functions
namespace Subd{

  map< int, complex<double> > subduce_lg_boson(const irrep_label& irrep, const string& little_group);
  map< int, complex<double> > subduce_lg_fermion(const irrep_label& irrep, const string& little_group);
  map< int, complex<double> > subduce_oct(const irrep_label& irrep);
  MatrixXcd Subduce_all(double& mom_sq, double& mass_sq,  int& twoJ, const irrep_label& irrep,
                                                   const string& little_group,double R1_phi, double R1_theta, double R1_psi);
  int find_n_subduced_embeddings(const string& group, const string& irrep, int twoJ, int eta_tilde);
  const double PI = (atan(double(1)) * double(4.0));
  map< int, Eigen::MatrixXcd > Subduce_with_phases(double& mom_sq, double& mass_sq,  int& twoJ, const irrep_label& irrep,
                                                         const string& little_group,double R1_phi, double R1_theta, double R1_psi);
    
}




//Utility functions

/* J's for string stream Jo's code */

namespace {

string J_name( int twoJ ){
  stringstream ss;

  if( twoJ % 2 ){ //fermion
    ss << twoJ << "o2";
  }

  else{ //boson
    ss << (twoJ/2);
  }

  return ss.str();
  }

/* For sign of eta tilde  */


string sign(int x){

  if(x == -1){ return "-"; }
  else{ return "+"; }
  }

}



#endif
