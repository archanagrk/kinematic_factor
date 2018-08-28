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



//structures
struct irrep_label{
  string irrep;   /* irrep name */
  int row;        /* row of the irrep - 1 based */

  int twoS;       /* total spin of the two particle system (*2) */
  int ell;        /* orb ang mom of the two particle system     */
  int twoJ;       /* total ang mom subduced from (*2)           */
  int n;          /* embedding of spin-J into irrep "irrep" - 1 based    */

  int P1P2;       /* product of the intrinsic parities of the two particles */

  bool operator<(const irrep_label &rhs) const; /*usual map label problems */
  //bool operator!=(const irrep_label &rhs) const; /*usual map label problems */
};

#define PI 3.14159

#endif
