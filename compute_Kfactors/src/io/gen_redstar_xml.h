#ifndef __HADRON3PT__
#define __HADRON3PT__

#include </usr/local/Eigen/Dense>

//adat
#include "io/adat_io.h"
#include "io/adat_xml_group_reader.h"
#include "hadron/hadron_sun_npart_npt_corr.h"
#include "hadron/particle_op.h"
#include "hadron/hadron_sun_npart_irrep.h"
#include "hadron/hadron_sun_npart_irrep_op.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include <adat/handle.h>
#include "hadron/cgc_irrep_mom.h"
#include "hadron/cgc_su3.h"

//kfac includes
#include "../lib/kfac_utils.h"

using namespace std;
using namespace ADAT;
using namespace ADATXML;
using namespace Hadron;


//structs


struct hadron{

  string name;
  Array1dO<string> levels;
  int twoJ;
  int P;
  int ell;
  double max_mom;
  ADAT::Array1dO<string> elab;

};

//function declaration
  void write_ei( XMLWriter& xml, const std::string& path, const Eigen::Vector3d& input);


#endif



