#ifndef __HADRON3PT__
#define __HADRON3PT__

#include <Eigen/Dense>

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


struct had_npt_layout{
    
    /* Layout Options */
    
    
    int Nt_corr;
    int t_origin;
    int bc_spec;  /* boundary condition -1 */
    bool convertUDtoL;
    bool convertUDtoS;
    bool average_1pt_diagrams;
    bool zeroUnsmearedGraphsP;
    string ensemble;
    int decayDir;
    Array1dO<int> lattSize;
    
    /* Contains the ensemble of Npoint Fncts */
    
    //Array1dO<Hadron::KeyHadronSUNNPartNPtCorr_t> npointL;
    //Hadron::KeyHadronSUNNPartNPtCorr_t n;
    
};

struct db{

  Array1dO<std::string> proj_op_xmls;
  std::string corr_graph_db;
  std::string noneval_graph_xml;
  std::string smeared_hadron_node_xml;
  std::string unsmeared_hadron_node_xml;
  std::string hadron_npt_graph_db;
  Array1dO<std::string>  hadron_node_dbs;
  std::string output_db;
};

struct hadron{

  int twoJ;
  int P;
  double msq;
  double max_mom;

};

//function declaration

  void write_had_layout( XMLWriter& xml, const std::string& path, const had_npt_layout& label);
  void write_ei( XMLWriter& xml, const std::string& path, const Eigen::Vector3d& input);
  void write_db_keys( XMLWriter& xml, const std::string& path, const db& label);

#endif



