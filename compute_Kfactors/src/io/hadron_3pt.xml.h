#ifndef __HADRON3PT__
#define __HADRON3PT__


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

using namespace std;
using namespace ADAT;
using namespace ADATXML;
using namespace Hadron;


//structs


struct had_npt_Llabel{
    
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
    Array1dO<int> latticeSize;
    
    /* Contains the ensemble of Npoint Fncts */
    
    Array1dO<Hadron::KeyHadronSUNNPartNPtCorr_t> npointL;
    
};

//function declaration

  void write_had_nptL( XMLWriter& xml, const std::string& path, const had_npt_Llabel& label);

#endif



