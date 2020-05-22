#ifndef __GENREDXML__
#define __GENREDXML__

#include "compute_kfactor_xml_read_write.h"

using namespace std;
using namespace ADAT;
using namespace ADATXML;
using namespace Hadron;

// fucntion to write out the xml that can be read by redstar
void gen_redstar_xml(vector<NPtCorr_t>& had, vector<NPtIrrepLam_t>& irreps, XMLWriter& red_xml);


#endif



