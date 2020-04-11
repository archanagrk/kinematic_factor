#ifndef __SUBDPOLVEC__
#define __SUBDPOLVEC__


//dependency
#include "subduction.h"
#include "pol_vec.h"


using namespace Subd;


//**********************************************************************************************************************

//functions
namespace SubdPol
{

/* Multiplies the subduced operators with the respective polarization tensors for J=1 */

sub_hel Subduce_all(double &mom_sq, double &mass_sq, int &twoJ, const irrep_label &irrep, const string &little_group, double R1_phi, double R1_theta, double R1_psi);


/* Map of helicity and subduction coefficients */
map<int, Eigen::MatrixXcd> Subduce_with_pol(double &mom_sq, double &mass_sq, int &twoJ, const irrep_label &irrep, const string &little_group, double R1_phi, double R1_theta, double R1_psi);

} // namespace SubPol


struct sub_hel
{
  Eigen::MatrixXcd sum;
  int two_hel;
};

#endif
