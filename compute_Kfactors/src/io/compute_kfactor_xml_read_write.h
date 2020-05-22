#ifndef __CKFRWXML__
#define __CKFRWXML__

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

//structs
struct IrrepLam_t{
  string irrep;                      // irrep
  int row;                           // row of the irrep
  Eigen::Vector3d mom;               // momenta
  double mom_sq;                     // square of momenta
  string lev;                        // level - proj0, proj1 etc
  int two_lam;                       // 2*helicity
};

struct NPtIrrepLam_t{
    vector<IrrepLam_t> Npt;          // Npt - for threept Npt = 3
    Array1dO<Complex> kfac;          // the value of prefac to for the particular combination of irreps
};


struct flavour{
  int threeY;                         // three * hypercharge
  int twoIz;                          // two * isospin z component
  int twoI;                           // two * isospin
};

struct NPtCorr_t{

  string name;                        //name e.g. Kneg
  Array1dO<string> levels;            // levels if projected proj0, proj1 etc
  int twoJ;                           // 2*continuum spin
  int P;                              // parity
  int ell;                            // L
  ADAT::Array1dO<string> elab;        // elab files

  double max_mom;                     // max mom of particle
  double min_mom;                     // min mom of particle
  bool canonical;                     // mom only in canonical direction
  ADAT::Array1dO<ADAT::Array1dO<int>> omit_mom;    // omit these moms

  flavour flavor;                     // flavor structure isospin, hypercharge and isospin in z dir
  bool projected;
  bool smearedP;                      // smeared or not
  int t_slice;                        // t_slice
  bool creation_op;                   // creation op or not

};

  //**********************************************************************************************************************
//fucntions for reading
  void read_xml_ini( XMLReader& xml_in, int& npt, int& L, double& Xi, double& XiE, std::vector<NPtCorr_t>& had, int& num_matrix_elem, string& matrix_type,  Array1dO<string>& matrix_name, 
                string& print_zero, string& compute_phase, string& elab_dir,bool& make_redstar_xml, string& redstar_xml);
                
  //**********************************************************************************************************************
//fucntions for writing
  void write_ei( XMLWriter& xml, const std::string& path, const Eigen::Vector3d& input);
  void write_irrep(XMLWriter&  xml_out, IrrepLam_t& irrep_lam);
  void write_xml_out( XMLWriter& xml_out, int& npt, int& L, double& Xi, double& XiE, std::vector<NPtCorr_t>& had, int& num_matrix_elem, string& matrix_type,  Array1dO<string>& matrix_name, 
                string& print_zero, string& compute_phase, string& elab_dir, bool& make_redstar_xml, vector<NPtIrrepLam_t>& irreps, int& count );  

#endif