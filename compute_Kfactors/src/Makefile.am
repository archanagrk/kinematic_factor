#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/src \
		@CXXFLAGS@ @ADAT_CXXFLAGS@ 

AM_LDFLAGS  = -L@top_builddir@/src \
               @LDFLAGS@ @ADAT_LDFLAGS@

LDADD = @ADAT_LIBS@ @LIBS@ 


#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.

bin_PROGRAMS =  compute_kfactors compute_kfactors_with_phases compute_kfactors_j0 compute_kfactors_with_phases_j0 hadron_3pt_xml redstar_intrfc 

check_PROGRAMS = tes_ph extract_ph

##############################################################

compute_kfactors_SOURCES = exe/compute_kfactors.cc  lib/irreps_and_rows.cc  lib/kfac_utils.cc  lib/little_group.cc  lib/subduction.cc  lib/pol_vec.cc  lib/rot_matrx.cc  lib/k_factor.cc  lib/levi_civita.cc  io/xml_tools.cc 

compute_kfactors_with_phases_SOURCES = exe/compute_kfactors_with_phases.cc  lib/irreps_and_rows.cc  lib/kfac_utils.cc  lib/little_group.cc  lib/subduction.cc  lib/pol_vec.cc  lib/rot_matrx.cc  lib/k_factor.cc  lib/levi_civita.cc  lib/phase.cc

compute_kfactors_with_phases_j0_SOURCES = exe/compute_kfactors_with_phases_j0.cc  lib/irreps_and_rows.cc  lib/kfac_utils.cc  lib/little_group.cc  lib/subduction.cc  lib/pol_vec.cc  lib/rot_matrx.cc  lib/k_factor.cc  lib/levi_civita.cc  lib/phase.cc

compute_kfactors_j0_SOURCES = exe/compute_kfactors_j0.cc  lib/irreps_and_rows.cc  lib/kfac_utils.cc  lib/little_group.cc  lib/subduction.cc  lib/pol_vec.cc  lib/rot_matrx.cc  lib/k_factor.cc  lib/levi_civita.cc

tes_ph_SOURCES = test/tes_ph.cc  lib/rot_matrx.cc lib/little_group.cc lib/kfac_utils.cc lib/subduction.cc    lib/k_factor.cc  lib/levi_civita.cc  lib/pol_vec.cc    lib/irreps_and_rows.cc

extract_ph_SOURCES = test/extract_ph.cc lib/irreps_and_rows.cc  lib/kfac_utils.cc  lib/little_group.cc  lib/subduction.cc  lib/pol_vec.cc  lib/rot_matrx.cc  lib/k_factor.cc  lib/levi_civita.cc  lib/phase.cc

hadron_3pt_xml_SOURCES = io/hadron_3pt.xml.cc io/xml_tools.cc

redstar_intrfc_SOURCES = io/redstar_intrfc.cc  lib/irreps_and_rows.cc  lib/kfac_utils.cc  lib/little_group.cc  lib/subduction.cc lib/pol_vec.cc  lib/rot_matrx.cc  lib/k_factor.cc  lib/levi_civita.cc  io/xml_tools.cc

##############################################################

#BUILT_SOURCES = levi_civita.cc irreps_and_rows.cc kfac_utils.cc euler_mat.cc irreps_and_rows.cc

##############################################################

#
# The following rule will always try to rebuild libadat.a when you 
# try to compile anything else
#
####DEPENDENCIES = build_lib ${lib_LIBRARIES}
DEPENDENCIES = ${lib_LIBRARIES}

${bin_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

# build lib is a target that goes to the build dir of the library and 
# does a make to make sure all those dependencies are OK. In order
# for it to be done every time, we have to make it a 'phony' target
.PHONY: build_lib
build_lib:
	cd @top_builddir@/lib ; $(MAKE)
