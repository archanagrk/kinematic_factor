#
# Boiler Plate stuff: Autoconf Substitutions etc
#
AM_CXXFLAGS = -I@top_srcdir@/src \
		@CXXFLAGS@ @ADAT_CXXFLAGS@ @SEMBLE_CXXFLAGS@ 

AM_LDFLAGS  = -L@top_builddir@/src \
               @LDFLAGS@ @ADAT_LDFLAGS@ @SEMBLE_LDFLAGS@

LDADD = @ADAT_LIBS@ @LIBS@ @SEMBLE_LIBS@


#
# The programs to build
# 
## NOTE: if you want the programs to be build by default when 'make' is run
##   add them to bin_PROGRAMS.  But, if you don't want them to be installed
##   in @prefix@/bin by 'make install' then add them to nodist_PROGRAMS.
##   Another option is to add it to check_PROGRAMS, which means it will only
##   be built by 'make check' and is not installed.


bin_PROGRAMS = compute_matrix_prefactor gen_redstar_xml


##############################################################

compute_matrix_prefactor_SOURCES = exe/compute_matrix_prefactor.cc lib/irreps_and_rows.cc  lib/kfac_utils.cc  lib/little_group.cc  lib/subduction.cc  lib/subd_pol_vec.cc  lib/pol_vec.cc  lib/rot_matrx.cc  lib/kfactor_pigammarho.cc lib/kfac_params.cc lib/k_factor_factory.cc  lib/levi_civita.cc  io/xml_tools.cc  lib/phase.cc

gen_redstar_xml_SOURCES 		 = io/gen_redstar_xml.cc  lib/kfac_utils.cc  io/xml_tools.cc

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
