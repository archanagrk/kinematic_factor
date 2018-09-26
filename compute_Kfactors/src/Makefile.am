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

bin_PROGRAMS =  compute_kfactors compute_kfactors_with_phases

check_PROGRAMS = tes_ph

##############################################################

compute_kfactors_SOURCES = exe/compute_kfactors.cc  irreps_and_rows.cc  kfac_utils.cc  little_group.cc  subduction.cc  pol_vec.cc  rot_matrx.cc  k_factor.cc  levi_civita.cc 

compute_kfactors_with_phases_SOURCES = exe/compute_kfactors_with_phases.cc  irreps_and_rows.cc  kfac_utils.cc  little_group.cc  subduction.cc  pol_vec.cc  rot_matrx.cc  k_factor.cc  levi_civita.cc  phase.cc

tes_ph_SOURCES = test/tes_ph.cc  rot_matrx.cc little_group.cc kfac_utils.cc subduction.cc    k_factor.cc  levi_civita.cc  pol_vec.cc    irreps_and_rows.cc 

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
