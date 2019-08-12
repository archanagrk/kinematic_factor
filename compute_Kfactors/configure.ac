AC_INIT(Kfactors,1.0,archanar@jlab.org)
AC_CONFIG_AUX_DIR(config)

AM_INIT_AUTOMAKE([subdir-objects])


AC_LANG_CPLUSPLUS
AC_PROG_CXX(g++ cxx CC)
AC_PROG_CC(gcc cc)
AC_PROG_RANLIB

dnl Initialise automake
AM_INIT_AUTOMAKE
dnl
dnl
dnl --with-xxxx and --enable-xxxx switches 
dnl

dnl --with-adat=DIR
AC_ARG_WITH(adat,
  AC_HELP_STRING(
    [--with-adat=DIR],
    [Building on top of ADAT where ADAT is installed in DIR]
  ),
  [ADAT_HOME="$with_adat"]
)


if test "X${ADAT_HOME}X" = "XX" ; then
  AC_PATH_PROG(ADAT_CONFIG, [adat-config], [])
else
  AC_PATH_PROG(ADAT_CONFIG, [adat-config], [], [${ADAT_HOME}/bin:${PATH}])
fi

if test "X${ADAT_CONFIG}X" = "XX" ; then
  AC_MSG_ERROR([ADAT configuration program adat-config not found.])
fi

AC_MSG_NOTICE([Found ADAT configuration program ${ADAT_CONFIG}])
AC_SUBST(CXX, "`${ADAT_CONFIG} --cxx`")
AC_MSG_NOTICE([ADAT compiler: ${CXX}])
AC_SUBST(ADAT_CXXFLAGS, "`${ADAT_CONFIG} --cxxflags`")
AC_MSG_NOTICE([ADAT compile flags: ${ADAT_CXXFLAGS}])
AC_SUBST(ADAT_LDFLAGS,  "`${ADAT_CONFIG} --ldflags`")
AC_MSG_NOTICE([ADAT linking flags: ${ADAT_LDFLAGS}])
AC_SUBST(ADAT_LIBS,     "`${ADAT_CONFIG} --libs`")
AC_MSG_NOTICE([ADAT libraries flags: ${ADAT_LIBS}])

dnl include the m4 macro PAC_ADAT_LINK_CXX_FUNC
m4_include([m4/pac_adat_link_cxx_func.m4])
dnl This function is defined in m4/pac_adat_link_cxx_func.m4
dnl it tries to try to compile the program 
AC_MSG_CHECKING([if we can compile/link a simple adat program])
PAC_ADAT_LINK_CXX_FUNC(
  ${ADAT_CXXFLAGS},
  ${ADAT_LDFLAGS},
  ${ADAT_LIBS},
  ,
  ,
  [adat_link_ok=yes],
  [adat_link_ok=no]
)

if test "X${adat_link_ok}X" = "XyesX" ; then 
  AC_MSG_RESULT(yes)
else 
  AC_MSG_RESULT(no)
  AC_MSG_ERROR([Cannot compile/link a program with adat.
    Use --with-adat=<dir> to select a working version.])
fi 

##############################
# SEMBLE INTERFACE
##############################

dnl --with-semble=DIR
AC_ARG_WITH(semble,
  AC_HELP_STRING(
    [--with-semble=DIR],
    [Building on top of semble where semble is installed in DIR]
  ),
  [SEMBLE_HOME="$with_semble"]
)

if test "X${SEMBLE_HOME}X" = "XX" ; then
  AC_PATH_PROG(SEMBLE_CONFIG, [semble-config], [])
else
  AC_PATH_PROG(SEMBLE_CONFIG, [semble-config], [], [${SEMBLE_HOME}/bin:${PATH}])
fi

if test "X${SEMBLE_CONFIG}X" = "XX" ; then
  AC_MSG_ERROR([semble configuration program semble-config not found.
  		Cannot compile/linke without semble.
		Use --with-semble=<dir> to select a working version.])
fi

AC_MSG_NOTICE([Found semble configuration program ${SEMBLE_CONFIG}])
AC_SUBST(CXX, "`${SEMBLE_CONFIG} --cxx`")
AC_MSG_NOTICE([semble compiler: ${CXX}])
AC_SUBST(SEMBLE_CXXFLAGS, "`${SEMBLE_CONFIG} --cxxflags`")
AC_MSG_NOTICE([semble compile flags: ${SEMBLE_CXXFLAGS}])
AC_SUBST(SEMBLE_LDFLAGS,  "`${SEMBLE_CONFIG} --ldflags`")
AC_MSG_NOTICE([semble linking flags: ${SEMBLE_LDFLAGS}])
AC_SUBST(SEMBLE_LIBS,     "`${SEMBLE_CONFIG} --libs`")
AC_MSG_NOTICE([semble libraries flags: ${SEMBLE_LIBS}])




AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(src/Makefile)

AC_OUTPUT
