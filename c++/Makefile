# $Id$ -*-makefile-*-

# Purpose: Makefile for C++ module c++

# Copyright (C) 1994--2017 Charlie Zender
# License: GNU General Public License (GPL) Version 3
# See http://www.gnu.org/copyleft/gpl.html for full license text

# Quickie test copies:
# scp ~/sw/c++/Makefile esmf.ess.uci.edu:sw/c++
# scp ~/sw/c++/Makefile goldhill.cgd.ucar.edu:sw/c++
# scp ~/sw/c++/Makefile soot.ess.uci.edu:sw/c++

# Usage (testing):
# NB: CPP token controls OpenMP activation so ccache may use non-OpenMP objects
# make CCACHE=N # Ensure OpenMP is appropriately (de-)activated
# make -W foo.c
# make --no-builtin-variables     <-- Turn off implicit, built-in variables
# make --warn-undefined-variables <-- Fail if Makefile attempts to use undefined variables

# Usage (Compilation):
# cd ~/sw/c++;make;cd -                  # Default build
# cd ~/sw/c++;make OPTS=D NETCDF4=Y;cd - # Build with netCDF4
# cd ~/sw/c++;make OPTS=D NETCDF4=Y SZ=Y;cd - # Build with netCDF4 and Szip
# 32-bit ABI on NCAR AIX systems:
# NETCDF_LIB=/usr/local/lib32/r4i4 make --jobs=1 ABI=32 OPTS=D;cd -
# 64-bit ABI netCDF3 on NCAR AIX systems (bluefire):
# cd ~/sw/c++;GSL_INC='/contrib/gsl-1.12/include' GSL_LIB='/contrib/gsl-1.12/lib -lgsl' NETCDF_LIB='/usr/local/lib64/r4i4'  make --jobs=1 OPTS=D ABI=64 all;cd -
# 64-bit ABI netCDF4 on NCAR AIX systems (bluefire):
# cd ~/sw/c++;CURL_LIB='/contrib/curl/7.21.2/lib' GSL_INC='/contrib/gsl-1.12/include' GSL_LIB='/contrib/gsl-1.12/lib -lgsl' HDF5_ROOT='/contrib/hdf5-1.8.7_seq' NETCDF_ROOT='/contrib/netcdf/4.1.3_seq' SZ_LIB='/contrib/szip/lib' make --jobs=1 ABI=64 NETCDF4=Y OPTS=D SZ=Y all;cd -

# Makefile conventions for optimization levels (OPTS levels)
# OPTS=O: Optimize: Produce symbol tables, use architecture-specific speed-ups
# OPTS=R: Regular: No custom switches
# OPTS=D: Debugging: Produce symbol tables, no bounds-check or profiling 
# OPTS=X: EXtreme debugging: Produce symbol tables, check bounds, enable profiling

# In rule 
# foobar.o : foo.c bar.c
# $@ is target = foobar.o
# $< is first dependency = foo.c
# $^ is name of all dependencies = foo.c bar.c
# $* is stem of match, and so much more

# Recursively expanded variables are defined with `=':
# Values are installed verbatim and expanded when substitution occurs
# Variable may be referenced before it is defined and expansion will be correct
# However variable may not be modified (e.g., appended) after definition

# Simply expanded variables are defined with `:=':
# Values are immediately expanded at the time variable is defined
# Variables must be defined before being used
# Variables may be modified without causing infinite recursion

# Top-level tokens defining directory structure
# These tokens may be over-ridden by environment variables or when invoking make, e.g. make DODS=Y
ifndef MY_BIN_DIR
 MY_BIN_DIR := ../bin
endif
ifndef MY_LIB_DIR
 MY_LIB_DIR := ../lib
endif
ifndef MY_OBJ_DIR
 MY_OBJ_DIR := ../obj
endif
MY_BLD_DIR := .
MY_DAT_DIR := .
MY_DOC_DIR := .
MY_DPN_DIR := ${MY_OBJ_DIR}
MY_INC_DIR := ${HOME}/include

# Primary tokens which determine build options
# Specify non-default when invoking make, e.g. make DODS=Y
ifndef ${ABI}
# 32- vs. 64-bit ABI: 32=32-bit mode, 64=64-bit mode (default) if available
 ABI := 64
endif # endif ABI
ifndef ${CCACHE}
 CCACHE := N
endif # CCACHE
ifndef ${DBG}
# Debugging token N=No (default) Y=Yes
 DBG := N
endif # endif DBG
ifndef GCC_RCH_ARG
 GCC_RCH_ARG :=
endif # endif GCC_RCH_ARG
ifdef GSL_INC
 GSL_INC := -I${GSL_INC}
else
 GSL_INC := $(shell gsl-config --cflags)
endif # endif GSL_INC
ifdef GSL_LIB
 GSL_LIB := -L${GSL_LIB}
else
 GSL_LIB := $(shell gsl-config --libs)
endif # endif GSL_LIB
ifndef I18N
 I18N := Y
 MY_SHR_DIR := ${HOME}/share
 MY_ES_DIR := ${MY_SHR_DIR}/locale/es/LC_MESSAGES
 MY_FR_DIR := ${MY_SHR_DIR}/locale/fr/LC_MESSAGES
endif # I18N
ifndef ICC_RCH_ARG
 ICC_RCH_ARG :=
endif # ICC_RCH_ARG
ifndef MK_DPN
 MK_DPN = ${CPP} -M # NB: Recursive expansion required
 MK_DPN_CXX = ${CXX} -M # NB: Recursive expansion required
endif # MK_DPN
ifndef ${LAMMPI}
 LAMMPI := Y
endif # LAMMPI
ifndef ${MPI}
 MPI := N
endif # MPI
ifndef ${NETCDF4} # netCDF4 support
 NETCDF4 := N
endif # endif NETCDF4
ifndef NETCDF_ROOT
 NETCDF_ROOT := /usr/local
endif
ifndef NETCDF4_ROOT
 NETCDF4_ROOT := ${NETCDF_ROOT}
endif
ifndef NETCDF_INC
 NETCDF_INC := ${NETCDF_ROOT}/include # Directory containing netcdf.h
endif
ifndef NETCDF_LIB
 NETCDF_LIB := ${NETCDF_ROOT}/lib # Directory containing libnetcdf.a
endif
ifndef ${OMP} # OpenMP
 OMP := Y
endif # endif OMP
ifndef OPTS
 OPTS := O
endif # endif OPTS
ifndef PGI_RCH_ARG
 PGI_RCH_ARG :=
endif # endif PGI_RCH_ARG
ifndef ${PRC}
# Precision token S=Single, D=Double (default), L=Long Double
 PRC := D
endif # endif PRC
ifndef PSC_RCH_ARG
 PSC_RCH_ARG :=
endif # endif PSC_RCH_ARG
ifndef PVM_ARCH
 PVM_ARCH := $(shell ${HOME}/sh/pvmgetarch)
endif
ifndef STC
# Created statically linked executable
 STC := N
endif
ifndef SZ
# Link to Szip library
 SZ := N
endif
ifndef UNAMES
 UNAMES := $(shell uname -s)
endif
ifndef USR_TKN
 USR_TKN := -DABORT_ON_ERROR
endif # endif USR_TKN
ifndef VRS_SNG
 VRS_SNG := $(shell date +%Y%m%d)
endif # endif VRS_SNG

# Derived-tokens based on primary tokens
# Internationalize NCO with i18n features
ifeq (${I18N},Y)
 MY_SHR_DIR := ${HOME}/share
 MY_ES_DIR := ${MY_SHR_DIR}/locale/es/LC_MESSAGES
 MY_FR_DIR := ${MY_SHR_DIR}/locale/fr/LC_MESSAGES
endif

ifeq (${NETCDF4},Y)
 CPPFLAGS += -DENABLE_NETCDF4
endif
# endif NETCDF4

ifeq (${PRC},D)
 USR_TKN := -DPRC_DBL ${USR_TKN}
else
ifeq (${PRC},L)
 USR_TKN := -DPRC_LDB ${USR_TKN}
else
 USR_TKN := -DPRC_FLT ${USR_TKN}
endif # endif PRC
endif # endif PRC

# NB: CPP_TKN is recursively expanded variable, define only when components are ready
CPP_TKN = ${USR_TKN} -D${PVM_ARCH} -DNO_NETCDF_2 -DVERSION='${VRS_SNG}' -DHOSTNAME='${HOST}' -DUSER='${USER}'

# Default build. DAP, NETCDF4, and PNETCDF clauses may overwrite this
NC_LDFLAGS := -L${NETCDF_LIB}
NC_LIBS := -lnetcdf 

ifeq (${NETCDF4},Y)
# Enable netCDF4 functionality
ifndef HDF5_ROOT
 HDF5_ROOT := ${NETCDF4_ROOT}
endif
ifndef HDF5_INC
 HDF5_INC := ${HDF5_ROOT}/include # Directory containing hdf5.h
endif
ifndef HDF5_LIB
 HDF5_LIB := ${HDF5_ROOT}/lib # Directory containing libhdf5.a
endif
 NETCDF_INC := ${NETCDF4_ROOT}/include
 NETCDF_LIB := ${NETCDF4_ROOT}/lib
 NC_LDFLAGS := -L${NETCDF_LIB} -L${HDF5_LIB}
 NC_LIBS := -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl
ifdef CURL_LIB
 NC_LDFLAGS += -L${CURL_LIB}
endif # end 
ifeq (${SZ},Y)
ifndef SZ_LIB
 SZ_LIB := /sopt/lib # Directory containing libsz.a
endif
 NC_LDFLAGS += -L${SZ_LIB}
 NC_LIBS += -lsz
endif # !SZ
ifneq (${HDF5_INC},${NETCDF_INC})
 NETCDF_INC += -I${HDF5_INC}
endif # end 
endif # endif NETCDF4

# These tokens should not be altered by hand
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
# Decide among the plethora of Linux compilers
ifndef LINUX_CXX
# C++ compiler for Linux
 LINUX_CXX := g++
 #LINUX_CXX := clang
 #LINUX_CXX := como
 #LINUX_CXX := icpc
 #LINUX_CXX := insure
 #LINUX_CXX := pathCC
 #LINUX_CXX := pgCC
endif # endif LINUX_CXX
ifndef LINUX_CC
# C compiler for Linux
 LINUX_CC := gcc -std=c99 -pedantic -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := clang
 #LINUX_CC := como --c99
 #LINUX_CC := icc -std=c99 -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := insure
 #LINUX_CC := pathcc -std=c99 -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := pgcc -c9x -D_BSD_SOURCE -D_POSIX_SOURCE
endif # endif LINUX_CC
ifndef LINUX_FC
# Fortran compiler for Linux
 LINUX_FC := g95
 #LINUX_FC := gfortran
 #LINUX_FC := lf95
 #LINUX_FC := ifort
 #LINUX_FC := pathf95
 #LINUX_FC := pgf90
endif # endif LINUX_CC
endif # endif LINUX

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := ${null} ${null}
comma := ,
newline := \n
# '/' and '+' appear in filenames ('/' is directory separator)
# Operating on these with Perl is problematic since they are special Rx characters
# We replace `/' and '+' by non-special Rx's, call perl, then demangle
# Unique character(s) to substitute for '/' and '+' before passing to perl Rx
slash_rx := cszzsc
plus_rx := xdikmj
# Unique character(s) to replace by ${slash_rx} after processing by Perl
slash := /
plus := +
MY_OBJ_DIR_RX := $(subst ${slash},${slash_rx},${MY_OBJ_DIR})
MY_DPN_DIR_RX := $(subst ${slash},${slash_rx},${MY_DPN_DIR})

# Source file names with directories removed
# List source targets alphabetically by "category":
MDL_SRC := ccc.cc gsl.cc i18n.cc tst.cc cls.cc snw.cc slr_zen_ngl_cos # Non-library
MDL_SRC += bit.cc cln.cc mmr.cc mth.cc omp.cc spc_fnc.cc utl.cc vec.cc xtr.cc # libcsz_c++
MDL_SRC += blm.cc chm.cc fio.cc flx_sfc.cc htrn_c++.cc idx_rfr.cc idx_rfr_H2O.cc lbl.cc mie_sln.cc mie_MaS99.cc mie_Sor01.cc mie_Wis79.cc tdy.cc aer.cc flp.cc mnr_dst.cc nco.cc phz_fnc.cc psd.cc spc_bbd.cc spc_slr.cc ssl.cc sz_grd.cc wvl_grd.cc # libcsm_c++
# Directories to search for source files
MDL_PTH := . ${HOME}/include
# Dependency list for executable
MDL_OBJ := $(addprefix ${MY_OBJ_DIR}/,$(addsuffix .o, $(basename ${MDL_SRC}))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix ${MY_DPN_DIR}/,$(addsuffix .d, $(basename ${MDL_SRC}))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.d file
VPATH := $(subst ${space},:,${MDL_PTH})
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,${MDL_PTH},-I${dir})

# Variables having to do with binary executables created by module
MDL_BIN_TRG := ccc gsl i18n tst snw slr_zen_ngl_cos # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
# List header targets alphabetically by "category":
MDL_INC_TRG := getopt.h getopt_bsd.h dbg.hh lsm.hh phys_cst.hh # Raw (no functions)
MDL_INC_TRG += libcsz_c++.hh a2d.hh bit.hh cln.hh cln_cst.hh mmr.hh mth.hh mth_cst.hh omp.hh spc_fnc.hh utl.hh vec.hh xtr.hh # libcsz_c++
MDL_INC_TRG += libcsm_c++.hh blm.hh chm.hh fio.hh flx_sfc.hh htrn_c++.hh idx_rfr.hh idx_rfr_H2O.hh lbl.hh mie_sln.hh mie_MaS99.hh mie_Sor01.hh mie_Wis79.hh tdy.hh aer.hh flp.hh mnr_dst.hh nco.hh pdf.hh phz_fnc.hh psd.hh rt.hh spc_bbd.hh spc_slr.hh ssl.hh sz_grd.hh wvl_grd.hh # libcsm_c++
MDL_INC_SYM_LNK := # Symbolic links
MDL_INC_STB = ${MDL_INC_TRG} ${MDL_INC_SYM_LNK} # All header files in ${MY_INC_DIR}
MDL_INC = $(addprefix ${MY_INC_DIR}/,${MDL_INC_STB}) # dst_cln removes these files

# Redefine default C and C++ pattern rules
${MY_OBJ_DIR}/%.o : %.c
	${CC} ${CPPFLAGS} ${CFLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)
${MY_OBJ_DIR}/%.o : %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)

# Default Fortran pattern rules: CRAY and RS6K must override these rules
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<

# Rules for installing header files
${MY_INC_DIR}/%.h : %.h
	cp -f -p $(notdir $@) $@
${MY_INC_DIR}/%.hh : %.hh
	cp -f -p $(notdir $@) $@

# Rules for installing i18n files
%.po : %.cc
	xgettext --default-domain=$* --join-existing $<
${MY_ES_DIR}/%.mo : %.po
#	Linux version accepts more arguments than Solaris version
#	msgfmt --output-file=$@ --statistics $<
	msgfmt -o $@ $<

# Automatically generate dependency file for each source file
# $* is stem, e.g., f
# $@ is filename of target, e.g., f.d
# Linux gcc may return an extra `.F' on Fortran names, e.g., `hello.F.o: hello.F'
# (.F)? gets rid of this extra `.F'
${MY_DPN_DIR}/%.d : %.F
# Following command makes, e.g., f.d begin "f.o f.d : f.F ..."
# Since f.o is not preceded by ${MY_OBJ_DIR}, objects are not recompiled when sources are touched.
#	${MK_DPN} ${CPPFLAGS} $< | perl -p -e 's/$*\.F\.o/$*.o $@/g;' > $@
# Following command makes, e.g., f.d begin "/home/zender/obj/LINUX/f.o f.d : f.F ..."
# This works fairly well, but is a hack
# First pattern substitutes MY_OBJ_DIR_RX, which has placeholders for slashes
# Second pattern substitutes slashes for the placeholders
	${MK_DPN} ${CPPFLAGS} $< | perl -p -e 's/$*(\.F)?\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@
# Following command makes, e.g., f.d begin "${MY_OBJ_DIR}/f.o f.d : f.F ..."
# This would be ideal but I cannot get dollar sign to cooperate
#	${MK_DPN} ${CPPFLAGS} $< | perl -p -e 's/$*\.F\.o/\${dollar}MY_OBJ_DIR\/$*.o $@/g;' > $@

${MY_DPN_DIR}/%.d : %.f
	${MK_DPN} ${CPPFLAGS} $< | perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.c
#	${MK_DPN} ${CPPFLAGS} $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	${MK_DPN} ${CPPFLAGS} $< | perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.cc
# NB: Use ${CXX} rather than ${CPP} on C++ files for now because, e.g., SUNMP cpp does not pre-process .cc files quite correctly
# Extra hack to allow C++ filenames to contain '+' character
# $(subst ${plus},${plus_rx},${*}) is filename stub with an Rx in place of '+'
	${MK_DPN_CXX} ${CPPFLAGS} $< | perl -p -e 's/\${plus}/${plus_rx}/g;s/$(subst ${plus},${plus_rx},${*})\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g;s/${plus_rx}/\${plus}/g' > $@

# First LDFLAGS is for typical C programs with netCDF, math, and networking
# Second LDFLAGS, when present, enables C/Fortran linking

# Manually define autotools tokens normally defined in HAVE_CONFIG_H
# Initialize OS-specific tokens to empty
CPP_TKN_OS := 
CPP_TKN_OS += -DHAVE_C99
ifneq (${null},$(findstring ${PVM_ARCH},LINUXALPHALINUXAMD64LINUXARMFREEBSDWIN32))
# HAVE_INLINE is GNU standard to activate inline functions (used by GSL)
 CPP_TKN_OS += -DHAVE_INLINE
endif # !LINUX
# Message Passing Interface
ifeq (${MPI},Y)
 CPP_TKN_OS += -DENABLE_MPI
endif # !MPI
ifneq (${PVM_ARCH},ALPHA)
# HAVE_LONG_LONG
# long long is ISO C99 standard but is neither ISO C++ nor ANSI C standard
# All compilers except Tru64 support it.
 CPP_TKN_OS += -DHAVE_LONG_LONG
endif # ALPHA
ifneq (${null},$(findstring ${PVM_ARCH},LINUXALPHALINUXAMD64LINUXARMFREEBSDWIN32))
# fxm: GCC numeric_limits do not yet work
# HAVE_NUMERIC_LIMITS
# numeric_limits<T> templates can replace all of <cfloat> and <climits>
# CPP_TKN_OS += -DHAVE_NUMERIC_LIMITS
endif # LINUX
# fxm: Define HAVE_LIBINTL, HAVE_LOCALE_H, HAVE_GETTEXT, HAVE_OMP_H

# Works on AIX and AIX46K
ifneq (${null},$(findstring AIX,${PVM_ARCH}))
# 20030804: Always use re-entrant (_r) compilers---Jim Edwards NCAR/IBM 
# -qlanglvl=extended: Include extensions to standard
# -qlanglvl=extc99: Include orthogonal extensions to C99 standard
CC := xlc_r -qlanglvl=extc99
#CC := gcc -std=c99 -pedantic -D_BSD_SOURCE
CXX := xlC_r
# CXX := g++
# fxm: 20040311 should this be cpp? gcc -M?
CPP := xlc -c -qlanglvl=extc99 -qsuppress=1501-245
ifneq (${null},$(findstring xl,${CC}))
# Visual Age compiler headers must occur before g++ headers
 CPP_TKN_OS += -I/usr/vacpp/include
endif # xlC compilers
# Add /usr/local/include for libintl.h explicitly until netCDF is moved there
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC} -I/usr/local/include
FC := xlf95_r
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ ${NC_LDFLAGS} ${NC_LIBS}
LEX := flex
LINT := lint
YACC := bison
# AIX VA Compiler Collection
ifneq (${null},$(findstring xl,${CC}))
# Visual Age compiler defaults specified in /etc/vac.cfg
# Additional switch to fix compiler warnings on csz.c
# -qarch=auto : Automatically detect architecture of compiling machine and assume execution on same machine
# -qlonglong allow long long integers (and strtoll(), strtoull()) (default on in C not C++). Redundant with -qlanglvl=extc99
# -qmaxmem=num Limit memory used by space intensive optimizations to <num> kilobytes
# -qspill=size Size in B of register allocation spill area, mie needs > 1040 B
# -qsrcmsg prints transgressing source line with finger
# -qsuppress=1501-245 : Suppress RLIM_INFINITY memory message due to ulimits
# -qtune=auto : Optimize executable for architecture detected during compilation
 CFLAGS := -qlonglong -qmaxmem=8192 -qspill=2048 -qsrcmsg -qsuppress=1501-245
 FFLAGS := -NS2000 -qfixed=132 -qsrcmsg
# -bh:5 suppresses annoying messages from xlC linker WARNING: Duplicate symbol: ...
 LDFLAGS += -bh:5 -qsuppress=1501-245
#LDFLAGS += -lxlf90 # Required for linking Fortran objects
ifeq (${OMP},Y)
# -qsmp=noauto : Turn on SMP/OMP code generation but do no automatic parallelization
# -qsmp=omp : Use industry standard OMP without IBM extensions
 OMP_FLG_C := -qsmp=omp
 OMP_FLG_F := -qsmp=omp
else
 CPP_DFN += -U_OPENMP
 FC := xlf95
endif # endif OMP
ifeq (${OPTS},O)
# -O : -O3 is safe, -O5 is dangerous
# -qstrict: Ensure that -O3 optimizations do not alter program semantics
# -Q : Inline all appropriate subprograms
 CFLAGS += -O3 -g -qstrict -Q
 FFLAGS += -O3 -g -qstrict -Q
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g
endif
ifeq (${OPTS},X)
# -qcheck=all generates run-time bounds, divide-by-zero, and null-dereference checks
# -qflttrap generates instructions for floating point exceptions
# -qidirfirst uses headers found in -I directories first
# -qlonglong allow long long integers
# -qmakedep creates .d file
# -qwarn64 check for possible long-to-integer or pointer-to-integer truncation
# -qhalt=e stop compiler if error severity equals or exceeds i, w, e, s, u
 CFLAGS += -g -qflttrap -qidirfirst -qlonglong -qwarn64 -qcheck=all -qhalt=s
 FFLAGS += -g
endif
ifeq (${ABI},64)
 AR := ar -X 64
 CFLAGS += -q64
 FFLAGS += -q64
 LDFLAGS += -q64 
else
 CPPFLAGS += -D_LARGE_FILES
endif # endif ABI
# Additional flags for AIX:
# -M Generate information to be included in a "make" description file; output goes to .u file
# -c Do not send object files to the linkage editor
# -P Preprocess but do not compile; output goes to .i file
# Using -P causes additional warning messages about lm 
# Not using -P causes *.o files to be created twice
${MY_DPN_DIR}/%.d : %.c
	${MK_DPN} ${CPPFLAGS} $< ;perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' $*.u > $@ ; \
	rm -f $*.i $*.o $*.u;
${MY_DPN_DIR}/%.d : %.cc
	${MK_DPN} ${CPPFLAGS} $< ;perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' $*.u > $@ ; \
	rm -f $*.i $*.o $*.u;
endif # endif AIX VA Compiler Collection
# GNU Compiler Collection
ifneq (${null},$(findstring gcc,${CC}))
 CFLAGS := -Wall -Wunused
ifeq (${OPTS},O)
 CFLAGS += -O -g ${GCC_RCH_ARG}
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS += 
endif
ifeq (${OPTS},X)
 CFLAGS += -g -O
 LDFLAGS += /usr/local/lib/ccmalloc-g++.o -L/usr/local/lib -lccmalloc -ldl
endif
ifeq (${ABI},64)
 CC += -maix64 
 CXX += -maix64 
endif # endif ABI
 CXXFLAGS := ${CFLAGS}
endif # endif GNU Compiler Collection
# -q64: Select 64-bit compiler mode (required for accessing large files)
# -qwarn64: Warn on possible long-to-integer or pointer-to-integer truncation
CXXFLAGS := ${CFLAGS}
ifeq (${OMP},Y)
 CFLAGS += ${OMP_FLG_C}
 CXXFLAGS += ${OMP_FLG_C}
 FFLAGS += ${OMP_FLG_F}
# LD behavior assumes C source code
 LDFLAGS := ${OMP_FLG_C} ${LDFLAGS}
endif # endif OMP
endif
# endif AIX

ifeq (${PVM_ARCH},ALPHA)
ifeq (${OMP},Y)
 OMP_FLG := -omp
endif # endif OMP
CXX := cxx -std strict_ansi
#CXX := cxx -std ansi -D__NO_USE_STD_IOSTREAM
CC := cc
CFLAGS := ${OMP_FLG}
# OSF1 is for CCM compatibility
CPP_TKN_OS += -DOSF1
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := f90
FFLAGS := -c ${OMP_FLG} -automatic -fpe3
FIXEDFLAGS := -extend_source
FREEFLAGS := -DHIDE_SHR_MSG -free
LD := ld
LDFLAGS += ${OMP_FLG} -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ ${NC_LDFLAGS} ${NC_LIBS}
LEX := flex
LINT := lint
YACC := bison
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
endif
ifeq (${OPTS},O)
 CFLAGS += -O2 -g -ansi_alias
 FFLAGS += -O3 -g -inline speed
endif
ifeq (${OPTS},D)
 CFLAGS += -g -check_bounds -check -check_omp
 FFLAGS += -g3 -C
endif
ifeq (${OPTS},X)
 CFLAGS := -g -N 132
 FFLAGS := -g -check bounds -check omp_bindings -check overflow -check underflow
endif
CXXFLAGS := ${CFLAGS}
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif ALPHA

ifeq (${PVM_ARCH},CRAY)
CXX := g++
CC := cc
CPP := cpp
CPP_TKN_OS +=
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := f90
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ ${NC_LDFLAGS} ${NC_LIBS}
#LDFLAGS += -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O -g
 FFLAGS := -N 132
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS := -g -N 132
endif
ifeq (${OPTS},X)
 CFLAGS += -g -N 132
 FFLAGS := -g -N 132 -e i
endif
# 1997/10/21 Added -P to suppress #line # directives on Fortran files
${MY_OBJ_DIR}/%.o : %.F
	${CPP} -P ${CPPFLAGS} $< > $(patsubst %.F,%.f,$(notdir $<))
	${FC} -c ${FFLAGS} $(patsubst %.F,%.f,$(notdir $<)) 
	-mv -f $(notdir $@) ${MY_OBJ_DIR}
	rm -f $(patsubst %.F,%.f,$(notdir $<)) 
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} $<
	mv -f $(notdir $@) ${MY_OBJ_DIR}
endif
# endif CRAY

# Works on LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, and FREEBSD and MACOSX
ifneq (${null},$(findstring ${PVM_ARCH},LINUXALPHALINUXAMD64LINUXARMFREEBSDMACOSX))
CXX := ${LINUX_CXX}
CC := ${LINUX_CC}
CPP := ${CXX}
CPP_TKN_OS +=
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := ${LINUX_FC}
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ ${NC_LDFLAGS} ${NC_LIBS}
ifeq (${PVM_ARCH},MACOSX)
 LDFLAGS += -lresolv
endif # endif MACOSX
LEX := flex
LINT := lint
YACC := bison
# Comeau C Compiler
ifeq (${CXX},como)
# Diagnose and suppress warnings and errors
# 177: "declared but not referenced" warnings
# 550: "set but not used" warnings
 CFLAGS := --diag_suppress 177,550 --display_error_number
 CPPFLAGS := $(filter-out -DHAVE_LONG_LONG,${CPPFLAGS}) -I/usr/lib/gcc-lib/i486-linux/3.3.2/include
 LDFLAGS := ${COMOROOT}/libcomo/libcomo.a ${LDFLAGS}
ifeq (${OPTS},O)
 CFLAGS += -O -g
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS +=
endif
ifeq (${OPTS},X)
 CFLAGS += -g
endif
 CXXFLAGS := ${CFLAGS}
endif # endif Comeau C Compiler
# GNU Compiler Collection or LLVM
# 20140204: gcc and clang should receive identical options
ifeq (g++,$(firstword ${CXX}))
 GCC_OR_CLANG := Y
endif
ifeq (clang,$(firstword ${CXX}))
 GCC_OR_CLANG := Y
endif
ifeq (${GCC_OR_CLANG},Y)
# Enable Large File Support (LFS) by default
 CPPFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
 CFLAGS := -Wall -Wunused
# Compilation flags for numerical routines recommended by GSL 1.3 manual, p. 397
# CFLAGS += -Werror -Wall -W -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common -g -O4
# Compilation flags recommended by Coyote Gulch: -march=pentium4 -mfpmath=sse -fomit-frame-pointer -ffast-math -O3
# Compilation flags recommended by GSL that I like and use:
# -D_BSD_SOURCE: Support 4.3 BSD Unix extensions to ANSI C (prevents nameser.h warnings)
# -D_POSIX_SOURCE: Support POSIX.1 standard additions to ANSI C (prevents fileno warnings)
# -pedantic: Disallow non-ISO constructs (including type long long) (sometimes useful)
# -Werror: Consider warnings as errors
# -W: Extra warnings, including missing return values, comparison of signed with unsigned
# -Wall: Warn about common programming problems
# -Wcast-align: Warn if casting pointer to type of different size
# -Wcast-qual: Warn if const qualifier removed from pointer
# -Werror: Consider warnings as errors
# -Wmissing-prototypes: Warn if missing prototypes
# -Wpointer-arith: Warn if pointer arithmetic on types without size, e.g., void
# -Wshadow: Warn if local variable has same name as other local variable
# -Wswitch: Warn if switch statement has enumerated index and case label outside enumeration range
# -Wunused: Warn on unused functions, labels, parameters, values, and variables
# -Wwrite-strings: Apply const-qualifier to string constants, die if overwritten
# -fno-common: Prevent global variables from being simultaneously defined in different files
# -g: Put debugging symbols in executable
# -pg: Enable profiling, generate gmon.out output files (also needed by linker)
# -O4: Turn on optimization so uninitialized variables are flagged. Downside: optimizes-out many debugging symbols
# Unused compilation flags that may help:
# -fexpensive-optimizations:
# -ffast-math: Implement non-IEEE shortcuts
# -fno-inline: Valgrind recommends this
# -fstrict-aliasing: Assume objects of different types never reside at same address. Makes union access rules important. Enabled by -O2
# -fomit-frame-pointer: Avoid extra frame instructions, reduce debuggability, enabled by -O*
# -funroll-loops: Unroll loops if iterations known at compile time
# -m64: Compile in 64-bit mode for x86_64
# -mfpmath=sse,387: 
# -mmmx:
# -msse:
# -msse2: 
# -m3dnow:
# -m3dnowext:
# -pipe:
# Compilation flags recommended by GSL that I do not like and do not use (give reason why):
# -ansi: Support only strict ANSI C. Equivalent to -std=c89, conflicts with -std=c99
# -Waggregate-return: Warn if functions return aggregates like structures or unions
# -Wconversion: Warn if converting signed to unsigned. Intended for obsolete, non-prototyped code. Triggers fabsf(), sqrtf(), warnings.
# -Wnested-externs: Warn if extern is encountered within function. C only?
# -Wstrict-prototypes: Warn if inconsistent prototypes. C only?
# -Wtraditional: Warn if constructs differ between traditional and ANSI C. C only?
# -Dinline=: inline is not an ANSI keyword, must undefine inline to work with -ansi
# -fshort-enums: Make enums as short as possible, ususally non-int. Do not ever invoke this! This breaks ABI and causes subtle problems
ifeq (${OMP},Y)
ifneq (clang,$(firstword ${CXX}))
# 20140526 clang does not recognize/utilize -fopenmp
 OMP_FLG_C := -fopenmp
 OMP_FLG_F := -fopenmp
 LDFLAGS += -lgomp -lpthread
 OMP_FLG_LD := -lgomp -lpthread
endif # endif clang
endif # endif OMP
ifeq (clang,$(firstword ${CXX}))
# 20140531 avoid clang error that -O4 is ame as -O3
 CFLAGS := $(subst -O4,-O3,${CFLAGS})
 LDFLAGS += -lstdc++
endif # endif clang
ifeq (${OPTS},O)
 CFLAGS += -O -g ${GCC_RCH_ARG}
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS += 
endif
ifeq (${OPTS},X)
 CFLAGS += -g -O -pg -fno-inline
 LDFLAGS += -pg
# LDFLAGS += /usr/local/lib/ccmalloc-g++.o -L/usr/local/lib -lccmalloc -ldl
endif
ifneq (${null},$(findstring AMD64,${PVM_ARCH}))
ifeq (${ABI},64)
 CFLAGS += -m64
 FFLAGS += -m64
 LDFLAGS += -m64
endif # endif ABI
endif # endif LINUXAMD64
 CXXFLAGS := ${CFLAGS}
endif # endif GNU Compiler Collection or LLVM
# Intel (Kai) C++ Compiler
ifeq (icpc,$(firstword ${CXX}))
# -cxxlib-gcc: link using C++ run-time libraries provided with GCC (default) (generated code is binary compatible with g++)
# -lcprts : Eliminates most undefined references when linking with -cxxlib-gcc
# Remaining unresolved symbols are
# -cxxlib-icc: link using C++ run-time libraries provided by Intel (generated code is not binary compatible with g++)
# -fast: enable -xP -O3 -ipo -static
# -ipo[n]: enable multi-file IP optimizations (between files)
# -no-gcc: do not define  __GNUC__, __GNUC_MINOR__, and __GNUC_PATCHLEVEL__ macros
# -static: prevents linking with shared libraries
# -std=c99: Enable C99 support for C programs
# -aX<codes>: generate generic IA32 code and code specialized for processors specified by <codes> where codes are one or more of following:
# B: Intel Pentium M and compatible Intel processors
# K: Intel Pentium III and compatible Intel processors
# N: Intel Pentium 4 and compatible Intel processors
# P: Intel Pentium 4 processors with SSE3 extensions
# W: Intel Pentium 4 and compatible Intel processors
# -x<codes>: same as -aX only no generic code, run exclusively on <codes>
# -Wall: enable all warnings
# -Werror:force warnings to be reported as errors
# -w0: display errors (same as -w)
# -w1: display warnings and errors (DEFAULT)
# -w2: display remarks, warnings, and errors
# -wd<L1>[,<L2>,...] disable diagnostics L1 through LN
# warning #274: declaration is not visible outside of function
# remark #981: operands are evaluated in unspecified order
# remark #810: conversion from "double" to "float" may lose significant bits
# remark #1572: floating-point equality and inequality comparisons are unreliable
 CFLAGS := -w1
 CPPFLAGS += -I${ICC_ROOT}/include -I${ICC_ROOT}/include/c++ -no-gcc
 LDFLAGS :=  ${LDFLAGS} -lsvml
ifeq (${OMP},Y)
 OMP_FLG_C := -openmp -openmp_report0
endif # endif OMP
ifeq (${OPTS},O)
 CFLAGS += -O3 -g ${ICC_RCH_ARG}
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS +=
endif
ifeq (${OPTS},X)
 CFLAGS += -g -Wall -wd810,981,1572 -inline_debug_info
endif
 CXXFLAGS := ${CFLAGS}
endif # endif Intel (Kai) C++ Compiler
# Intel (Kai) Fortran Compiler
ifeq (${FC},ifort)
# -132: specify 132 column lines for fixed form sources (same as -extend_source)
# -axW: generate generic IA32 code and specialized code exclusively for Pentium 4 processors (requires linking to libsvml.a). Note that -axW dramatically and unpleasantly increases verbosity of compilation output
# -diag-disable 5268: quiet warnings that comments exceed allowed length
# -e95: die on warnings for non-standard Fortran95 (which includes OpenMP)
# -fast: enable -xP -O3 -ipo -static
# -fpconstant: extend precision of single precision constants assigned to double precision variables to double precision
# -fpp2: necessary, but not sufficient, for OpenMP
# -ftz: enable/disable flush denormal results to zero
# -ip: enable single-file IP optimizations
# -ipo: enable multi-file IP optimizations
# -module: specifies path to place .mod files
# -mp: maintain floating point precision (disables some optimizations) (deprecated 20090129)
# -nus: do not append an underscore to external subroutine names
# -vec_report1: indicate vectorized loops (DEFAULT)
# -tune pn4: optimize for Pentium(R) 4 processor (DEFAULT)
# -us: append an underscore to external subroutine names
# -vms: enable VMS and DEC statement extensions
# -w95: supress warnings for non-standard fortran (which includes OpenMP) (deprecated 20140129)
 FFLAGS := -diag-disable 5268 -extend_source -fpconstant -fpp2 -implicitnone
# -lsvml needed with -[a]xW
 LDFLAGS_F += -lsvml
 OMP_FLG_F := -openmp
ifeq (${PRC},D)
 FFLAGS += -i4 -r8
else
 FFLAGS += -i4
endif
ifeq (${OPTS},O)
 FFLAGS += -O3 -g -vec_report1 -ftz ${ICC_RCH_ARG}
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
# -C: extensive runtime diagnostics
# -e95 issues warnings for non-standard fortran (which includes OpenMP)
 FFLAGS += -g -C -e95
endif
endif # end Intel (Kai) Fortran Compiler
# Pathscale (QLogic) C++ Compiler
ifeq (pathCC,$(firstword ${CXX}))
# pathcc -show-defaults
# shows that pathcc automatically sets many hardware-specific options
# man -k pathscale for full listing
# Pathscale supports many GCC and PGI-like options including
# -m32
# -m64
# -march=auto Automatically detect and configure target based on current hardware
# -march=opteron
# -march=pentium4
# -march=xeon
# -msse2
# -msse3
# -m3dnow
# man -k pathscale for full listing
# -O2 = -O: Default optimization
# -Ofast = -O3:
 CFLAGS := 
 CPPFLAGS += 
 LDFLAGS += 
 OMP_FLG_C := -apo -mp
ifeq (${OPTS},O)
 CFLAGS += -O3 -g ${PSC_RCH_ARG}
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS += -O2
endif
ifeq (${OPTS},X)
 CFLAGS += -g 
endif
 CXXFLAGS := ${CFLAGS}
endif # end Pathscale (QLogic) C++ Compiler
# Portland Group C++ Compiler
ifeq (pgCC,$(firstword ${CXX}))
# Enable Large File Support (LFS) by default
 CFLAGS := -Mlfs
# Pass kludgy PGI identifier to flag for broken parallel loop variables, crypt()...
 CPPFLAGS := -DPGI_CXX ${CPPFLAGS}
 LDFLAGS += -Mlfs
 OMP_FLG_C := -mp
ifeq (${OPTS},O)
 CFLAGS += -fast ${PGI_RCH_ARG}
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS += 
endif
ifeq (${OPTS},X)
 CFLAGS += -g -Mbounds
endif
 CXXFLAGS := ${CFLAGS}
endif # endif Portland Group C++ Compiler
# Portland Group Fortran Compiler
ifeq (${FC},pgf90)
 FFLAGS := -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign -Ktrap=fp -Mlfs
 OMP_FLG_F := -mp
ifeq (${PRC},D)
 FFLAGS += -Mr8 -Mi4
endif
ifeq (${OPTS},O)
 FFLAGS += -fast
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
 FFLAGS += -g -Mbounds
endif
endif # endif Portland Group Fortran Compiler
# G77 Fortran compiler
ifeq (${FC},g77)
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore
ifeq (${OPTS},O)
 FFLAGS += -O -g
endif
ifeq (${OPTS},D)
 FFLAGS += -g -fdebug-kludge
endif
ifeq (${OPTS},R)
 FFLAGS += -fdebug-kludge
endif
ifeq (${OPTS},X)
 FFLAGS := -g -O -fdebug-kludge -fbounds-check
endif
endif # endif G77 Fortran compiler
ifeq (${OMP},Y)
 CFLAGS += ${OMP_FLG_C}
 CXXFLAGS += ${OMP_FLG_C}
 FFLAGS += ${OMP_FLG_C}
# LD behavior assumes C source code
 LDFLAGS := ${OMP_FLG_C} ${LDFLAGS}
endif # endif OMP
endif
# endif LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, and FREEBSD and MACOSX

ifeq (${PVM_ARCH},NECSX)
ifeq (${OMP},Y)
 OMP_FLG := -Popenmp
endif # endif OMP
CXX := c++
#CC := c++ -Xa
CC := cc
CPP := c++ -E
#CPP := /usr/lib/cpp
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := f90
LD := ld
LDFLAGS += $(TMP_LDFLAGS) $(TMP_LIBS) -lm
ifeq ($(USE_FORTRAN_ARITHMETIC),Y)
# LDFLAGS += -L/lib -lf
endif
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -h2 -hmath vector -hxint
# CFLAGS += -Cvopt -math vector -xint
 FFLAGS = -Cvopt -f3
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS = -g -f3
endif
ifeq (${OPTS},X)
 CFLAGS += -h0 -g  -hstack=nan
# CFLAGS += -Cdebug  -init stack=nan
 FFLAGS = -Cdebug -eR -f3 -Wf"-init stack=nan heap=nan"
endif
MK_DPN = /usr/local/bin/mkdep.perl /usr/lib/cpp # NECSX try this
${MY_DPN_DIR}/%.d : %.c
	${MK_DPN} ${CPPFLAGS} $< | perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@
endif
# endif NECSX

# SGI6, SGI64, SGIMP64
ifneq (${null},$(findstring SGI,${PVM_ARCH}))
ifeq (${OMP},Y)
 OMP_FLG := -mp -mpio
endif # endif OMP
CXX := CC -LANG:std
CC := cc -c99
# 20000302: -w suppresses warnings which will swamp linker
#CXX := g++ -w
#CC := gcc
CPP_TKN_OS +=
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
ifdef $(MIPSPRO_SGI)
# SGIs like dataproc keep omp.h in special location determined by module MIPSpro
 CPPFLAGS := -I$(MIPSPRO_SGI)/usr/include ${CPPFLAGS}
endif
FC := f90 -cpp
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ ${NC_LDFLAGS} ${NC_LIBS}
LEX := flex
LINT := lint
YACC := bison
ifeq (${PVM_ARCH},SGI6)
 GCC_ABI_FLG := -mabi=32
 GCC_LDFLAGS_SZ_SPC := ${GCC_ABI_FLG} -mips3 
 SGI_ABI_FLG := -n32 -mips3 ${OMP_FLG}
else # SGI64, SGIMP64
ifeq (${ABI},64)
 GCC_ABI_FLG := -mabi=64
 GCC_LDFLAGS_SZ_SPC := ${GCC_ABI_FLG} -mips4 -L/usr/local/lib/mabi=64
 SGI_ABI_FLG := -64 -mips4 ${OMP_FLG}
else # ABI=32
 GCC_ABI_FLG := -mabi=32
 GCC_LDFLAGS_SZ_SPC := ${GCC_ABI_FLG} -mips4
 SGI_ABI_FLG := -n32 -mips4 ${OMP_FLG}
endif # endif ABI
endif # endif SGI64, SGIMP64
ifeq (gcc,$(firstword ${CC}))
 LDFLAGS += $(GCC_LDFLAGS_SZ_SPC) ${LDFLAGS}
 CFLAGS := ${GCC_ABI_FLG} -Wall
ifeq (${OPTS},O)
 CFLAGS += -O2 -g
endif
ifeq (${OPTS},R)
 CFLAGS += 
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS += -g -O
endif
 CXXFLAGS := ${CFLAGS}
endif
# endif CC=gcc
ifeq (CC,$(firstword ${CXX}))
# SGI Native C++ headers are in /usr/include/CC
 LDFLAGS := ${SGI_ABI_FLG} ${LDFLAGS}
 CFLAGS := ${SGI_ABI_FLG}
ifeq (${OPTS},O)
 CFLAGS += -O2 -g3
endif
ifeq (${OPTS},R)
 CFLAGS += 
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS += -g -trapuv
endif
 CXXFLAGS := ${CFLAGS}
endif
# endif CC=cc
# Fortran flags
FFLAGS := ${SGI_ABI_FLG} -extend_source
ifeq (${OPTS},O)
 FFLAGS += -O2 -g
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 FFLAGS += -g -check_bounds -trapuv
endif
# end fortran flags
endif
# endif SGI6, SGI64, SGIMP64

ifeq (${UNAMES},SunOS) 
#CXX := CC
#CC := cc
CXX := g++
CC := gcc
CFLAGS := -Wall
CPP := cpp
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := f90 -DHIDE_SHR_MSG
#FFLAGS := -xs -stackvar -e -Qoption f90comp -r8const
FFLAGS := -xs -stackvar -e 
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ ${NC_LDFLAGS} ${NC_LIBS} -lsunmath -lsocket -lnsl
#LDFLAGS += -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq (${ABI},64)
# Sun compiler
ifeq (${CC},cc)
 CFLAGS := -xarch=v9
 FFLAGS := -xarch=v9
 LDFLAGS += -xarch=v9
endif # endif Sun compiler
# GNU Compiler Collection
ifeq (gcc,$(firstword ${CC}))
 CFLAGS += -mabi=64
 FFLAGS += -mabi=64
 LDFLAGS += -mabi=64
endif # endif GNU Compiler Collection
endif # endif ABI
ifeq (${OPTS},O)
 CFLAGS += -O2 -g
 FFLAGS += -fast
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS += -g 
 FFLAGS += -g
# NB: 19980601 -C (range-checking) is not supported by Sun f90
ifeq (${FC},f77)
 FFLAGS += -C
endif
endif
 CXXFLAGS := ${CFLAGS}
endif
# endif SunOS=SUN4SOL2,SUNMP

ifeq (${PVM_ARCH},WIN32)
CXX := g++
CC := gcc
CPP := cpp
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPP_TKN_OS += -Di386
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := g77
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ ${NC_LDFLAGS} ${NC_LIBS}
LDFLAGS += 
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS += -O -g
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
endif
# endif WIN32

# Default to MPICC, MPICXX to CC, CXX and cross your fingers
MPICC := ${CC}
MPICXX := ${CXX}
# Manipulate CC, CPP, CXX as appropriate for MPI-enabled operators
ifeq (${MPI},Y)
ifneq (${null},$(findstring xl,${CC}))
# NB: AIX is not debugged yet
 MPICC := $(subst xlc,mpcc,${CC})
 MPICXX := $(subst xlC,mpCC,${CXX})
endif # endif AIX VA Compiler Collection
ifneq (${null},$(findstring gcc,${CC}))
 MPICC := $(subst gcc,mpicc,${CC})
ifeq (${LAMMPI},Y)
 MPICXX := $(subst g++,mpic++,${CXX})
else
 MPICXX := $(subst g++,mpicxx,${CXX})
endif # !LAMMPI
endif # endif GNU Compiler Collection
ifeq (icpc,$(firstword ${CXX}))
 MPICC := $(subst icc,mpicc,${CC})
ifeq (${LAMMPI},Y)
 MPICXX := $(subst icpc,mpic++,${CXX})
else
 MPICXX := $(subst icpc,mpicxx,${CXX})
endif # !LAMMPI
endif # endif Intel (Kai) C++ Compiler
 CPP := ${MPICC}
endif # endif MPI
# Use MPI modifications, if any, to build all objects
CC := ${MPICC}
CXX := ${MPICXX}

# Disable OpenMP on platforms that automatically support it
ifeq (${OMP},N)
ifneq (${null},$(findstring SGI,${PVM_ARCH}))
 CFLAGS := $(filter-out -mp -mpio,${CFLAGS})
 LDFLAGS := $(filter-out -mp -mpio,${LDFLAGS})
endif # endif SGI
 CPPFLAGS += -U_OPENMP
endif # endif OMP

ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
ifeq (${CCACHE},Y)
# Prefix CC and CXX with ccache
 CC := ccache ${CC}
 CXX := ccache ${CXX}
endif # !CCACHE
endif # !LINUX

ifeq (${STC},Y)
# Created statically linked executable
 LDFLAGS := -static ${LDFLAGS}
endif # endif STC

# Define any remaining variables
libcsz_c++ := ${MY_LIB_DIR}/libcsz_c++
libcsm_c++ := ${MY_LIB_DIR}/libcsm_c++
libnco_c++ := ${MY_LIB_DIR}/libnco_c++

# Default targets
all: ${MDL_BIN_TRG} lib
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes that matter
.SUFFIXES: .cc .c .o .F .d
# Delete targets that were not successfully made
.DELETE_ON_ERROR:
lib_dir :
	- if ! test -d ${MY_LIB_DIR}; then mkdir -p ${MY_LIB_DIR}; fi
bin_dir :
	- if ! test -d ${MY_BIN_DIR}; then mkdir -p ${MY_BIN_DIR}; fi
obj_dir :
	- if ! test -d ${MY_OBJ_DIR}; then mkdir -p ${MY_OBJ_DIR}; fi

# Targets in bin
CCC_OBJ := ${MY_OBJ_DIR}/ccc.o ${MY_OBJ_DIR}/cls.o ${MY_OBJ_DIR}/getopt_bsd.o
#ifeq (${null},$(findstring LINUX,${PVM_ARCH}))
# CCC_OBJ += ${MY_OBJ_DIR}/alloca.o ${MY_OBJ_DIR}/getopt.o ${MY_OBJ_DIR}/getopt1.o
#endif
ccc:	${MY_BIN_DIR}/ccc
${MY_BIN_DIR}/ccc:	lib $(CCC_OBJ)
ifneq (${null},$(findstring MACOSX,${PVM_ARCH}))
	${CXX} -o $@ $(CCC_OBJ) ${LDFLAGS} ${GSL_LIB}
else
	${CXX} -o $@ $(CCC_OBJ) ${LDFLAGS} ${GSL_LIB} -lcrypt
endif
	chmod 755 $@

GSL_OBJ := ${MY_OBJ_DIR}/gsl.o ${MY_OBJ_DIR}/cls.o ${MY_OBJ_DIR}/getopt_bsd.o
gsl:	${MY_BIN_DIR}/gsl
${MY_BIN_DIR}/gsl:	lib $(GSL_OBJ)
	${CXX} -o $@ $(GSL_OBJ) ${LDFLAGS} ${GSL_LIB}
	chmod 755 $@

# i18n needs only -lintl and ABI flags
I18N_OBJ := ${MY_OBJ_DIR}/i18n.o
# 20190420: i18n no longer needs -lintl on Fedora
#I18N_LDFLAGS := -L${NETCDF_ROOT}/lib -lintl
I18N_LDFLAGS := 
ifneq (${null},$(findstring AIX,${PVM_ARCH}))
 I18N_LDFLAGS := -L/usr/local/lib -lintl -qsuppress=1501-245
endif # endif AIX
ifneq (${null},$(findstring SGI,${PVM_ARCH}))
 I18N_LDFLAGS := ${LDFLAGS} -lintl
endif # endif SGI
ifneq (${null},$(findstring icpc,${CXX}))
 I18N_LDFLAGS := 
ifeq (${OMP},Y)
 I18N_LDFLAGS += -openmp -openmp_report0
endif # endif OMP
endif # endif icpc
i18n:	${MY_BIN_DIR}/i18n
${MY_BIN_DIR}/i18n:	${I18N_OBJ} ${MY_SHR_DIR}/locale/es/LC_MESSAGES/i18n.mo
	${CXX} -o $@ ${I18N_OBJ} ${I18N_LDFLAGS}
	chmod 755 $@

mie :
	-cd ${HOME}/$@; ${MAKE} $@

SNW_OBJ := ${MY_OBJ_DIR}/snw.o ${MY_OBJ_DIR}/getopt_bsd.o
snw:	${MY_BIN_DIR}/snw
${MY_BIN_DIR}/snw:	${SNW_OBJ}
#	${CXX} -mdivide-breaks -ftrapping-math -o $@ ${SNW_OBJ} ${LDFLAGS}
	${CXX} -ftrapping-math -o $@ ${SNW_OBJ} ${LDFLAGS}
	chmod 755 $@

SLR_ZEN_NGL_COS_OBJ := ${MY_OBJ_DIR}/slr_zen_ngl_cos.o ${MY_OBJ_DIR}/getopt_bsd.o
slr_zen_ngl_cos:	${MY_BIN_DIR}/slr_zen_ngl_cos
${MY_BIN_DIR}/slr_zen_ngl_cos:	${SLR_ZEN_NGL_COS_OBJ}
#	${CXX} -mdivide-breaks -ftrapping-math -o $@ ${SLR_ZEN_NGL_COS_OBJ} ${LDFLAGS}
	${CXX} -ftrapping-math -o $@ ${SLR_ZEN_NGL_COS_OBJ} ${LDFLAGS}
	chmod 755 $@

# Target-specific variable values syntax TARGET ... : VARIABLE-ASSIGNMENT
tst: LDFLAGS := ${OMP_FLG_LD}
tst: GSL_LIB := ${GSL_LIB}
TST_OBJ := ${MY_OBJ_DIR}/tst.o
tst:	${MY_BIN_DIR}/tst
${MY_BIN_DIR}/tst:	${TST_OBJ}
	${CXX} -o $@ ${TST_OBJ} ${LDFLAGS} ${GSL_LIB}
	chmod 755 $@

bin : ${MDL_BIN_TRG}
bin_cln: 
	rm -f ${MDL_BIN}
strip: 
	- printf "Before stripping...\n";ls -l ${MDL_BIN};strip ${MDL_BIN};printf "After stripping...\n";ls -l ${MDL_BIN}
# Targets in bld
bld_cln :
	cd ${MY_BLD_DIR}; rm -f TAGS
# Targets in dat
dat_cln :
	-cd ${MY_DAT_DIR}; rm -f foo.nc
# Targets in dpn
dpn : ${MDL_DPN}
dpn_cln :
	rm -f ${MDL_DPN}
# Targets in inc
inc : ${MDL_INC}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
lib : inc ${libcsz_c++}.a ${libcsm_c++}.a ${libnco_c++}.a
${libnco_c++}.a :
	-cd ${HOME}/nco/src/nco_c++; ${MAKE} -f Makefile.old lib
${libcsz_c++}.a : ${libcsz_c++}.a(${MY_OBJ_DIR}/bit.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/cln.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/fio.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/mmr.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/mth.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/omp.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/spc_fnc.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/utl.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/vec.o) \
	${libcsz_c++}.a(${MY_OBJ_DIR}/xtr.o)
${libcsm_c++}.a : ${libcsm_c++}.a(${MY_OBJ_DIR}/blm.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/chm.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/flx_sfc.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/htrn_c++.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/idx_rfr.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/idx_rfr_H2O.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/lbl.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/mie_sln.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/mie_MaS99.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/mie_Sor01.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/mie_Wis79.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/tdy.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/aer.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/flp.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/mnr_dst.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/nco.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/phz_fnc.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/psd.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/rt.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/spc_bbd.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/spc_slr.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/ssl.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/sz_grd.o) \
	${libcsm_c++}.a(${MY_OBJ_DIR}/wvl_grd.o)
# NB: Do not place comments in the middle of a rule, or have \ at the end
lib_cln :
	rm -f ${libcsz_c++}.a ${libcsz_c++}.la ${libcsz_c++}.so*
	rm -f ${libcsm_c++}.a ${libcsm_c++}.la ${libcsm_c++}.so*
	rm -f ${libnco_c++}.a ${libnco_c++}.la ${libnco_c++}.so*
# Targets in obj
obj : ${MDL_OBJ}
obj_cln :
	rm -f ${MDL_OBJ}

# Housekeeping
clean : cln
cln : dat_cln dpn_cln obj_cln
debug: dbg
dbg:
	@printf "ABI = ${ABI}\n"
	@printf "AR = ${AR}\n"
	@printf "CC = ${CC}\n"
	@printf "CCACHE = ${CCACHE}\n"
	@printf "CFLAGS = ${CFLAGS}\n"
	@printf "CPP = ${CPP}\n"
	@printf "CPPFLAGS = ${CPPFLAGS}\n"
	@printf "CPP_PTH = ${CPP_PTH}\n"
	@printf "CPP_TKN = ${CPP_TKN}\n"
	@printf "CXX = ${CXX}\n"
	@printf "CXXFLAGS = ${CXXFLAGS}\n"
	@printf "DBG = ${DBG}\n"
	@printf "FC = ${FC}\n"
	@printf "FFLAGS = ${FFLAGS}\n"
	@printf "GCC_RCH_ARG = ${GCC_RCH_ARG}\n"
	@printf "GSL_INC = ${GSL_INC}\n"
	@printf "GSL_LIB = ${GSL_LIB}\n"
	@printf "HOST = ${HOST}\n"
	@printf "HOSTNAME = ${HOSTNAME}\n"
	@printf "ICC_RCH_ARG = ${ICC_RCH_ARG}\n"
	@printf "ICC_ROOT = ${ICC_ROOT}\n"
	@printf "LAMMPI = ${LAMMPI}\n"
	@printf "LAMMPICC = ${LAMMPICC}\n"
	@printf "LAMMPICXX = ${LAMMPICXX}\n"
	@printf "LDFLAGS = ${LDFLAGS}\n"
	@printf "LINUX_CC = ${LINUX_CC}\n"
	@printf "LINUX_CXX = ${LINUX_CXX}\n"
	@printf "LINUX_FC = ${LINUX_FC}\n"
	@printf "MDL_DPN = ${MDL_DPN}\n"
	@printf "MDL_OBJ = ${MDL_OBJ}\n"
	@printf "MDL_PTH = ${MDL_PTH}\n"
	@printf "MDL_SRC = ${MDL_SRC}\n"
	@printf "MPICC = ${MPICC}\n"
	@printf "MPICXX = ${MPICXX}\n"
	@printf "MPICH_CC = ${MPICH_CC}\n"
	@printf "MPICH_CXX = ${MPICH_CXX}\n"
	@printf "MY_BIN_DIR = ${MY_BIN_DIR}\n"
	@printf "MY_BLD_DIR = ${MY_BLD_DIR}\n"
	@printf "MY_DAT_DIR = ${MY_DAT_DIR}\n"
	@printf "MY_DOC_DIR = ${MY_DOC_DIR}\n"
	@printf "MY_DPN_DIR_RX = ${MY_DPN_DIR_RX}\n"
	@printf "MY_ES_DIR = ${MY_ES_DIR}\n"
	@printf "MY_INC_DIR = ${MY_INC_DIR}\n"
	@printf "MY_LIB_DIR = ${MY_LIB_DIR}\n"
	@printf "MY_SHR_DIR = ${MY_SHR_DIR}\n"
	@printf "MY_OBJ_DIR = ${MY_OBJ_DIR}\n"
	@printf "MY_OBJ_DIR_RX = ${MY_OBJ_DIR_RX}\n"
	@printf "NETCDF4 = ${NETCDF4}\n"
	@printf "NETCDF4_ROOT = ${NETCDF4_ROOT}\n"
	@printf "NETCDF_INC = ${NETCDF_INC}\n"
	@printf "NETCDF_LIB = ${NETCDF_LIB}\n"
	@printf "NETCDF_ROOT = ${NETCDF_ROOT}\n"
	@printf "OMP = ${OMP}\n"
	@printf "OPTS = ${OPTS}\n"
	@printf "PGI_RCH_ARG = ${PGI_RCH_ARG}\n"
	@printf "PRC = ${PRC}\n"
	@printf "PSC_RCH_ARG = ${PSC_RCH_ARG}\n"
	@printf "PVM_ARCH = ${PVM_ARCH}\n"
	@printf "STC = ${STC}\n"
	@printf "UNAMES = ${UNAMES}\n"
	@printf "VPATH = ${VPATH}\n"
	@printf "VRS_SNG = ${VRS_SNG}\n"
distclean : dst_cln
dst_cln : cln bld_cln lib_cln
install : inc 
tags : 
	etags *.cc *.c *.com *.F *.F90 *.h *.hh Makefile
uninstall : dst_cln bin_cln inc_cln lib_cln

# Create dependency files only if they will not be immediately deleted
INCLUDE_DPN := TRUE
GOALS_WHICH_DELETE_DEPENDENCY_FILES := cln clean dir distclean dst_cln dpn_cln tags uninstall
ifeq (${null},$(findstring $(MAKECMDGOALS),${GOALS_WHICH_DELETE_DEPENDENCY_FILES}))
 INCLUDE_DPN := TRUE
else
 INCLUDE_DPN := FALSE
endif
ifeq (${INCLUDE_DPN},TRUE)
-include ${MDL_DPN}
endif
