# $Id$

# Purpose: Makefile for C module c

# Copyright (C) 1994--2017 Charlie Zender
# License: GNU General Public License (GPL) Version 3
# See http://www.gnu.org/copyleft/gpl.html for full license text

# Usage:
# make
# cd ~/c;make OPTS=D;cd - # Build with netCDF4
# make --no-builtin-variables <-- Turn off implicit, built-in variables
# make --warn-undefined-variables <-- Fail if Makefile attempts to use undefined variables
# 32-bit ABI on NCAR AIX systems:
# NETCDF_LIB=/usr/local/lib32/r4i4 make --jobs=1 ABI=32 OPTS=D;cd -

# Distribution:
# scp ~/c/Makefile goldhill.cgd.ucar.edu:c
# scp ~/c/Makefile esmf.ess.uci.edu:c

# In rule 
# foobar.o : foo.c bar.c
# $@ is target = foobar.o
# $< is first dependency = foo.c
# $^ is name of all dependencies = foo.c bar.c
# $* is stem of match, and so much more

# Recursively expanded variables are defined with `=':
# Values are installed verbatim and expanded when substitution occurs
# Variable may be used before it is defined and expansion will be correct
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
endif # !ABI
ifndef ${DBG}
# Debugging token N=No (default) Y=Yes
 DBG := N
endif # !DBG
ifndef I18N
 I18N := Y
 MY_SHR_DIR := ${HOME}/share
 MY_ES_DIR := ${MY_SHR_DIR}/locale/es/LC_MESSAGES
 MY_FR_DIR := ${MY_SHR_DIR}/locale/fr/LC_MESSAGES
endif # !I18N
ifdef GSL_INC
 GSL_INC := -I${GSL_INC}
else
 GSL_INC := $(shell gsl-config --cflags)
endif # !GSL_INC
ifdef GSL_LIB
 GSL_LIB := -L${GSL_LIB}
else
 GSL_LIB := $(shell gsl-config --libs)
endif # !GSL_LIB
ifndef MK_DPN
 MK_DPN = ${CPP} -M # NB: Recursive expansion required
 MK_DPN_CXX = ${CXX} -M # NB: Recursive expansion required
endif # !MK_DPN
ifndef NETCDF_ROOT
 NETCDF_ROOT := /usr/local
endif
ifndef NETCDF4_ROOT
 NETCDF4_ROOT := ${NETCDF_ROOT}
endif
ifndef NETCDF_INC
ifdef INC_NCAR
 NETCDF_INC := ${INC_NCAR} # NCAR module path
else # endelse INC_NCAR
 NETCDF_INC := ${NETCDF_ROOT}/include # Directory containing netcdf.h
endif # !INC_NCAR
endif # !NETCDF_INC
ifndef NETCDF_LIB
ifdef LIB_NCAR
 NETCDF_LIB := ${LIB_NCAR} # NCAR module path
else # endelse LIB_NCAR
 NETCDF_LIB := ${NETCDF_ROOT}/lib # Directory containing libnetcdf.a
endif # !LIB_NCAR
endif # !NETCDF_LIB
ifndef ${OMP} # OpenMP
 OMP := Y
endif # !OMP
ifndef OPTS
 OPTS := O
endif # !OPTS
ifndef ${PRC}
# Precision token S=Single D=Double (default) L=Long Double
 PRC := D
endif # !PRC
ifndef PVM_ARCH
 PVM_ARCH := $(shell ${HOME}/sh/pvmgetarch)
endif
ifndef STC
# Created statically linked executable
 STC := N
endif
ifndef UNAMES
 UNAMES := $(shell uname -s)
endif
ifndef USR_TKN
 USR_TKN := -DABORT_ON_ERROR
endif # !USR_TKN
ifndef VRS_SNG
 VRS_SNG := $(shell date +%Y%m%d)
endif # !VRS_SNG

# Derived-tokens based on primary tokens
# These tokens should not be altered by hand
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
# Decide among the plethora of Linux compilers
ifndef LINUX_CXX
# C++ compiler for Linux
 LINUX_CXX := g++
 #LINUX_CXX := como
 #LINUX_CXX := icpc
 #LINUX_CXX := insure
 #LINUX_CXX := pathCC
 #LINUX_CXX := pgCC
endif # !LINUX_CXX
ifndef LINUX_CC
# C compiler for Linux
 LINUX_CC := gcc -std=c99 -pedantic -D_BSD_SOURCE
 #LINUX_CC := como --c99
 #LINUX_CC := icc -std=c99 -D_BSD_SOURCE
 #LINUX_CC := insure
 #LINUX_CC := pathcc -std=c99 -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := pgcc -c9x -D_BSD_SOURCE -D_POSIX_SOURCE
endif # !LINUX_CC
ifndef LINUX_FC
# Fortran compiler for Linux
 LINUX_FC := g95
 #LINUX_FC := gfortran
 #LINUX_FC := ifort
 #LINUX_FC := lf95
 #LINUX_FC := pathf95
 #LINUX_FC := pgf90
endif # !LINUX_CC
ifeq (${PVM_ARCH},LINUXAMD64)
ifndef LINUX_CXX
# GCC 3.3 fails in atomicity.h with -O, -O0, and -O2 so use 3.4 instead 
 LINUX_CXX := g++
endif # !LINUX_CXX
ifndef LINUX_CC
 LINUX_CC := gcc -std=c99 -pedantic -D_BSD_SOURCE
endif # !LINUX_CC
ifndef LINUX_FC
 LINUX_FC := g95
endif # !LINUX_FC
endif # !LINUXAMD64
ifeq (${PVM_ARCH},LINUXALPHA)
 LINUX_FC := g77
endif # !LINUXALPHA
endif # !LINUX
ifeq (${PRC},D)
 USR_TKN := -DPRC_DBL ${USR_TKN}
else
ifeq (${PRC},L)
 USR_TKN := -DPRC_LDB ${USR_TKN}
else
 USR_TKN := -DPRC_FLT ${USR_TKN}
endif # !PRC
endif # !PRC

# NB: CPP_TKN is recursively expanded variable, define only when components are ready
CPP_TKN = ${USR_TKN} -D${PVM_ARCH} -DNO_NETCDF_2 -DVERSION='${VRS_SNG}' -DHOSTNAME='${HOST}' -DUSER='${USER}'

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
MDL_SRC := c.c csz_c.c date_time.c rvr.c slr_crd.c utsname.c makdep.c
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
MDL_BIN_TRG := makdep c hecht rvr slr_crd utsname # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := date_time.c dbg.h # Include targets
MDL_INC_SYM_LNK := # Symbolic links
MDL_INC_STB := ${MDL_INC_TRG} ${MDL_INC_SYM_LNK} # All header files in MY_INC_DIR
MDL_INC := $(addprefix ${MY_INC_DIR}/,${MDL_INC_STB}) # dst_cln removes these files

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
ifneq (${PVM_ARCH},ALPHA)
# HAVE_LONG_LONG
# long long is not ISO-ANSI standard. All compilers except Tru64 support it.
 CPP_TKN_OS += -DHAVE_LONG_LONG
endif # ALPHA

# Works on AIX and AIX46K
ifneq (${null},$(findstring AIX,${PVM_ARCH}))
# 20030804: Always use re-entrant (_r) compilers---Jim Edwards NCAR/IBM 
CC := xlc_r -qlanglvl=extc99
#CC := gcc -std=c99 -pedantic -D_BSD_SOURCE
CXX := xlC_r
# CXX := g++
# fxm: 20040311 should this be cpp? gcc -M?
CPP := xlc -c -qlanglvl=extc99
ifneq (${null},$(findstring xl,${CC}))
# Visual Age compiler headers must occur before g++ headers
 CPP_TKN_OS += -I/usr/vacpp/include
endif # xlC compilers
# Add /usr/local/include for libintl.h explicitly until netCDF is moved there
CPPFLAGS := ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC} -I/usr/local/include
FC := xlf95_r
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
# AIX VA Compiler Collection
ifneq (${null},$(findstring xl,${CC}))
# Visual Age compiler defaults specified in /etc/vac.cfg
# Additional switch to fix compiler warnings on csz.c
# -qmaxmem=num Limit memory used by space intensive optimizations to <num> kilobytes
# -qspill=size Size in B of register allocation spill area, mie needs > 1040 B
 CFLAGS := -qmaxmem=8192 -qspill=2048
 FFLAGS := -NS2000 -qfixed=132
# -bh:5 suppresses annoying messages from xlC linker WARNING: Duplicate symbol: ...
 LDFLAGS += -bh:5
#LDFLAGS += -lxlf90 # Required for linking Fortran objects
ifeq (${OMP},Y)
# -qsmp=noauto : Turn on SMP/OMP code generation but do no automatic parallelization
# -qsmp=omp : Use industry standard OMP without IBM extensions
 OMP_FLG_C := -qsmp=omp
 OMP_FLG_F := -qsmp=omp
else
 CPP_DFN += -U_OPENMP
 FC := xlf95
endif # !OMP
ifeq (${OPTS},O)
# -O : -O3 is safe, -O5 is dangerous
# -qstrict: Ensure that -O3 optimizations do not alter program semantics
# -Q : Inline all appropriate subprograms
 CFLAGS := -O3 -qstrict -Q
 FFLAGS += -O3 -qstrict -Q
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
endif # !ABI
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
endif # !AIX VA Compiler Collection
# GNU Compiler Collection
ifneq (${null},$(findstring gcc,${CC}))
 CFLAGS := -Wall -Wunused
ifeq (${OPTS},O)
 CFLAGS += -O
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
endif # !ABI
 CXXFLAGS := ${CFLAGS}
endif # !GNU Compiler Collection
# -q64: Select 64-bit compiler mode (required for accessing large files)
# -qwarn64: Warn on possible long-to-integer or pointer-to-integer truncation
CXXFLAGS := ${CFLAGS}
ifeq (${OMP},Y)
 CFLAGS += ${OMP_FLG_C}
 CXXFLAGS += ${OMP_FLG_C}
 FFLAGS += ${OMP_FLG_F}
# LD behavior assumes C source code
 LDFLAGS := ${OMP_FLG_C} ${LDFLAGS}
endif # !OMP
endif
# !AIX

ifeq (${PVM_ARCH},ALPHA)
ifeq (${OMP},Y)
 OMP_FLG := -omp
endif # !OMP
CXX := cxx -std strict_ansi
#CXX := cxx -std ansi -D__NO_USE_STD_IOSTREAM
CC := cc
CFLAGS := ${OMP_FLG}
# OSF1 is for CCM compatibility
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} -DOSF1 ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := f90
FFLAGS := -c ${OMP_FLG} -automatic -fpe3
FIXEDFLAGS := -extend_source
FREEFLAGS := -DHIDE_SHR_MSG -free
LD := ld
LDFLAGS += ${OMP_FLG} -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
endif
ifeq (${OPTS},O)
 CFLAGS += -O2 -ansi_alias
 FFLAGS += -O3 -inline speed
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
# !ALPHA

ifeq (${PVM_ARCH},CRAY)
CXX := g++
CC := cc
CPP := cpp
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := f90
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf
#LDFLAGS += -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O
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
# !CRAY

# Works on LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, and FREEBSD and MACOS
ifneq (${null},$(findstring ${PVM_ARCH},LINUXALPHALINUXAMD64LINUXARMFREEBSDMACOS))
CXX := ${LINUX_CXX}
CC := ${LINUX_CC}
CPP := ${CC}
CPP_TKN_OS +=
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := ${LINUX_FC}
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsm_c++ -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf
ifeq (${PVM_ARCH},MACOS)
 LDFLAGS += -lresolv
endif # !MACOS
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
 CFLAGS += -O
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
endif # !Comeau C Compiler
# GNU Compiler Collection or LLVM
# 20140204: gcc and clang should receive identical options
ifeq (gcc,$(firstword ${CC}))
 GCC_OR_CLANG := Y
endif
ifeq (clang,$(firstword ${CC}))
 GCC_OR_CLANG := Y
endif
ifeq (${GCC_OR_CLANG},Y)
# Enable Large File Support (LFS) by default
 CPPFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
 CFLAGS := -Wall -Wunused
# Compilation flags for numerical routines recommended by GSL 1.3 manual, p. 397
# CFLAGS += -Werror -Wall -W -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common -g -O4
# Compilation flags recommended by GSL that I like and use:
# -D_BSD_SOURCE: Support 4.3 BSD Unix extensions to ANSI C (prevents nameser.h warnings)
# -D_POSIX_SOURCE: Support POSIX.1 standard additions to ANSI C (prevents fileno warnings)
# -pedantic: Disallow non-ISO constructs (including type long long) (sometimes useful)
# -Werror: Consider warnings as errors
# -W: Extra warnings, including missing return values, comparison of signed with unsigned
# -Wall: Warn about common programming problems
# -Wmissing-prototypes: Warn if missing prototypes
# -Wshadow: Warn if local variable has same name as other local variable
# -Wpointer-arith: Warn if pointer arithmetic on types without size, e.g., void
# -Wcast-qual: Warn if const qualifier removed from pointer
# -Wcast-align: Warn if casting pointer to type of different size
# -Wwrite-strings: Apply const-qualifier to string constants, die if overwritten
# -fno-common: Prevent global variables from being simultaneously defined in different files
# -g: Put debugging symbols in executable
# -O4: Turn on optimization so uninitialized variables are flagged. Downside: optimizes-out many debugging symbols
# Unused compilation flags that may help:
# -ffast-math: Implement non-IEEE shortcuts
# -funroll-loops: Unroll loops if iterations known at compile time
# -fomit-frame-pointer: Avoid extra frame instructions, reduce debuggability, enabled by -O*
# -m64: Compile in 64-bit mode for x86_64
# -mfpmath=sse:
# Compilation flags recommended by GSL that I do not like and do not use:
# -ansi: Support only strict ANSI C. Equivalent to -std=c89, conflicts with -std=c99
# --no-alias? -fstrict-aliasing
# -Waggregate-return: Warn if functions return aggregates like structures or unions
# -Wconversion: Warn if converting signed to unsigned. Intended for obsolete, non-prototyped code. Triggers fabsf(), sqrtf(), warnings.
# -Wnested-externs: Warn if extern is encountered within function. C only?
# -Wstrict-prototypes: Warn if inconsistent prototypes. C only?
# -Wtraditional: Warn if constructs differ between traditional and ANSI C. C only?
# -Dinline=: inline is not an ANSI keyword, must undefine inline to work with -ansi
# -fshort-enums: Make enums as short as possible, ususally non-int. Do not ever invoke this! This breaks ABI and causes subtle problems
ifeq (${OMP},Y)
ifneq (clang,$(firstword ${CC}))
# 20140526 clang does not recognize/utilize -fopenmp
 OMP_FLG_C := -fopenmp
 OMP_FLG_F := -fopenmp
 LDFLAGS += -lgomp -lpthread
endif # !clang
endif # !OMP
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
 CFLAGS += -g -O
 LDFLAGS += /usr/local/lib/ccmalloc-g++.o -L/usr/local/lib -lccmalloc -ldl
endif
ifneq (${null},$(findstring AMD64,${PVM_ARCH}))
ifeq (${ABI},64)
 CFLAGS += -m64
 FFLAGS += -m64
 LDFLAGS += -m64
endif # !ABI
endif # !LINUXAMD64
 CXXFLAGS := ${CFLAGS}
endif # !GNU Compiler Collection or LLVM
# Intel (Kai) C++ Compiler
ifeq (icc,$(firstword ${CC}))
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
 CPPFLAGS += -I${ICC_ROOT}/include
 LDFLAGS := ${LDFLAGS} -lsvml
 OMP_FLG_C := -openmp -openmp_report0
ifeq (${OPTS},O)
 CFLAGS += -O
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS +=
endif
ifeq (${OPTS},X)
 CFLAGS += -g -Wall
endif
 CXXFLAGS := ${CFLAGS}
endif # !Intel (Kai) C++ Compiler
# Intel (Kai) Fortran Compiler
ifeq (${FC},ifc)
# -e95 issues warnings for non-standard fortran
# -fpp2 necessary, but not sufficient, for  OpenMP
 FFLAGS := -extend_source -implicitnone -vms -e95 -fpp2
# -lVaxlib needed for iargc_, getarg_
# LDFLAGS += -lVaxlib
ifeq (${OMP},Y)
 OMP_FLG_F := -openmp -openmp_report0
endif # !OMP
ifeq (${PRC},D)
 FFLAGS += -i4 -r8 -doubletemps
else
 FFLAGS += -i4
endif
ifeq (${OPTS},O)
 FFLAGS += -O2
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
 FFLAGS += -g -C -e95
endif
endif # !Intel (Kai) Fortran Compiler
# Pathscale (QLogic) C++ Compiler
ifeq (pathcc,$(firstword ${CC}))
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
ifeq (${CC},pgcc)
# Enable Large File Support (LFS) by default
 CFLAGS := -Mlfs
# Pass kludgy PGI identifier to flag for broken C99 designated initializers etc.
 CPPFLAGS := $(filter-out -DHAVE_C99,${CPPFLAGS})
 CPPFLAGS += -DPGI_CC
 LDFLAGS += -Mlfs
 OMP_FLG_C := -mp
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},O)
 CFLAGS += -fast ${PGI_RCH_ARG}
endif
ifeq (${OPTS},R)
 CFLAGS += 
endif
ifeq (${OPTS},X)
 CFLAGS += -g -Mbounds
endif
 CXXFLAGS := ${CFLAGS}
endif # !Portland Group C++ Compiler
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
endif # !Portland Group Fortran Compiler
# G77 Fortran compiler
ifeq (${FC},g77)
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore
ifeq (${OPTS},O)
 FFLAGS += -O
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
endif # !G77 Fortran compiler
ifeq (${OMP},Y)
 CFLAGS += ${OMP_FLG_C}
 CXXFLAGS += ${OMP_FLG_C}
 FFLAGS += ${OMP_FLG_C}
# LD behavior assumes C source code
 LDFLAGS := ${OMP_FLG_C} ${LDFLAGS}
endif # !OMP
endif
# !LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, FREEBSD, MACOS

ifeq (${PVM_ARCH},NECSX)
ifeq (${OMP},Y)
 OMP_FLG := -Popenmp
endif # !OMP
CXX := c++
#CC := c++ -Xa
CC := cc
CPP := c++ -E
#CPP := /usr/lib/cpp
CPPFLAGS := ${CPP_TKN} ${CPP_TKN_OS} -I${NETCDF_INC} ${GSL_INC}
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
# !NECSX

# SGI6, SGI64, SGIMP64
ifneq (${null},$(findstring SGI,${PVM_ARCH}))
ifeq (${OMP},Y)
 OMP_FLG := -mp -mpio
endif # !OMP
#CXX := CC -LANG:std
#CC := cc
# 20000302: -w suppresses warnings which will swamp linker
CXX := g++ -w
CC := gcc
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := f90 -cpp
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${PVM_ARCH},SGI6)
# SGI6
GCC_ABI_FLG := -mabi=32
GCC_LDFLAGS_SZ_SPC := ${GCC_ABI_FLG} -mips3 
SGI_ABI_FLG := -n32 -mips3 ${OMP_FLG}
# !SGI6
else
# SGI64, SGIMP64
GCC_ABI_FLG := -mabi=64
GCC_LDFLAGS_SZ_SPC := ${GCC_ABI_FLG} -mips4 -L/usr/local/lib/mabi=64
SGI_ABI_FLG := -64 -mips4 ${OMP_FLG}
# !SGI64, SGIMP64
endif
ifeq (${CC},gcc)
 LDFLAGS := $(GCC_LDFLAGS_SZ_SPC) ${LDFLAGS}
 CFLAGS := ${GCC_ABI_FLG} -Wall
ifeq (${OPTS},O)
 CFLAGS += -O2
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
# !CC=gcc
ifeq (${CC},cc)
# SGI Native C++ headers are in /usr/include/CC
 LDFLAGS := ${SGI_ABI_FLG} ${LDFLAGS}
 CFLAGS := ${SGI_ABI_FLG}
ifeq (${OPTS},O)
 CFLAGS += -O2
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
endif
# !CC=cc
# Fortran flags
FFLAGS := ${SGI_ABI_FLG} -extend_source
ifeq (${OPTS},O)
 FFLAGS += -O2
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
# !SGI6, SGI64, SGIMP64

ifeq (${UNAMES},SunOS) 
CXX := g++
CC := gcc
CFLAGS := -Wall
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_TKN_OS} ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := f90 -DHIDE_SHR_MSG
#FFLAGS := -xs -stackvar -e -Qoption f90comp -r8const
FFLAGS := -xs -stackvar -e 
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf -lsunmath -lsocket -lnsl
#LDFLAGS += -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O2
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
# !SunOS=SUN4SOL2,SUNMP

ifeq (${PVM_ARCH},WIN32)
C_ONLY := Y
CXX := g++
CC := gcc -std=gnu99
CPP := cpp
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS += ${CPP_TKN} ${CPP_TKN_OS} -Di386 ${CPP_PTH} -I${NETCDF_INC} ${GSL_INC}
FC := g77
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf
LDFLAGS += 
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS += -O
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
endif
# !WIN32

# Link to Electric Fence library, see man efence
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
ifneq (${null},$(findstring ${OPTS},DX))
# Electric Fence works on malloc(), not new(), so it is no good on most C++ programs
# However, ccc contains a malloc() checker which is useful to use with Electric Fence
# DDD/gdb may fail to run code compiled with efence, and complain about libpthread.so.0
# LDFLAGS += -lefence
endif
endif
# !LINUX

ifeq ($(STC),Y)
# Created statically linked executable
 LDFLAGS := -static ${LDFLAGS}
endif # !STC

# Disable OpenMP on platforms that automatically support it
ifeq (${OMP},N)
ifneq (${null},$(findstring SGI,${PVM_ARCH}))
 CFLAGS := $(filter-out -mp -mpio,${CFLAGS})
 LDFLAGS := $(filter-out -mp -mpio,${LDFLAGS})
endif # !SGI
 CPPFLAGS += -U_OPENMP
endif # !OMP

# Define any remaining variables
libcsz_c := ${MY_LIB_DIR}/libcsz_c

# Default targets
all: ${MDL_BIN_TRG} lib ${MY_OBJ_DIR}/date_time.o
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes that matter
.SUFFIXES: .cc .c .o .F .d
# Delete targets that were not successfully made
.DELETE_ON_ERROR:
lib_dir :
	- if ! test -d ${MY_LIB_DIR}; then mkdir ${MY_LIB_DIR}; fi
bin_dir :
	- if ! test -d ${MY_BIN_DIR}; then mkdir ${MY_BIN_DIR}; fi
obj_dir :
	- if ! test -d ${MY_OBJ_DIR}; then mkdir ${MY_OBJ_DIR}; fi

# Targets in bin
C_OBJ := ${MY_OBJ_DIR}/c.o ${MY_OBJ_DIR}/getopt_bsd.o
c:	${MY_BIN_DIR}/c
# Target-specific variable values syntax TARGET ... : VARIABLE-ASSIGNMENT
# Rules begin in leftmost column else interpreted as commands
c : LDFLAGS := $(filter-out -lcsm_c++,${LDFLAGS})
c : LDFLAGS := $(filter-out -lcsz_c++,${LDFLAGS})
c : LDFLAGS := $(filter-out -lnco_c++,${LDFLAGS})
${MY_BIN_DIR}/c:	${C_OBJ} 
	${CC} -o $@ $(C_OBJ) ${LDFLAGS} ${GSL_LIB}
	chmod 755 $@

hello:	${MY_BIN_DIR}/hello
${MY_BIN_DIR}/hello:	${MY_OBJ_DIR}/hello.o
	${CC} -o $@ ${MY_OBJ_DIR}/$@.o ${LDFLAGS}
	chmod 755 $@

RVR_OBJ := ${MY_OBJ_DIR}/rvr.o ${MY_OBJ_DIR}/getopt_bsd.o
rvr:	${MY_BIN_DIR}/rvr
${MY_BIN_DIR}/rvr:	${RVR_OBJ} 
	${CC} -o $@ ${RVR_OBJ} ${LDFLAGS}
	chmod 755 $@

HECHT_OBJ := ${MY_OBJ_DIR}/hecht.o
hecht:	${MY_BIN_DIR}/hecht
${MY_BIN_DIR}/hecht:	${HECHT_OBJ} 
	${CC} -o $@ ${HECHT_OBJ} ${LDFLAGS}
	chmod 755 $@

SLR_CRD_OBJ := ${MY_OBJ_DIR}/slr_crd.o ${MY_OBJ_DIR}/getopt_bsd.o
slr_crd:		${MY_BIN_DIR}/slr_crd
${MY_BIN_DIR}/slr_crd:	${MY_OBJ_DIR}/slr_crd.o
	${CC} -o $@ ${SLR_CRD_OBJ} ${LDFLAGS} -lm
	chmod 755 $@

UTSNAME_OBJ :=	${MY_OBJ_DIR}/utsname.o ${MY_OBJ_DIR}/getopt_bsd.o
utsname:	${MY_BIN_DIR}/utsname
${MY_BIN_DIR}/utsname:	${UTSNAME_OBJ} 
	${CC} -o $@ ${UTSNAME_OBJ} ${LDFLAGS}
	chmod 755 $@

hname_and_IP:	${MY_BIN_DIR}/hname_and_IP
${MY_BIN_DIR}/hname_and_IP:	${MY_OBJ_DIR}/hname_and_IP.o 
	${CC} -o ${MY_BIN_DIR}/hname_and_IP ${MY_OBJ_DIR}/$@.o -lnsl -lsocket -lresolv
	chmod 755 ${MY_BIN_DIR}/hname_and_IP

MAKDEP_OBJ :=	${MY_OBJ_DIR}/makdep.o
# Use system getopt() rather than BSD getopt to eliminate unnecessary dependencies
CPPFLAGS:= -I/usr/include ${CPPFLAGS}
# Target-specific variable values syntax TARGET ... : VARIABLE-ASSIGNMENT
# makdep depends on no libraries
makdep : LDFLAGS := ${null}
makdep:	${MY_BIN_DIR}/makdep
${MY_BIN_DIR}/makdep:	${MAKDEP_OBJ} 
	${CC} -o $@ ${OMP_FLG_C} ${MAKDEP_OBJ} ${LDFLAGS}
	chmod 755 $@

bin : ${MDL_BIN_TRG}
bin_cln: 
	rm -f ${MDL_BIN}
# Targets in bld
bld_cln :
	cd ${MY_BLD_DIR}; rm -f TAGS
clean : cln
cln : dat_cln dpn_cln obj_cln
debug: dbg
dbg:
	@printf "CC = ${CC}\n"
	@printf "CFLAGS = ${CFLAGS}\n"
	@printf "CPP = ${CPP}\n"
	@printf "CPPFLAGS = ${CPPFLAGS}\n"
	@printf "CPP_PTH = ${CPP_PTH}\n"
	@printf "CXX = ${CXX}\n"
	@printf "CXXFLAGS = ${CXXFLAGS}\n"
	@printf "DBG = ${DBG}\n"
	@printf "FC = ${FC}\n"
	@printf "FFLAGS = ${FFLAGS}\n"
	@printf "GSL_LIB = ${GSL_LIB}\n"
	@printf "HOST = ${HOST}\n"
	@printf "HOSTNAME = ${HOSTNAME}\n"
	@printf "LDFLAGS = ${LDFLAGS}\n"
	@printf "MDL_DPN = ${MDL_DPN}\n"
	@printf "MDL_OBJ = ${MDL_OBJ}\n"
	@printf "MDL_PTH = ${MDL_PTH}\n"
	@printf "MDL_SRC = ${MDL_SRC}\n"
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
	@printf "NETCDF_INC = ${NETCDF_INC}\n"
	@printf "NETCDF_LIB = ${NETCDF_LIB}\n"
	@printf "OMP = ${OMP}\n"
	@printf "OPTS = ${OPTS}\n"
	@printf "PRC = ${PRC}\n"
	@printf "PVM_ARCH = ${PVM_ARCH}\n"
	@printf "UNAMES = ${UNAMES}\n"
	@printf "VPATH = ${VPATH}\n"
distclean : dst_cln
dst_cln : cln bld_cln
install : inc 
tags : 
	etags *.cc *.c *.com *.F *.hh *.h Makefile
uninstall : dst_cln bin_cln inc_cln lib_cln
# Targets in dat
dat_cln :
	-cd ${MY_DAT_DIR}; rm -f foo.nc
# Targets in dpn
dpn : ${MDL_DPN}
dpn_cln :
	rm -f ${MDL_DPN}
# Targets in inc
inc : $(MY_INC_TRG)
	cp -f -p ${MDL_INC_TRG} ${MY_INC_DIR}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
lib : $(libcsz_c).a inc
$(libcsz_c).a : $(libcsz_c).a(${MY_OBJ_DIR}/csz_c.o)
lib_cln :
	rm -f $(libcsz_c).a
# Targets in obj
obj_cln :
	rm -f ${MDL_OBJ}

# Create dependency files only if they will not be immediately deleted
INCLUDE_DPN := TRUE
GOALS_WHICH_DELETE_DEPENDENCY_FILES := cln clean distclean dst_cln dpn_cln tags uninstall
ifeq (${null},$(findstring $(MAKECMDGOALS),${GOALS_WHICH_DELETE_DEPENDENCY_FILES}))
 INCLUDE_DPN := TRUE
else
 INCLUDE_DPN := FALSE
endif
ifeq (${INCLUDE_DPN},TRUE)
-include ${MDL_DPN}
endif
