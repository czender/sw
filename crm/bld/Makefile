# $Id$ -*-makefile-*-

# Purpose: Makefile for CRM
# Allow for precedence of user-modified files, different search paths
# This Makefile requires GNU make, AT&T Berkeley-style make will not parse it

# Known problems: 
# 19990907: GNU make 3.77 has a bug in wildcard function on Solaris
# Solution: Use a different version of GNU make on that platform

# 19990907: Solaris f90 does not support -r8 -i4 so netCDF output is garbage
# Workaround: Use f77 on Solaris and #undef CRM_SRB in misc.h to remove Fortran90 code

# Usage: 
# make dir;make
# Unset environment variables before compiling:
# unset MY_BIN_DIR MY_INC_DIR MY_OBJ_DIR MY_LIB_DIR;make

# Quickie test copies:
# scp dust.ess.uci.edu:crm/bld/Makefile ~/crm/bld/Makefile 
# scp givre.ess.uci.edu:crm/bld/Makefile ~/crm/bld/Makefile 

# Usage:
# make
# cd ~/crm/bld;make OPTS=D NETCDF4=Y;cd - # netCDF4
# cd ~/crm/bld;make OPTS=D MPI=Y NETCDF4=Y;cd - # greenplanet
# cd ~/crm/bld;NETCDF_ROOT='/usr' NETCDF_INC='/usr/include' NETCDF_LIB='/usr/lib' make OPTS=D NETCDF4=N NETCDFF=Y; cd - # netCDF3 from /usr
# 32-bit ABI on NCAR AIX systems:
# cd ~/crm/bld;NETCDF_LIB=/usr/local/lib32/r4i4 make lib_cln cln all ABI=32;cd -
# 64-bit ABI netCDF4 on NCAR AIX systems (bluefire):
# cd ~/crm/bld;CURL_LIB='/contrib/curl/lib' HDF5_ROOT='/contrib/hdf5-1.8.4-patch1_sequential' NETCDF_ROOT='/contrib/netcdf-4.1.1_sequential' SZ_LIB='/contrib/szip/lib' make SZ=Y all;cd -

MY_BIN_DIR := ../bin
MY_BLD_DIR := ../bld
MY_DAT_DIR := ../data
MY_DOC_DIR := ../doc
MY_DPN_DIR := ../obj
MY_INC_DIR := .
ifndef MY_LIB_DIR
 MY_LIB_DIR := ../lib
endif
MY_OBJ_DIR := ${MY_DPN_DIR}
MY_SRC_DIR := ../src

# Primary tokens which determine build options
# Specify non-default when invoking make, e.g. make DODS=Y
ifndef ${ABI}
# 32- vs. 64-bit ABI: 32=32-bit mode, 64=64-bit mode (default) if available
 ABI := 64
endif # endif ABI
ifndef ${DBG}
# Debugging token N=No (default) Y=Yes
 DBG := N
endif # endif DBG
ifndef DPN_SFX
 DPN_SFX := d
endif # endif DPN_SFX
ifndef GCC_RCH_ARG
 GCC_RCH_ARG :=
endif # endif GCC_RCH_ARG
ifndef ICC_RCH_ARG
 ICC_RCH_ARG :=
endif # endif ICC_RCH_ARG
ifndef MDL_SFX
 MDL_SFX := mod
endif # endif MDL_SFX
ifndef ${MPI} # MPI
 MPI := N
endif # endif MPI
ifndef ${NETCDF4} # netCDF4 support
 NETCDF4 := Y
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
ifndef ${NETCDFF} # Link to libnetcdff.a
 NETCDFF := Y
endif # endif NETCDFF
ifndef ${OMP} # OpenMP
 OMP := Y
endif # endif OMP
ifndef OPTS
 OPTS := O
endif
ifndef ${PRC}
# Precision token D=Double (default) S=Single
 PRC := D
endif # endif PRC
ifndef PVM_ARCH
 PVM_ARCH := $(shell ../bld/pvmgetarch)
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
 USR_TKN := 
endif # endif USR_TKN

ifneq (${null},$(findstring AIX,${PVM_ARCH}))
 CCM_TRG_FLG := RS6K
endif
ifneq (${null},$(findstring CRAY,${PVM_ARCH}))
 CCM_TRG_FLG := CRAY
endif
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
 CCM_TRG_FLG := SUN
endif
ifneq (${null},$(findstring SGI,${PVM_ARCH}))
 CCM_TRG_FLG := SGI
endif
ifneq (${null},$(findstring SUN,${PVM_ARCH}))
 CCM_TRG_FLG := SUN
endif
ifeq (${NETCDF4},Y)
 CPPFLAGS += -DENABLE_NETCDF4
 USR_TKN += -DENABLE_NETCDF4
endif # endif NETCDF4

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
CPP_TKN = ${USR_TKN} -D${PVM_ARCH} -DHOSTNAME='${HOST}' -DUSER='${USER}' -D${CCM_TRG_FLG}

# Default build. DAP, NETCDF4, and PNETCDF clauses may overwrite this
NC_LDFLAGS := -L${NETCDF_LIB}
NC_LIBS := -lnetcdf -lcurl

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
# fxm: Temporary kludge: Machines with MPI_ROOT set are assumed to have parallel filesystems for NCO builds
ifdef MPI_ROOT
 HDF5_INC += -I${MPI_ROOT}/include # Directory containing mpi.h
 HDF5_LIB += ${MPI_ROOT}/lib64 # Directory containing libmpi.a
endif # !MPI_ROOT
 NETCDF_INC := ${NETCDF4_ROOT}/include $(shell nc-config --fflags)
 NETCDF_LIB := ${NETCDF4_ROOT}/lib
 NC_LDFLAGS := -L${NETCDF_LIB} $(addprefix -L,${HDF5_LIB})
# NC_LIBS := -lnetcdf -lhdf5_hl -lhdf5 -lz ${NC_LIBS}
 NC_LIBS := $(shell nc-config --flibs)
ifdef MPI_ROOT
 NC_LIBS += -lmpi
endif # !MPI_ROOT
ifdef CURL_LIB
 NC_LDFLAGS += -L${CURL_LIB}
endif # end 
ifeq (${SZ},Y)
ifndef SZ_LIB
 SZ_LIB := /usr/lib # Directory containing libsz.a
endif
 NC_LDFLAGS += -L${SZ_LIB}
 NC_LIBS += -lsz
endif # !SZ
ifneq (${HDF5_INC},${NETCDF_INC})
 NETCDF_INC += -I${HDF5_INC}
endif # end 
endif # endif NETCDF4

# Link to libnetcdff.a
ifeq (${NETCDFF},Y)
 NC_LIBS := -lnetcdff ${NC_LIBS}
endif # endif NETCDFF

# Derived-tokens based on primary tokens
# These tokens should not be altered by hand
# Works on LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, FREEBSD, and WIN32
ifneq (${null},$(findstring ${PVM_ARCH},LINUXALPHALINUXAMD64LINUXARMFREEBSDWIN32))
# Decide among the plethora of Linux compilers
ifndef LINUX_CXX
# C++ compiler for Linux
 LINUX_CXX := g++
 #LINUX_CXX := como
 #LINUX_CXX := icc
 #LINUX_CXX := insure
 #LINUX_CXX := pathCC
 #LINUX_CXX := pgCC
endif # endif LINUX_CXX
ifndef LINUX_CC
# C compiler for Linux
 LINUX_CC := gcc -std=c99 -pedantic -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := como --c99
 #LINUX_CC := icc -std=c99 -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := insure
 #LINUX_CC := pathcc -std=c99 -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := pgcc
endif # endif LINUX_CC
ifndef LINUX_FC
# Fortran compiler for Linux
 #LINUX_FC := g95
 LINUX_FC := gfortran
 #LINUX_FC := lf95
 #LINUX_FC := ifort
 #LINUX_FC := pathf95
 #LINUX_FC := pgf95
endif # endif LINUX_FC
endif # endif LINUX

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := ${null} ${null}
comma := ,
newline := \n

# Directories to search for source files
MDL_PTH := \
	${MY_SRC_DIR}/crm \
	${MY_SRC_DIR}/ccmlsm_share \
	${MY_SRC_DIR}/control \
	${MY_SRC_DIR}/csm_share \
	${MY_SRC_DIR}/dom \
	${MY_SRC_DIR}/eul \
	${MY_SRC_DIR}/physics \
	${MY_SRC_DIR}/srchutil
# Find all C, and Fortran files in a given directory
FIND_FNC = $(wildcard ${dir}/*.c ${dir}/*.F)
# Assemble source files from all directories
SRC_LST = $(foreach dir, ${MDL_PTH},$(FIND_FNC))
# Source file names with directories removed
MDL_SRC := $(notdir $(SRC_LST))
# Dependency list for executable
MDL_OBJ := $(addprefix ${MY_OBJ_DIR}/,$(addsuffix .o, $(basename ${MDL_SRC}))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix ${MY_DPN_DIR}/,$(addsuffix .${DPN_SFX}, $(basename ${MDL_SRC}))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.${DPN_SFX} file
VPATH := $(subst ${space},:,${MDL_PTH})
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,${MDL_PTH},-I${dir})

# Variables having to do with binary executables created by module
MDL_BIN_TRG := crm # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All executables in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := # Include targets
MDL_INC_SYM_LNK := # Symbolic links
MDL_INC_STB := ${MDL_INC_TRG} ${MDL_INC_SYM_LNK} # All header files in MY_INC_DIR
MDL_INC := $(addprefix ${MY_INC_DIR}/,${MDL_INC_STB}) # dst_cln removes these files

# C and C++ rules
${MY_OBJ_DIR}/%.o : %.c
	${CC} ${CPPFLAGS} ${CFLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)
${MY_OBJ_DIR}/%.o : %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)

# Fortran rules
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
 %.${MDL_SFX} : %.F90
	- if test -f ${NETCDF_LIB}/netcdf.${MDL_SFX}; then ln -s -f ${NETCDF_LIB}/netcdf.${MDL_SFX} .; else if test -f ${NETCDF_INC}/netcdf.${MDL_SFX} ; then ln -s -f ${NETCDF_INC}/netcdf.${MDL_SFX} .; fi; fi
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f90
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
%.f : %.F
	${CPP} -P ${CPPFLAGS} -o $@ $<
%.f90 : %.F90
	${CPP} -P ${CPPFLAGS} -o $@ $<

# Dependency rules
# Automatically generate a dependency file for each source file
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.${DPN_SFX}
DPN_GNR := makdep
${MY_DPN_DIR}/%.${DPN_SFX} : %.c
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.cc
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.F
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.F90
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.f
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.f90
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@

# Rules for installing header files
${MY_INC_DIR}/%.com : %.com
	cp -f -p $(notdir $@) $@
${MY_INC_DIR}/%.F90 : %.F90
	cp -f -p $(notdir $@) $@

# First LDFLAGS is for typical C programs with netCDF, math, and networking
# Second LDFLAGS enables C/Fortran linking

# Works on AIX and AIX46K
ifneq (${null},$(findstring AIX,${PVM_ARCH}))
ifeq (${OMP},Y)
 CC := xlc_r -qlanglvl=extc99
 FC := xlf95_r
# -qsmp=noauto : Turn on SMP/OMP code generation but do no automatic parallelization
# -qsmp=omp : Use industry standard OMP without IBM extensions
 OMP_FLG := -qsmp=omp
else
 CC := xlc
 CPP_DFN += -U_OPENMP
 FC := xlf95
endif # endif OMP
CXX := xlC
CFLAGS := ${OMP_FLG} -qmaxmem=-1
CPP := /lib/cpp -P
# AIX xlf requires passing all CPP tokens of type -DFOO and -UFOO with -WF flag
CPP_DFN := ${CPP_TKN}
CPP_DFN_F := -WF,$(subst $(space),$(comma),$(CPP_DFN))
CPP_DFN_F := $(subst $(comma)$(comma),$(comma),$(CPP_DFN_F))
CPPINC := ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
# -qarch=auto : Automatically detect architecture of compiling machine and assume execution on same machine
# -qsrcmsg prints transgressing source line with finger
# -qsuppress=1501-245 : Suppress RLIM_INFINITY memory message due to ulimits
# -qtune=auto : Optimize executable for architecture detected during compilation
FFLAGS := ${OMP_FLG} -qarch=auto -qmaxmem=-1 -qspillsize=2500 -qsuppress=1501-245 -qtune=auto
FIXEDFLAGS := -qsuffix=f=f:cpp=F -qfixed=132
FREEFLAGS := -qsuffix=f=f90:cpp=F90
LD := ld
# Heap is ~eight times larger than stack
# Heap can grow to ~2 GB with 32-bit addressing
# Local arrays or arrays whose size is known at compile time will be allocated from the heap
# Arrays whose size is not known at compile time will always be allocated from the stack
# -bmaxdata: Overrides ulimits stack and heap settings
# MATCH adds to LDFLAGS: -bmaxdata:2147483648 -bmaxstack:268435456
# Obtain maximum allowable: -bmaxdata:0x80000000 -bmaxstack:0x10000000
#LDFLAGS := ${OMP_FLG} -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS} -bnoquiet
LDFLAGS := ${OMP_FLG} -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS} -bnoquiet -qsuppress=1501-245 
LEX := lex
LINT := lint
YACC := yacc
ifeq (${PRC},D)
# -qdpc=e : Promote all floating-point constants to double precision
 FFLAGS += -qrealsize=8 -qintsize=4 -qdpc=e
else
 FFLAGS += -qrealsize=4 -qintsize=4
endif
ifeq (${OPTS},O)
# -qstrict: Ensure that -O3 optimizations do not alter program semantics
# -Q : Inline all appropriate subprograms
 CFLAGS += -O3 -g -qstrict -Q
 FFLAGS += -O3 -g -qstrict -Q
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=inv:ov:zero:en -qfullpath -qsigtrap -qextchk -qwarn64
endif
ifeq (${OPTS},X)
 CFLAGS += -g
# -qcheck=all generates run-time bounds, divide-by-zero, and null-dereference checks
# -qthreaded: Ensure all internally-generated temporaries have automatic storage
# -qwarn64: Detect truncation of 8 byte integers to 4 bytes
# -qflttrap=... # Detect and trap floating point arithmetic problems
# :und: underflow # Underflow NB: Underflow is often OK and too hard to get rid of
# :imp: imprecise # Check for specified exceptions on subprogram entry and exit only
# :inv: invalid # Invalid operations
# :ov: overflow # Overflow
# :zero: zero # Divide-by-zero
# :en: enable # Enable specified exception trapping
# -qfullpath # Show full path back to floating point exceptions
# -qhalt=e stop compiler if error severity equals or exceeds i, w, e, s, u
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=inv:ov:und:zero:en -qfullpath -qsigtrap -qextchk -qthreaded -qwarn64 -C
endif
# -q64: Select 64-bit compiler mode (required for accessing large files)
# -qwarn64: Warn on possible long-to-integer or pointer-to-integer truncation
ifeq (${ABI},64)
 AR := ar -X 64
 CFLAGS += -q64
 FFLAGS += -q64
 LDFLAGS += -q64 
else
 CPPFLAGS += -D_LARGE_FILES
endif # endif ABI
CXXFLAGS := ${CFLAGS}
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif AIX

ifeq (${PVM_ARCH},CRAY)
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS} -lm
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O
 FFLAGS = -N 132
endif
ifeq (${OPTS},D)
 CFLAGS := -g 
 FFLAGS = -g -N 132
endif
ifeq (${OPTS},X)
 CFLAGS := -g -N 132
 FFLAGS = -g -N 132 -e i
endif
# 1997/10/21 Added -P to suppress #line # directives on Fortran files
# 1998/09/07 Added -N to prevent CPP from producing .f that complain about single quotes ('s)
# 1998/09/08 Added -I${NETCDF_INC} to ${FC} for CCM code (e.g., initext.F) that uses Fortran includes
# NB: -N and -P are not part of CPPFLAGS because they should only used for Fortran files
${MY_OBJ_DIR}/%.o : %.F
	${CPP} -N -P ${CPPFLAGS} $< > ${MY_OBJ_DIR}/$(patsubst %.F,%.f,$(notdir $<))
	${FC} -I${NETCDF_INC} ${FFLAGS} -b $@ ${MY_OBJ_DIR}/$(patsubst %.F,%.f,$(notdir $<))
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} $<
	mv -f $(notdir $@) ${MY_OBJ_DIR}
endif
# endif CRAY

# Works on LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, and FREEBSD
ifneq (${null},$(findstring ${PVM_ARCH},LINUXALPHALINUXAMD64LINUXARMFREEBSD))
CXX := ${LINUX_CXX}
CC := ${LINUX_CC}
CPP := ${CC}
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := ${LINUX_FC}
LD := ld
LDFLAGS := -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS}
LEX := flex
LINT := lint
YACC := bison
# GNU Compiler Collection
ifeq (${CC},gcc)
 CFLAGS := -Wall
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
endif
 CXXFLAGS := ${CFLAGS}
endif # end GNU Compiler Collection
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
 FFLAGS += -g -O -fdebug-kludge -fbounds-check
endif
endif # end G77 Fortran compiler
# GNU Fortran9x Compiler
ifeq (${FC},g95)
# -DF2K causes sng to use compiler's F2K-compatible command line argument routines
 CPPFLAGS += -DF2K
# FFLAGS := -ffixed-line-length-132 -fno-second-underscore -Wline-truncation -Waliasing -Wimplicit-none -Wunused-label -Wall
# FFLAGS := -ffixed-line-length-132 -fno-second-underscore -Wimplicit-none -Wunused-label -Wall
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore -Wall
ifeq (${PRC},D)
 FFLAGS += -r8
else
 FFLAGS += 
endif
ifeq (${OPTS},O)
 FFLAGS += -O -g
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 FFLAGS += -g -O -fbounds-check -pedantic
endif
endif # end GNU Fortran9x Compiler
# GNU Fortran Compiler
ifeq (gfortran,$(firstword ${FC}))
 FFLAGS := -fcray-pointer -ffixed-line-length-132 -Wall
ifeq (${OMP},Y)
 OMP_FLG_C := -fopenmp
 OMP_FLG_F := -fopenmp
 LDFLAGS += -lgomp -lpthread
endif # endif OMP
ifeq (${PRC},D)
 FFLAGS += -fdefault-real-8
else
 FFLAGS += 
endif
ifeq (${OPTS},O)
 FFLAGS += -O -g
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 FFLAGS += -g -O -fbounds-check -pedantic -Wconversion -Werror -Wimplicit-interface -Wunderflow
endif
endif # end GNU Fortran Compiler
# Intel (Kai) C Compiler
ifeq (${CC},icc)
 CFLAGS := 
 OMP_FLG_C := -openmp
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
endif # end Intel (Kai) C Compiler
ifeq (${FC},ifort)
# -132 specify 132 column lines for fixed form sources (same as -extend_source)
# -fast: enable -xP -O3 -ipo -static
# -fpp2 necessary, but not sufficient, for  OpenMP
# -ftz enable/disable flush denormal results to zero
# -ip enable single-file IP optimizations
# -ipo enable multi-file IP optimizations
# -module specifies path to place .mod files
# -mp maintain floating point precision (disables some optimizations) (deprecated 20100708)
# -axW generate generic IA32 code and specialized code exclusively for Pentium 4 processors
# -vec_report[n] n=1 indicate vectorized loops (DEFAULT)
# -tune pn4 optimize for Pentium(R) 4 processor (DEFAULT)
# -w95 supresses warnings for non-standard fortran (which includes OpenMP)
 FFLAGS := -extend_source -fpconstant -fpp2 -implicitnone -w95
# -lPEPCF90 needed for iargc_, getarg_
# LDFLAGS += -lPEPCF90
 OMP_FLG_F := -openmp
ifeq (${PRC},D)
 FFLAGS += -i4 -r8
else
 FFLAGS += -i4
endif
ifeq (${OPTS},O)
 FFLAGS += -O3 -g -axW -vec_report1 -ftz
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
# Lahey-Fujitsu Fortran Compiler
ifeq (${FC},lf95)
# --wide: Allow wide source format for fixed form
# NB: Convert to/from big-endian with runtime option: a.out -Wl,-T[unit_nbr]
 FFLAGS := -DLinux -DHIDE_SHR_MSG --mod .:${NETCDF_INC} --wide
 OMP_FLG_F := --openmp
ifeq (${PRC},D)
# Passing --dbl confuses resolution of overloaded float/double functions
# Ensure such functions are enclosed in #ifndef PRC_DBL
 FFLAGS += --dbl
else
 FFLAGS += 
endif
ifeq ($(STC),Y)
 LDFLAGS := --staticlink ${LDFLAGS}
endif # endif STC
ifeq (${OPTS},O)
# --tpp: Pentium Pro/Pentium II code generation
# --[n]sav: [Do not] save local variables
 FFLAGS += -O3 -g --tpp --nchk --ntrace --nsav --nwarn
endif
ifeq (${OPTS},D)
 FFLAGS += -g --chk a,e,s --pca --trace --trap
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
 FFLAGS += -g --chk a,e,s,u,x --chkglobal --pca --trace --trap
endif
endif # end Lahey-Fujitsu Fortran Compiler
# Portland Group C++ Compiler
ifeq (${CXX},pgCC)
# Enable Large File Support (LFS) by default
 CFLAGS := -Mlfs
 LDFLAGS += -Mlfs
 OMP_FLG := -mp
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},O)
 CFLAGS += -fast
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
 FFLAGS += -r8 -i4
else
 FFLAGS += -r4 -i4
endif
ifeq ($(STC),Y)
 LDFLAGS := -Wl,-Bstatic ${LDFLAGS}
endif # endif STC
ifeq (${OPTS},O)
 FFLAGS += -g -fast
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 FFLAGS += -g -Mbounds
endif
endif # end Portland Group Fortran Compiler
ifeq (${OMP},Y)
 CFLAGS += ${OMP_FLG_C}
 FFLAGS += ${OMP_FLG_F}
 LDFLAGS := ${OMP_FLG_F} ${LDFLAGS}
endif # endif OMP
endif
# endif LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, FREEBSD

ifeq (${PVM_ARCH},RS6K)
C++ := g++
CC := gcc -ansi
CPP := /lib/cpp -P
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := xlf
LD := ld
LDFLAGS := -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS}
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O2
 PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS := -O -NS2000 -qfixed=132
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS := -g -NS2000 -qfixed=132
endif
 FFLAGS += -qREALSIZE=8 -qINTSIZE=4
${MY_OBJ_DIR}/%.o : %.F
	${CPP} ${CPPFLAGS} $< ${MY_OBJ_DIR}/$(basename $<).f 
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) ${MY_OBJ_DIR}/$(basename $<).f
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif RS6K

# SGI6, SGI64, SGIMP64
ifneq (${null},$(findstring SGI,${PVM_ARCH}))
C++ := CC -LANG:std
CC := cc
# 20000302: -w suppresses warnings which will swamp linker
#C++ := g++ -w
#CC := gcc
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
ifdef $(MIPSPRO_SGI)
# SGIs like dataproc keep omp.h in special location determined by module MIPSpro
 CPPFLAGS := -I$(MIPSPRO_SGI)/usr/include ${CPPFLAGS}
endif
FC := f90 -cpp
LD := ld
LEX := flex
LINT := lint
YACC := bison
ifeq (${PVM_ARCH},SGI6)
# SGI6
GCC_ABI_FLG := -mabi=32
GCC_LDFLAGS_SZ_SPC := ${GCC_ABI_FLG} -mips3 
SGI_ABI_FLG := -n32 -mips3
# endif SGI6
else
# SGI64, SGIMP64
GCC_ABI_FLG := -mabi=64
GCC_LDFLAGS_SZ_SPC := ${GCC_ABI_FLG} -mips4 -L/usr/local/lib/mabi=64
SGI_ABI_FLG := -64 -mips4
#SGI_ABI_FLG += -mp -mpio # Add multithreading
# endif SGI64, SGIMP64
endif
ifeq (${CC},gcc)
 LDFLAGS += $(GCC_LDFLAGS_SZ_SPC) -lm -L${NETCDF_LIB} -lnetcdf
ifeq (${OPTS},O)
 CFLAGS += ${GCC_ABI_FLG} -O2 -Wall
endif
ifeq (${OPTS},R)
 CFLAGS += ${GCC_ABI_FLG} -Wall
endif
ifeq (${OPTS},D)
 CFLAGS += ${GCC_ABI_FLG} -g -Wall
endif
ifeq (${OPTS},X)
 CFLAGS += ${GCC_ABI_FLG} -g -O -Wall
endif
 CXXFLAGS := ${CFLAGS}
endif
# endif CC=gcc
ifeq (${CC},cc)
 LDFLAGS += ${SGI_ABI_FLG} -lm -L${NETCDF_LIB} -lnetcdf
ifeq (${OPTS},O)
 CFLAGS += ${SGI_ABI_FLG} -O2
endif
ifeq (${OPTS},R)
 CFLAGS += ${SGI_ABI_FLG}
endif
ifeq (${OPTS},D)
 CFLAGS += ${SGI_ABI_FLG} -g
endif
ifeq (${OPTS},X)
 CFLAGS += ${SGI_ABI_FLG} -g -trapuv
endif
endif
# endif CC=cc
# Fortran flags
ifeq (${OPTS},O)
 FFLAGS := ${SGI_ABI_FLG} -O2 -extend_source
endif
ifeq (${OPTS},R)
 FFLAGS := ${SGI_ABI_FLG} -extend_source
endif
ifeq (${OPTS},D)
 FFLAGS := ${SGI_ABI_FLG} -g -extend_source
endif
ifeq (${OPTS},X)
 FFLAGS := ${SGI_ABI_FLG} -g -extend_source -check_bounds -trapuv
endif
 FFLAGS += -r8 -i4
# end fortran flags
endif
# endif SGI6, SGI64, SGIMP64

ifeq (${PVM_ARCH},SUN4)
C++ := g++
CC := acc
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS}
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O2
 FFLAGS := -fast -e -Nl99
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS := -g -e -Nl99
endif
 FFLAGS += -r8 -i4
endif
# endif SUN4

ifeq (${PVM_ARCH},SUN4SOL2)
C++ := g++
CC := gcc -ansi
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS}
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O2
 FFLAGS := -O -e
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS := -g -e
endif
ifeq (${OPTS},X)
 CFLAGS := -g 
 FFLAGS := -g -e
# NB: 1998/06/01 -C (range-checking) is not supported by Sun f90
ifeq (${FC},f77)
 FFLAGS += -C -r8 -i4
endif
endif
# NB: 1999/09/07 -r8 -i4 are not supported by Sun f90, Sun f77 now uses -dbl
ifeq (${FC},f77)
 FFLAGS += -r8 -i4
endif
endif
# endif SUN4SOL2

ifeq (${PVM_ARCH},SUNMP)
C++ := g++
CC := gcc -ansi
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS}
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O2
 FFLAGS := -O -e
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS := -g -e
endif
ifeq (${OPTS},X)
 CFLAGS := -g 
 FFLAGS := -g -e
# 1998/06/01 -C (range-checking) is not supported by Sun f90
ifeq (${FC},f77)
 FFLAGS += ${FFLAGS} -C
endif
endif
# 1999/09/07 -r8 -i4 are not supported by Sun f90
ifeq (${FC},f77)
 FFLAGS += -r8 -i4
else
 FFLAGS += -xtypemap=real:64,double:64,integer:32
endif
endif
# endif SUNMP

# Link to Electric Fence library, see man efence
ifeq (${PVM_ARCH},LINUX)
ifneq (${null},$(findstring ${OPTS},DX))
 LDFLAGS += -lefence
endif
endif
# endif LINUX

# Define any remaining variables

# Default targets
all : ${MDL_BIN_TRG}
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug 
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .cc .c .o .F .d

# Targets in bin
CRM_OBJ := ${MDL_OBJ}
crm:	$(CRM_OBJ)
	${FC} -o ${MY_BIN_DIR}/$@ ${MDL_OBJ} ${LDFLAGS}
	chmod 755 ${MY_BIN_DIR}/$@
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
	@printf "ABI = ${ABI}\n"
	@printf "AR = ${AR}\n"
	@printf "CC = ${CC}\n"
	@printf "CFLAGS = ${CFLAGS}\n"
	@printf "CPP = ${CPP}\n"
	@printf "CPPFLAGS = ${CPPFLAGS}\n"
	@printf "CPPINC = ${CPPINC}\n"
	@printf "CPP_DFN = ${CPP_DFN}\n"
	@printf "CPP_DFN_F = ${CPP_DFN_F}\n"
	@printf "CPP_PTH = ${CPP_PTH}\n"
	@printf "CPP_TKN = ${CPP_TKN}\n"
	@printf "CURL_LIB = ${CURL_LIB}\n"
	@printf "CXX = ${CXX}\n"
	@printf "CXXFLAGS = ${CXXFLAGS}\n"
	@printf "LDFLAGS = ${LDFLAGS}\n"
	@printf "LEX = $(LEX)\n"
	@printf "MDL_DPN = ${MDL_DPN}\n"
	@printf "MDL_OBJ = ${MDL_OBJ}\n"
	@printf "MDL_PTH = ${MDL_PTH}\n"
	@printf "MDL_SRC = ${MDL_SRC}\n"
	@printf "MY_BIN_DIR = ${MY_BIN_DIR}\n"
	@printf "MY_BLD_DIR = ${MY_BLD_DIR}\n"
	@printf "MY_DAT_DIR = ${MY_DAT_DIR}\n"
	@printf "MY_DOC_DIR = ${MY_DOC_DIR}\n"
	@printf "MY_DPN_DIR = ${MY_DPN_DIR}\n"
	@printf "MY_INC_DIR = ${MY_INC_DIR}\n"
	@printf "MY_LIB_DIR = ${MY_LIB_DIR}\n"
	@printf "MY_OBJ_DIR = ${MY_OBJ_DIR}\n"
	@printf "MY_SRC_DIR = ${MY_SRC_DIR}\n"
	@printf "NETCDF_INC = ${NETCDF_INC}\n"
	@printf "NETCDF_LIB = ${NETCDF_LIB}\n"
	@printf "PVM_ARCH = ${PVM_ARCH}\n"
	@printf "SRC_LST = $(SRC_LST)\n"
	@printf "VPATH = ${VPATH}\n"
	@printf "YACC = $(YACC)\n"
distclean : dst_cln
dst_cln : cln bld_cln
install : inc 
tags : 
	etags $(SRC_LST) Makefile
uninstall : dst_cln bin_cln inc_cln lib_cln
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
lib : 
	-cd ${MY_LIB_DIR}
lib_cln :
	-cd ${MY_LIB_DIR}
# Targets in obj
obj : ${MDL_OBJ}
obj_cln :
	rm -f ${MDL_OBJ}

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

