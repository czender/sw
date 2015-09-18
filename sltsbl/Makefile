# $Id$

# Purpose: Makefile for saltation-sandblasting module sltsbl

# Copyright (C) 1994--2005 Charlie Zender
# This software is distributed under the terms of the GNU General Public License
# See http://www.gnu.ai.mit.edu/copyleft/gpl.html for full license text

# Usage:
# Debugging             : cd ~/sltsbl;make OPTS=D sltsbl;cd - 
# Double precision reals: cd ~/sltsbl;make OPTS=D PRC=D sltsbl;cd -
# Single precision reals: cd ~/sltsbl;make OPTS=D PRC=S sltsbl;cd -

# Use separate directory for sltsbl objects because they are often 
# compiled with -r8 -i4 and thus incompatible with, e.g., map
ifndef PVM_ARCH
 PVM_ARCH := $(shell pvmgetarch)
endif
ifndef MY_OBJ_DIR
 MY_OBJ_DIR := ./obj/${PVM_ARCH}
endif
MY_BLD_DIR := .
MY_DAT_DIR := .
MY_DOC_DIR := .
MY_DPN_DIR := ${MY_OBJ_DIR}
MY_HTTP_DIR := /var/www/html/dead
MY_INC_DIR := .
MY_SRC_DIR := .
ifndef MY_BIN_DIR
 MY_BIN_DIR := ./bin/${PVM_ARCH}
endif
ifndef MY_LIB_DIR
 MY_LIB_DIR := ../lib
endif

# Primary tokens which determine build options
# Specify non-default when invoking make, e.g. make DODS=Y
ifndef ${ABI}
# 32 vs. 64 bit ABI: 32=32 bit mode, 64=64 bit mode (default) if available
 ABI := 64
endif # endif ABI
ifndef ${DBG}
# Debugging token N=No (default) Y=Yes
 DBG := N
endif # endif DBG
ifndef DPN_SFX
 DPN_SFX := d
endif # endif DPN_SFX
ifndef MDL_SFX
 MDL_SFX := mod
endif # endif MDL_SFX
ifndef NETCDF_INC
 NETCDF_INC := /usr/local/include 
endif
ifndef NETCDF_LIB
 NETCDF_LIB := /usr/local/lib
endif
ifndef ${OMP} # OpenMP
 OMP := Y
endif # endif OMP
ifndef OPTS
 OPTS := X
endif
ifndef ${PRC}
# Precision token D=Double (default) S=Single
 PRC := S
endif # endif PRC
ifndef STC
# Created statically linked executable
 STC := N
endif
ifndef UNAMES
 UNAMES := $(shell uname -s)
endif
ifndef USR_TKN
 USR_TKN := -DAlG01
endif # endif USR_TKN
ifndef VERSION
 VERSION := $(shell cat ${MY_DOC_DIR}/VERSION)
endif

# Derived-tokens based on primary tokens
# These tokens should not be altered by hand
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
# Decide among the plethora of Linux compilers
ifndef LINUX_CXX
# C++ compiler for Linux
 LINUX_CXX := g++
 #LINUX_CXX := como
 #LINUX_CXX := icc
 #LINUX_CXX := insure
 #LINUX_CXX := pgCC
endif # endif LINUX_CXX
ifndef LINUX_CC
# C compiler for Linux
 LINUX_CC := gcc -std=c99 -pedantic -D_BSD_SOURCE
 #LINUX_CC := como --c99
 #LINUX_CC := icc
 #LINUX_CC := insure
 #LINUX_CC := pgcc
endif # endif LINUX_CC
ifndef LINUX_FC
# Fortran compiler for Linux
 LINUX_FC := g95
 #LINUX_FC := gfortran
 #LINUX_FC := lf95
 #LINUX_FC := ifort
 #LINUX_FC := pgf90
endif # endif LINUX_FC
endif # endif LINUX
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
CPP_TKN = ${USR_TKN} -DDST -DBXM -DCCM -D${PVM_ARCH} ${PRC_TKN}

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := ${null} ${null}
comma := ,
newline := \n

# Directories to search for source files
MDL_PTH := ${MY_SRC_DIR}
# Find all Fortran source files in directory
FIND_FNC = $(wildcard ${dir}/*.F ${dir}/*.F90)
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
MDL_BIN_TRG := sltsbl # Binary targets
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
	-ln -s -f ${NETCDF_INC}/netcdf.${MDL_SFX} .
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f90
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
%.f : %.F
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

# The first LDFLAGS is for typical C programs with netCDF, math, and networking
# The second LDFLAGS enables C/Fortran linking

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
CFLAGS := -qmaxmem=-1
CPP := /lib/cpp -P
# AIX xlf requires wrapping all CPP tokens of type -DFOO and -UFOO in -WF flag
CPP_DFN := ${CPP_TKN}
CPP_DFN_F := -WF,$(subst $(space),$(comma),$(CPP_DFN))
CPP_DFN_F := $(subst $(comma)$(comma),$(comma),$(CPP_DFN_F))
CPPINC := ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
# -qarch=auto : Automatically detect architecture of compiling machine and assume execution on same machine
# -qtune=auto : Optimize executable for architecture detected during compilation
FFLAGS := ${OMP_FLG} -qmaxmem=-1 -qspillsize=2500 -qarch=auto -qtune=auto
FIXEDFLAGS := -qsuffix=f=f:cpp=F -qfixed=132
FREEFLAGS := -qsuffix=f=f90:cpp=F90
LD := ld
# MATCH adds to LDFLAGS: -bmaxdata:2147483648 -bmaxstack:268435456
LDFLAGS := ${OMP_FLG} -L${NETCDF_LIB} -lnetcdf
LDFLAGS := ${OMP_FLG} -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf -bnoquiet
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
 CFLAGS := -O2
# -qstrict: Ensure that -O3 optimizations do not alter program semantics
# -Q : Inline all appropriate subprograms
 FFLAGS += -O3 -qstrict -Q
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -qsigtrap -qextchk -qwarn64
endif
ifeq (${OPTS},X)
 CFLAGS += -g
# -qthreaded: Ensure all internally-generated temporaries have automatic storage
# -qwarn64: Detect truncation of 8 byte integers to 4 bytes
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -qsigtrap -qextchk -qthreaded -qwarn64 -C
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

ifeq (${PVM_ARCH},ALPHA)
ifeq (${OMP},Y)
 OMP_FLG := -omp
endif # endif OMP
CXX := cxx -std strict_ansi
CC := cc
CFLAGS := ${OMP_FLG}
CPP := cpp
CPPFLAGS := -DOSF1 ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f95
FFLAGS := -c ${OMP_FLG} -automatic -fpe3
FIXEDFLAGS := -extend_source
FREEFLAGS := -DHIDE_SHR_MSG -free
LD := ld
LDFLAGS := ${OMP_FLG} -L${NETCDF_LIB} -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
else
 FFLAGS += -r4 -i4
endif # PRC != D
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS += -O3 -inline speed
endif
ifeq (${OPTS},D)
 CFLAGS += -g -check_bounds -check -check_omp
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS := -g -N 132
 FFLAGS := -g -check bounds -check omp_bindings -check overflow -check underflow
endif
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
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f90
LD := ld
LDFLAGS := -L${NETCDF_LIB} -lnetcdf
LDFLAGS += -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O
# 1998/09/09: Added -F -O inline1 since BEE uses them on MATCH
 FFLAGS = -N 132 -F -O inline1
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS = -g -N 132 -F -O inline1
endif
ifeq (${OPTS},X)
 CFLAGS := -g -N 132
 FFLAGS = -g -N 132 -e i -F -O inline1
endif
# 1997/10/21 Added -P to suppress #line # directives on Fortran files
# 1999/05/24 Added -C for cpp to retain comments (else cpp thinks Fortran string concatenator // is a C++ comment and removes it)
${MY_OBJ_DIR}/%.o : %.F
	${CPP} -C -P ${CPPFLAGS} $< > $(patsubst %.F,%.f,$(notdir $<))
	${FC} -c ${FFLAGS} $(patsubst %.F,%.f,$(notdir $<)) 
#	${FC} -e p -c ${FFLAGS} $(notdir $<) 
#	${FC} -e p -Wp"${CPPFLAGS}" -c ${FFLAGS} $(notdir $<) 
	-mv -f $(notdir $@) ${MY_OBJ_DIR}
	rm -f $(patsubst %.F,%.f,$(notdir $<)) 
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
LDFLAGS := -mp -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
# GNU Compiler Collection
ifeq (${CC},gcc)
 CFLAGS := -Wall
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
endif
 CXXFLAGS := ${CFLAGS}
endif
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
# FFLAGS := -ffixed-line-length-132 -fno-second-underscore -Wline-truncation -Waliasing -Wimplicit-none -Wunused-label -Wall
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore -Wimplicit-none -Wunused-label -Wall
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
# FFLAGS += -g -O -fbounds-check -pedantic
 FFLAGS += -g -O -fbounds-check
endif
endif # end GNU Fortran9x Compiler
# Intel (Kai) C Compiler
ifeq (${CC},icc)
 CFLAGS := 
 OMP_FLG_C := -openmp
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
endif
# Intel (Kai) Fortran Compiler
ifeq (${FC},ifort)
# -132: specify 132 column lines for fixed form sources (same as -extend_source)
# -axW: generate generic IA32 code and specialized code exclusively for Pentium 4 processors (requires linking to libsvml.a)
# -fast: enable -xP -O3 -ipo -static
# -fpconstant: extend precision of single precision constants assigned to double precision variables to double precision
# -fpp2: necessary, but not sufficient, for  OpenMP
# -ftz: enable/disable flush denormal results to zero
# -ip: enable single-file IP optimizations
# -ipo: enable multi-file IP optimizations
# -module: specifies path to place .mod files
# -mp: maintain floating point precision (disables some optimizations)
# -nus: do not append an underscore to external subroutine names
# -vec_report1: indicate vectorized loops (DEFAULT)
# -tune pn4: optimize for Pentium(R) 4 processor (DEFAULT)
# -us: append an underscore to external subroutine names
# -vms: enable VMS and DEC statement extensions
# -w95: supress warnings for non-standard fortran (which includes OpenMP)
 FFLAGS := -extend_source -fpconstant -fpp2 -implicitnone -us -w95 -mp
# -lsvml needed with -[a]xW
 LDFLAGS += -lsvml
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
 FFLAGS := -DLinux -DHIDE_SHR_MSG --mod .:${NETCDF_INC} --wide
 OMP_FLG_F := --openmp
ifeq (${PRC},D)
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
 FFLAGS += -O3 --tpp --nchk --ntrace --nsav --nwarn
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
# Portland Group Fortran Compiler
ifeq (${FC},pgf90)
# 20020403: Omit -ktrap=fp so MATCH works in single precision
# FFLAGS := -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign -ktrap=fp
 FFLAGS := -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign
 OMP_FLG_F := -mp
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
else
 FFLAGS += -r4 -i4
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
endif # end Portland Group Fortran Compiler
ifeq (${OMP},Y)
 CFLAGS += ${OMP_FLG_C}
 FFLAGS += ${OMP_FLG_F}
 LDFLAGS := ${OMP_FLG_F} ${LDFLAGS}
endif # endif OMP
endif
# endif LINUX

ifeq (${UNAMES},IRIX64)
CXX := CC
CC := cc
CFLAGS := -64 -mips4
CPP := /usr/lib/cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f90 -cpp
FFLAGS := -64 -mips4 -extend_source -mp -mpio
LD := ld
LDFLAGS := -64 -mp -mips4 -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${PRC},D)
 FFLAGS += -r8 -i4 
else
 FFLAGS += -r4 -i4
endif
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS += -O2
endif
ifeq (${OPTS},R)
 CFLAGS += 
 FFLAGS += 
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS += -g -trapuv
 FFLAGS += -g -check_bounds -ansi -uvar
endif
endif
# endif IRIX64=SGI64,SGIMP64

ifeq (${UNAMES},SunOS) 
CXX := g++
CC := gcc -ansi
CFLAGS := -Wall
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f90
#FFLAGS := -xs -stackvar -e -Qoption f90comp -r8const
# Adding -dalign may be required if netCDF library was compiled with it
# -M: Module directory
FFLAGS := -xs -stackvar -e -M${NETCDF_INC}
LD := ld
LDFLAGS := -L/opt/SUNWspro/lib -lf77compat -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS += -fast -inline=%auto
 LDFLAGS += -fast # 20010420 Yes, -fast tells Sun linker to add extra libraries
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 CFLAGS += -g 
 FFLAGS += -g
ifeq (${FC},f77)
# NB: 19980601 -C (range-checking) is not supported by Sun f90
 FFLAGS += -C
endif
endif
ifeq (${OMP},Y)
 FFLAGS += -openmp # Equivalent to -mp=openmp -explicitpar -stackvar -D_OPENMP
 LDFLAGS += -openmp # 
endif
ifeq (${PRC},D)
# Do not add -r8const flag unless all constants are given precision in the code
# FFLAGS += -r8const
# Using xtypemap appears to cause dead to fail
 FFLAGS += -xtypemap=real:64,double:64,integer:32
ifeq (${FC},f77)
 FFLAGS += -r8 -i4
endif
else
 FFLAGS += -xtypemap=real:32,double:64,integer:32
ifeq (${FC},f77)
 FFLAGS += -r4 -i4
endif
endif
endif
# endif SunOS=SUN4SOL2,SUNMP

ifeq (${DBG},Y)
# Enable Debugging token 
 CPPFLAGS += -DDST_DBG
endif # endif DBG

# Disable OpenMP on platforms that automatically support it
ifeq (${OMP},N)
 FFLAGS := $(filter-out -omp,${FFLAGS})
 FIXEDFLAGS := $(filter-out -omp,${FIXEDFLAGS})
 CFLAGS := $(filter-out -mpio,${CFLAGS})
 FFLAGS := $(filter-out -mpio,${FFLAGS})
 LDFLAGS := $(filter-out -mpio,${LDFLAGS})
 CFLAGS := $(filter-out -mp,${CFLAGS})
 FFLAGS := $(filter-out -mp,${FFLAGS})
 LDFLAGS := $(filter-out -mp,${LDFLAGS})
 ifneq (${UNAMES},SunOS)
  CPPFLAGS += -U_OPENMP
 endif # endif SunOS
endif # endif OMP

ifeq ($(STC),Y)
# Create statically linked executable
 LDFLAGS := -Wl,-Bstatic ${LDFLAGS}
endif # endif STC

# Define any remaining variables

# Default targets
all: dir ${MDL_BIN_TRG}
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .cc .c .${DPN_SFX} .F .F90 .f90 .${MDL_SFX} .o
# Delete targets which were not successfully made
.DELETE_ON_ERROR:
# Target directories which may not exist
dir: bin_dir obj_dir 
bin_dir:
	-mkdir -p ${MY_BIN_DIR}
obj_dir:
	-mkdir -p ${MY_OBJ_DIR}

# Targets in bin
sltsbl : ${MDL_OBJ}
	${FC} -o ${MY_BIN_DIR}/$@ ${MDL_OBJ} ${LDFLAGS}
	chmod 755 ${MY_BIN_DIR}/$@

MAKDEP_OBJ :=	${MY_OBJ_DIR}/makdep.o
# Use system getopt() rather than BSD getopt() to eliminate unnecessary dependencies
CPPFLAGS:= -I/usr/include ${CPPFLAGS}
makdep:	${MY_BIN_DIR}/makdep
${MY_BIN_DIR}/makdep:	$(MAKDEP_OBJ) 
	${CC} -o $@ $(MAKDEP_OBJ)
	chmod 755 $@

bin : ${MDL_BIN_TRG}
	rm -f ${MDL_BIN}
binclean: bin_cln
bin_cln: 
	rm -f ${MDL_BIN}
# Targets in bld
buildclean: bld_cln
bld_cln :
	cd ${MY_BLD_DIR}; rm -f TAGS
# Targets in dat
data:	dat
dat_cln :
	-cd ${MY_DAT_DIR};
# Targets in dpn
depend: dpn
dpn : ${MDL_DPN}
dpn_cln :
	rm -f ${MDL_DPN}
# Targets in dff
diff : dff
dff : 
	diff --exclude='.#*' --exclude='*~' --exclude='*.MOD' --exclude='*.mod' --exclude='.cvsignore' --exclude=ChangeLog  --exclude='CVS' --exclude='foo*' --exclude=Makefile --exclude=README --exclude=TAGS --exclude=TODO ${HOME}/dead ${HOME}/f
# Targets in inc
include: inc
inc : 
	cp -f -p ${MDL_INC_TRG} ${MY_INC_DIR}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
library: lib
lib :
	-cd ${HOME}/f; gmake lib
libclean: lib_cln
lib_cln :
	-cd ${MY_LIB_DIR};
# Targets in mod
mod : ${MDL_OBJ}
	mv -f $(wildcard ${MY_OBJ_DIR}/*.mod) .
mod_cln :
	rm -f $(wildcard *.mod) $(wildcard ${MY_OBJ_DIR}/*.mod)
# Targets in obj
object: obj
obj : ${MDL_OBJ}
objclean: obj_cln
obj_cln :
	rm -f ${MDL_OBJ}

# Housekeeping
clean : cln
cln : dat_cln dpn_cln lib_cln mod_cln obj_cln
debug: dbg
dbg:
	@printf "ABI = ${ABI}\n"
	@printf "AR = ${AR}\n"
	@printf "CXX = ${CXX}\n"
	@printf "CXXFLAGS = ${CXXFLAGS}\n"
	@printf "CC = ${CC}\n"
	@printf "CFLAGS = ${CFLAGS}\n"
	@printf "CPP = ${CPP}\n"
	@printf "CPPFLAGS = ${CPPFLAGS}\n"
	@printf "CPP_PTH = ${CPP_PTH}\n"
	@printf "CPP_TKN = ${CPP_TKN}\n"
	@printf "DBG = ${DBG}\n"
	@printf "FC = ${FC}\n"
	@printf "FFLAGS = ${FFLAGS}\n"
	@printf "FIXEDFLAGS = ${FIXEDFLAGS}\n"
	@printf "FREEFLAGS = ${FREEFLAGS}\n"
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
dst_cln : bld_cln cln lib_cln
DST_MFL := ChangeLog formatlab.m input_10000.txt input_1000.txt input_100.txt lookuptable_final.f Makefile README TODO VERSION # Distribution metafiles
dst : ${DST_MFL} ${SRC_LST}
	/bin/rm -r -f ${DATA}/sltsbl;mkdir -p ${DATA}/sltsbl;/bin/cp -f ${DST_MFL} ${SRC_LST} ${DATA}/sltsbl;
	cd ${DATA};tar cvzf ${DATA}/tmp/sltsbl-${VERSION}.tar.gz ./sltsbl
	scp -p ${DATA}/sltsbl/sltsbl.F90 ${DATA}/tmp/sltsbl-${VERSION}.tar.gz dust.ess.uci.edu:/var/www/html/dead
	-ssh dust.ess.uci.edu "cd /var/www/html/dead;ln -s -f sltsbl-${VERSION}.tar.gz sltsbl.tar.gz"
install : inc 
tags : 
	etags *.cc *.c *.com *.F *.F90 *.f90 *.hh *.h Makefile
uninstall : bin_cln dst_cln inc_cln lib_cln

# Create dependency files only if they will not be immediately deleted
INCLUDE_DPN := TRUE
GOALS_WHICH_DELETE_DEPENDENCY_FILES := cln clean dir distclean dst_cln dpn_cln tags uninstall
ifeq (${null},$(findstring $(MAKECMDGOALS),${GOALS_WHICH_DELETE_DEPENDENCY_FILES}))
 INCLUDE_DPN := TRUE
else
 INCLUDE_DPN := FALSE
endif
ifeq (${INCLUDE_DPN},TRUE)
# Following incorporates dependency files into Makefile rules
-include ${MDL_DPN}
endif

