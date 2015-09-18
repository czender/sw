# $Id$

# Purpose: Makefile for solar spectrum module slr_spc

MY_BLD_DIR := .
MY_DAT_DIR := ${DATA}/aca
MY_DOC_DIR := .
MY_DPN_DIR := ${MY_OBJ_DIR}
MY_INC_DIR := ${HOME}/include
ifndef MY_BIN_DIR
 MY_BIN_DIR := ../bin
endif
ifndef MY_LIB_DIR
 MY_LIB_DIR := ../lib
endif
ifndef MY_OBJ_DIR
 MY_OBJ_DIR := ../obj
endif
ifndef OPTS
 OPTS := O
endif
ifndef ${PRC}
# Precision token D=Double (default) S=Single
 PRC := S
endif # endif PRC
ifndef PVM_ARCH
 PVM_ARCH := $(shell ${MY_BLD_DIR}/pvmgetarch)
endif
ifndef STC
# Created statically linked executable
 STC := N
endif
ifndef UNAMES
 UNAMES := $(shell uname -s)
endif
ifndef VRS_SNG
 VRS_SNG=$(shell date +%Y%m%d)
endif # endif VRS_SNG
USR_TKN := -DABORT_ON_ERROR
ifeq (${PRC},D)
 USR_TKN := -DPRC_DBL ${USR_TKN}
else
 USR_TKN := -DPRC_FLT ${USR_TKN}
endif # endif PRC
# NB: CPP_TKN must be defined as a recursively expanded variable
CPP_TKN = ${USR_TKN} -D${PVM_ARCH} -DNO_NETCDF_2 -DVERSION='${VRS_SNG}' -DHOSTNAME='${HOST}' -DUSER='${USER}'

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := ${null} ${null}
comma := ,
newline := \n
# We replace `/' by `cszzsc', call perl, then replace `cszzxc' by `/'
# Unique character(s) to substitute for ${slash} before passing to perl regex
slash_rx := cszzsc
# Unique character(s) to replace by ${slash_rx} before passing to perl regex
slash := /
MY_OBJ_DIR_RX := $(subst ${slash},${slash_rx},${MY_OBJ_DIR})
MY_DPN_DIR_RX := $(subst ${slash},${slash_rx},${MY_DPN_DIR})

# Source file names with directories removed
MDL_SRC := slr_spc.cc flx_slr_frc.f
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
MDL_BIN_TRG := slr_spc # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := # Include targets
MDL_INC_SYM_LNK := # Symbolic links
MDL_INC_STB := ${MDL_INC_TRG} ${MDL_INC_SYM_LNK} # All header files in MY_INC_DIR
MDL_INC := $(addprefix ${MY_INC_DIR}/,${MDL_INC_STB}) # dst_cln removes these files

# Redefine default C and C++ pattern rules
${MY_OBJ_DIR}/%.o : %.c
	${CC} ${CPPFLAGS} ${CFLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)
${MY_OBJ_DIR}/%.o : %.cc
	${C++} ${CPPFLAGS} ${C++FLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)

# Default Fortran pattern rules: CRAY and RS6K must override these rules
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<

# Automatically generate a dependency file for each source file.
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.d
# Linux gcc may return an extra `.F' on Fortran names, e.g., `hello.F.o: hello.F'
# (.F)? gets rid of this extra `.F'
${MY_DPN_DIR}/%.d : %.F
# Following command makes, e.g., f.d begin "f.o f.d : f.F ..."
# Since f.o is not preceded by ${MY_OBJ_DIR}, objects are not recompiled when sources are touched.
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.F\.o/$*.o $@/g;' > $@
# Following command makes, e.g., f.d begin "/home/zender/obj/LINUX/f.o f.d : f.F ..."
# This works fairly well, but is a hack.
# First pattern substitutes MY_OBJ_DIR_RX, which has placeholders for slashes.
# Second pattern substitutes slashes for the placeholders.
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*(\.F)?\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@
# Following command makes, e.g., f.d begin "${MY_OBJ_DIR}/f.o f.d : f.F ..."
# This would be the ideal command but I can't get the dollar sign to show up.
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.F\.o/\${dollar}MY_OBJ_DIR\/$*.o $@/g;' > $@

# fxm: Experimental dependency generator for *.F90 files
# 20000705: gcc preprocessor does not understand *.F90 files so CCM uses makdep instead
# This works:
# gcc -E -M -DLINUX -I. -I/home/zender/include -I/usr/local/include f.F
# This does not:
# gcc -E -M -DLINUX -I. -I/home/zender/include -I/usr/local/include f.F90
# -x option to gcc explicitly specifies language
# gcc -x f -E -M -DLINUX -I. -I/home/zender/include -I/usr/local/include f.F90
${MY_DPN_DIR}/%.d : %.F90
	cp -f $< $*.F; ${CPP} -M ${CPPFLAGS} $*.F | perl -p -e 's/$*(\.F)?\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/$*.F/$*.F90/g;s/${slash_rx}/\${slash}/g' > $@;rm -f $*.F

${MY_DPN_DIR}/%.d : %.f
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*(\.f)?\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.c
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/$*.o $@/g;' > $@
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.cc
# NB: Use ${C++} rather than ${CPP} on C++ files for now because, e.g., SUNMP cpp does not 
# preprocess .cc files quite correctly
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/$*.o $@/g;' > $@
	${C++} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

# First LDFLAGS is for typical C programs with netCDF, math, and networking
# Second LDFLAGS enables C/Fortran linking
# Works on AIX and AIX46K
ifneq (,$(findstring AIX,${PVM_ARCH}))
C++ := xlC
CC := xlc
# Additional switch to fix compiler warnings on csz.c
# -qmaxmem=num Limit the amount of memory used by space intensive optimizations to <num> kilobytes
CFLAGS := -qmaxmem=8192
CPP := xlc -c
# Add /usr/local/include for libintl.h explicitly until netCDF is moved there
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I/usr/local/include
FC := xlf90
FFLAGS := -NS2000 -qfixed=132
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf
LDFLAGS += -lxlf90 # Required for linking Fortran objects
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O
 FFLAGS += -O
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g
endif
ifeq (${OPTS},X)
# -qflttrap generates instructions for floating point exceptions
# -qidirfirst uses headers found in -I directories first
# -qlonglong allow long long integers
# -qmakedep creates .d file
# -qwarn64 check for possible long-to-integer or pointer-to-integer truncation
# -qhalt=e stop compiler if error severity equals or exceeds i, w, e, s, u
 CFLAGS += -g -qflttrap -qidirfirst -qlonglong -qwarn64 -qcheck=all -qhalt=s
 FFLAGS += -g
endif
C++FLAGS := ${CFLAGS}
# Additional flags for AIX:
# -M Generate information to be included in a "make" description file; output goes to .u file
# -c Do not send object files to the linkage editor
# -P Preprocess but do not compile; output goes to .i file
# Using -P causes additional warning messages about lm 
# Not using -P causes *.o files to be created twice
${MY_DPN_DIR}/%.d : %.c
	${CPP} -M ${CPPFLAGS} $< ;perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' $*.u > $@ ; \
	rm -f $*.i $*.o $*.u;
${MY_DPN_DIR}/%.d : %.cc
	${CPP} -M ${CPPFLAGS} $< ;perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' $*.u > $@ ; \
	rm -f $*.i $*.o $*.u;
endif
# endif AIX

ifeq (${PVM_ARCH},ALPHA)
ifeq (${OMP},Y)
 OMP_FLG := -omp
endif # endif OMP
C++ := cxx -std strict_ansi
#C++ := cxx -std ansi -D__NO_USE_STD_IOSTREAM
CC := cc
CFLAGS := ${OMP_FLG}
# OSF1 is for CCM compatibility
CPPFLAGS += ${CPP_TKN} -DOSF1 ${CPP_PTH} -I${NETCDF_INC}
FC := f90
FFLAGS := -c ${OMP_FLG} -automatic -fpe3
FIXEDFLAGS := -extend_source
FREEFLAGS := -DHIDE_SHR_MSG -free
LD := ld
LDFLAGS += ${OMP_FLG} -L${MY_LIB_DIR} -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf
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
C++FLAGS := ${CFLAGS}
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif ALPHA

ifeq (${PVM_ARCH},CRAY)
NETCDF_INC := /usr/local/include
NETCDF_LIB := /usr/local/lib
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf
LDFLAGS += -L/lib -lf
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
# 97/10/21 Added -P to suppress #line # directives on Fortran files
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

ifeq (${PVM_ARCH},LINUX)
NETCDF_INC := /usr/local/include
NETCDF_LIB := /usr/local/lib
C++ := g++
# 1998/05/17 As of gcc-2.8.0 using -ansi switch breaks because resolv.h structure have changed
#CC := gcc -ansi
CC := gcc
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := pgf90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${CC},gcc)
ifeq (${OPTS},O)
 CFLAGS := -O -Wall
endif
ifeq (${OPTS},D)
 CFLAGS := -g -Wall
endif
ifeq (${OPTS},R)
 CFLAGS := -Wall
endif
ifeq (${OPTS},X)
 CFLAGS := -g -O -Wall
endif
 C++FLAGS := ${CFLAGS}
endif
ifeq (${FC},pgf90)
# Add PGI include directory to FFLAGS rather than CPPFLAGS so it is only searched by Fortran objects
# 2001008 fxm: this means *.d dependency files will not contain PGI headers
FFLAGS += -I$(PGI)/linux86/include 
# Libraries needed by g++ to link PGI objects
LDFLAGS += -L$(PGI)/linux86/lib -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lpgc
ifeq (${OPTS},O)
 FFLAGS += -fast -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign
endif
ifeq (${OPTS},D)
 FFLAGS += -g -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign
endif
ifeq (${OPTS},R)
 FFLAGS += -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign
endif
ifeq (${OPTS},X)
 FFLAGS += -g -Mbounds -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign
endif
endif
ifeq (${FC},g77)
LDFLAGS += -lg2c
ifeq (${OPTS},O)
 FFLAGS += -O -ffixed-line-length-132 -fno-second-underscore
endif
ifeq (${OPTS},D)
 FFLAGS += -g -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq (${OPTS},R)
 FFLAGS += -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq (${OPTS},X)
 FFLAGS += -g -O -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge -fbounds-check
endif
endif
endif
# endif LINUX

# SGI6, SGI64, SGIMP64
ifneq (,$(findstring SGI,${PVM_ARCH}))
#C++ := CC -LANG:std
#CC := cc
# 20000302: -w suppresses warnings which will swamp linker
C++ := g++ -w
CC := gcc
CPPFLAGS += ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := f90 -cpp
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf
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
# endif SGI64, SGIMP64
endif
ifeq (${CC},gcc)
 LDFLAGS := $(GCC_LDFLAGS_SZ_SPC) ${LDFLAGS}
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
 C++FLAGS := ${CFLAGS}
endif
# endif CC=gcc
ifeq (${CC},cc)
 LDFLAGS := ${SGI_ABI_FLG} ${LDFLAGS}
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
# end fortran flags
endif
# endif SGI6, SGI64, SGIMP64

ifeq (${UNAMES},SunOS) 
C++ := g++
CC := gcc
CFLAGS := -Wall
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} 
FC := f90 -DHIDE_SHR_MSG
#FFLAGS := -xs -stackvar -e -Qoption f90comp -r8const
FFLAGS := -xs -stackvar -e 
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf -lsunmath -lsocket -lnsl
LDFLAGS += -lF77 -lM77 -lresolv
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
 C++FLAGS := ${CFLAGS}
endif
# endif SunOS=SUN4SOL2,SUNMP

ifeq (${PVM_ARCH},WIN32)
NETCDF_INC := /usr/local/include
NETCDF_LIB := /usr/local/lib
C_ONLY := Y
C++ := g++
CC := gcc -ansi
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS := ${CPP_TKN} -Di386 ${CPP_PTH} -I${NETCDF_INC}
FC := g77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_c++ -lnco_c++ -L${NETCDF_LIB} -lnetcdf
LDFLAGS += 
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS := -O
endif
ifeq (${OPTS},D)
 CFLAGS := -g
endif
endif
# endif WIN32

# Link to Electric Fence library, see man efence
ifeq (${PVM_ARCH},LINUX)
ifneq (,$(findstring ${OPTS},DX))
# LDFLAGS += -lefence
endif
endif
# endif LINUX

# Define any remaining variables

# Default targets
all: slr_spc
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes that matter
.SUFFIXES: .cc .c .o .F .d
# Delete targets that were not successfully made
.DELETE_ON_ERROR:

# Targets in bin
SLR_SPC_OBJ=	${MY_OBJ_DIR}/slr_spc.o \
		${MY_OBJ_DIR}/getopt_bsd.o \
		${MY_OBJ_DIR}/flx_slr_frc.o
slr_spc: ${MY_BIN_DIR}/slr_spc
${MY_BIN_DIR}/slr_spc: lib $(SLR_SPC_OBJ)
	${C++} -o $@ $(SLR_SPC_OBJ) ${LDFLAGS} 
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
	@printf "C++ = ${C++}\n"
	@printf "C++FLAGS = ${C++FLAGS}\n"
	@printf "CC = ${CC}\n"
	@printf "CFLAGS = ${CFLAGS}\n"
	@printf "CPP = ${CPP}\n"
	@printf "CPPFLAGS = ${CPPFLAGS}\n"
	@printf "CPP_PTH = ${CPP_PTH}\n"
	@printf "CPP_TKN = ${CPP_TKN}\n"
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
	@printf "VRS_SNG = ${VRS_SNG}\n"
distclean : dst_cln
dst_cln : cln bld_cln
install : inc 
tags : 
	etags *.cc *.c *.com *.F *.hh *.h Makefile
test : tst
tst :
uninstall : dst_cln bin_cln inc_cln lib_cln
# Targets in dat
dat_cln :
	-cd ${MY_DAT_DIR}; rm -f out.nc
# Targets in dpn
dpn : ${MDL_DPN}
dpn_cln :
	rm -f ${MDL_DPN}
# Targets in inc
inc : $(MY_INC_TRG)
	cp -f ${MDL_INC_TRG} ${MY_INC_DIR}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
lib :
	-cd ${HOME}/c++; gmake lib
lib_cln :
	-cd ${MY_LIB_DIR};
# Targets in obj
obj_cln :
	rm -f ${MDL_OBJ}

include ${MDL_DPN}

