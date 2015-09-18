# $Id$

# Purpose: Makefile for time module time

MY_BLD_DIR := .
MY_DAT_DIR := .
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
ifndef PVM_ARCH
 PVM_ARCH := $(shell ${MY_BLD_DIR}/pvmgetarch)
endif
ifndef UNAMES
 UNAMES := $(shell uname -s)
endif
ifndef NETCDF_INC
 NETCDF_INC := /usr/local/include 
endif
ifndef NETCDF_LIB
 NETCDF_LIB := /usr/local/lib
endif

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := ${null} ${null}
comma := ,
newline := \n

# Source file names with directories removed
MDL_SRC := slr_crd_jjm.F slr_crd_utl.F time_diff.c time_gmt2unix.c time_tst.c time_unix2gmt.c
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
MDL_BIN_TRG := time_gmt2unix time_unix2gmt slr_crd_jjm # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All C++ files in MY_BIN_DIR
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
	${C++} ${CPPFLAGS} ${C++FLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)

# Fortran rules
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
%.mod : %.F90
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
DPN_GNR := makdep
${MY_DPN_DIR}/%.d : %.c
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.d : %.cc
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.d : %.F
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.d : %.F90
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@

# Rules for installing header files
${MY_INC_DIR}/%.com : %.com
	cp -f -p $(notdir $@) $@
${MY_INC_DIR}/%.F90 : %.F90
	cp -f -p $(notdir $@) $@

# Automatically generate a dependency file for each source file.
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.d

# The first LDFLAGS is for typical C programs with netCDF, math, and networking
# The second LDFLAGS enables C/Fortran linking

ifeq (${PVM_ARCH},AIX)
C++ := xlC
CC := xlc_r
CPP := /lib/cpp -P
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := xlf90_r
FFLAGS := -qarch=pwr3 -qrealsize=8 -O3 -qstrict -qdpc=e -qsmp=noauto -qspillsize=2500
FIXEDFLAGS := -qfixed
FREEFLAGS := -qsuffix=f=f90:cpp=F90
LD := ld
LDFLAGS := -qsmp=noauto -bmaxdata:0x80000000 -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lncaru -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS :=  -O2
 FFLAGS += -O3 -qstrict
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS += -g
endif
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif AIX

ifeq (${PVM_ARCH},ALPHA)
C++ := cxx
CC := cc
CPP := cpp
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
FFLAGS := -r8 -i4 -c -omp -automatic
FIXEDFLAGS := -extend_source -omp -automatic
FREEFLAGS := -DHIDE_SHR_MSG -free
LD := ld
LDFLAGS := -omp -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O2
 FFLAGS += -O3 -strict
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS += -g
endif
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif ALPHA

ifeq (${PVM_ARCH},CRAY)
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf
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
# 1997/10/21 Added -P to suppress #line # directives on Fortran files
# 1998/09/07 Added -N to prevent CPP from producing .f that complain about single quotes ('s)
# NB: -N and -P are not part of CPPFLAGS because they should only used for Fortran files
${MY_OBJ_DIR}/%.o : %.F
	${CPP} -N -P ${CPPFLAGS} $< > ${MY_OBJ_DIR}/$(patsubst %.F,%.f,$(notdir $<))
	${FC} ${FFLAGS} -b $@ ${MY_OBJ_DIR}/$(patsubst %.F,%.f,$(notdir $<))
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} $<
	mv -f $(notdir $@) ${MY_OBJ_DIR}
endif
# endif CRAY

# Works on LINUX and LINUXALPHA
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
C++ := g++
CC := gcc
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
ifeq (${PVM_ARCH},LINUX)
 FC := pgf90
else # LINUXALPHA
 FC := g77
endif
LD := ld
LDFLAGS := -mp -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf 
LEX := flex
LINT := lint
YACC := bison
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
 C++FLAGS := ${CFLAGS}
endif
ifeq (${FC},pgf90)
 FFLAGS := -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign -Ktrap=fp
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
endif
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
 FFLAGS += -g -O -fdebug-kludge -fbounds-check
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
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90 -cpp
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf
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

ifeq (${PVM_ARCH},SUN4)
C++ := g++
CC := acc
CPP := cpp
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf
LDFLAGS += -cg92 -L/opt/SUNWspro/SC3.0/lib/cg92 -lF77 -lM77 -lresolv
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
endif
# endif SUN4

ifeq (${UNAMES},SunOS) 
C++ := g++
CC := gcc -ansi
CFLAGS := -Wall
CPP := cpp
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC} 
FC := f90 -DHIDE_SHR_MSG
#FFLAGS := -xs -stackvar -e -Qoption f90comp -r8const
FFLAGS := -xs -stackvar -e 
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf -lsunmath -lsocket -lnsl -lm
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
C_ONLY := Y
C++ := g++
CC := gcc -ansi
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS := -D${PVM_ARCH} -Di386 ${CPP_PTH} -I${NETCDF_INC}
FC := g77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf
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
ifneq (${null},$(findstring ${OPTS},DX))
 LDFLAGS += -lefence
endif
endif
# endif LINUX

# Define any remaining variables
libcsz_c := ${MY_LIB_DIR}/libcsz_c

# Default targets
all: bin_dir obj_dir lib_dir ${MDL_BIN_TRG}
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .cc .c .o .F .d .F90
# Delete targets which were not successfully made
.DELETE_ON_ERROR:
lib_dir :
	- if ! test -d ${MY_LIB_DIR}; then mkdir ${MY_LIB_DIR}; fi
bin_dir :
	- if ! test -d ${MY_BIN_DIR}; then mkdir ${MY_BIN_DIR}; fi
obj_dir :
	- if ! test -d ${MY_OBJ_DIR}; then mkdir ${MY_OBJ_DIR}; fi

# Targets in bin
time_unix2gmt:	${MY_BIN_DIR}/time_unix2gmt
${MY_BIN_DIR}/time_unix2gmt:	${MY_OBJ_DIR}/time_unix2gmt.o
	${CC} -o $@ ${MY_OBJ_DIR}/time_unix2gmt.o ${LDFLAGS}
	chmod 755 $@
time_gmt2unix:	${MY_BIN_DIR}/time_gmt2unix
${MY_BIN_DIR}/time_gmt2unix:	${MY_OBJ_DIR}/time_gmt2unix.o
	${CC} -o $@ ${MY_OBJ_DIR}/time_gmt2unix.o ${LDFLAGS}
	chmod 755 $@
time_tst:	${MY_BIN_DIR}/time_tst
${MY_BIN_DIR}/time_tst:	${MY_OBJ_DIR}/time_tst.o
	${CC} -o $@ ${MY_OBJ_DIR}/time_tst.o ${LDFLAGS}
	chmod 755 $@
SLR_CRD_JJM_OBJ := 	${MY_OBJ_DIR}/slr_crd_jjm.o \
			${MY_OBJ_DIR}/slr_crd_utl.o
slr_crd_jjm: ${MY_BIN_DIR}/slr_crd_jjm
${MY_BIN_DIR}/slr_crd_jjm: $(SLR_CRD_JJM_OBJ) 
	${FC} -o $@ $(SLR_CRD_JJM_OBJ) ${LDFLAGS}
	chmod 755 $@

bin : ${MDL_BIN_TRG}
	rm -f ${MDL_BIN}
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
	@printf "FC = ${FC}\n"
	@printf "GSL_LIB = ${GSL_LIB}\n"
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
	@printf "OPTS = ${OPTS}\n"
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
	cp -f ${MDL_INC_TRG} ${MY_INC_DIR}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
# Targets in obj
obj : ${MDL_OBJ}
obj_cln :
	rm -f ${MDL_OBJ}

# Create dependency files only if they will not be immediately deleted
INCLUDE_DPN := TRUE
GOALS_WHICH_DELETE_DEPENDENCY_FILES := cln clean dir distclean dst_cln dpn_cln tags uninstall
ifeq (,$(findstring $(MAKECMDGOALS),${GOALS_WHICH_DELETE_DEPENDENCY_FILES}))
 INCLUDE_DPN := TRUE
else
 INCLUDE_DPN := FALSE
endif
ifeq (${INCLUDE_DPN},TRUE)
-include ${MDL_DPN}
endif
