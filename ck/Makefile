# $Id$

# Purpose: Makefile for Fortran module ck

# C-Fortran linking issues:
# pgf90 subroutines have one underscore by default
# pgf90 underscore behavior is altered by -Mnosecond_underscore
# g77 subroutines have two underscores by default
# g77 functions (e.g., newdate()) have one underscore by default 
# g77 underscore behavior is altered by -fno-second-underscore 

# Usage: 
# make

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

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := ${null} ${null}
comma := ,
newline := \n

# Source file names with directories removed
MDL_SRC := ck.F htrn.c
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
MDL_BIN_TRG := htrn # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := hitran.com hitran.h # Include targets
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
${MY_INC_DIR}/%.hh : %.hh
	cp -f -p $(notdir $@) $@
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
CC := xlc_r -qlanglvl=extc99
CPP := /lib/cpp -P
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := xlf90_r
FFLAGS := -qarch=pwr3 -qrealsize=8 -O3 -qstrict -qdpc=e -qsmp=noauto -qspillsize=2500
FIXEDFLAGS := -qfixed
FREEFLAGS := -qsuffix=f=f90:cpp=F90
LD := ld
LDFLAGS := -qsmp=noauto -bmaxdata:0x80000000 -L${MY_LIB_DIR} -lcsz_f90 -L${NETCDF_LIB} -lncaru -lnetcdf
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
LDFLAGS := -omp -L${MY_LIB_DIR} -lcsz_f90 -L${NETCDF_LIB} -lnetcdf
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
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f90 -L${NETCDF_LIB} -lnetcdf
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

# Works on LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, and FREEBSD and MACOSX
ifneq (${null},$(findstring ${PVM_ARCH},LINUXALPHALINUXAMD64LINUXARMFREEBSDMACOSX))
C++ := g++
CC := gcc
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
ifeq (${PVM_ARCH},LINUX)
 FC := pgf90
else # LINUXALPHA
 FC := g77
endif
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f90 -L${NETCDF_LIB} -lnetcdf
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

ifeq (${UNAMES},IRIX64)
C++ := CC
CC := cc
CFLAGS := -64 -mips4
CPP := cpp
# spoof SGI64
CPPFLAGS := -DSGI64 ${CPP_PTH} -I${NETCDF_INC}
FC := f90 -cpp -DHIDE_SHARE_MSG
FFLAGS := -64 -mips4 -extend_source -mp -mpio
LD := ld
# 20000706: Added -mp -mpio for OpenMP compliance
LDFLAGS := -64 -mp -mips4 -L${MY_LIB_DIR} -lcsz_f90 -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
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
 FFLAGS += -g -check_bounds -trapuv
endif
endif
# endif IRIX64=SGI64,SGIMP64

ifeq (${PVM_ARCH},SUN4)
C++ := g++
CC := acc
CPP := cpp
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f90 -L${NETCDF_LIB} -lnetcdf
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
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f90 -L${NETCDF_LIB} -lnetcdf -lsunmath -lsocket -lnsl
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
LDFLAGS := -L${MY_LIB_DIR} -lcsz_f90 -L${NETCDF_LIB} -lnetcdf
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

# Default targets
all: ${MDL_BIN_TRG} lib
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .cc .c .o .F .F90 .d .mod
# Delete targets which were not successfully made
.DELETE_ON_ERROR:
lib_dir :
	- if ! test -d ${MY_LIB_DIR}; then mkdir ${MY_LIB_DIR}; fi
bin_dir :
	- if ! test -d ${MY_BIN_DIR}; then mkdir ${MY_BIN_DIR}; fi
obj_dir :
	- if ! test -d ${MY_OBJ_DIR}; then mkdir ${MY_OBJ_DIR}; fi

# Targets in bin
CK_OBJ := 	${MY_OBJ_DIR}/ck.o ${MY_OBJ_DIR}/disort.o ${MY_OBJ_DIR}/d1mach.o ${MY_OBJ_DIR}/r1mach.o	${MY_OBJ_DIR}/linpak.o
ck:	${MY_BIN_DIR}/ck
${MY_BIN_DIR}/ck:	$(CK_OBJ) 
	${FC} -o $@ $(CK_OBJ) -lcsz_f90 ${LDFLAGS}
	chmod 755 $@

LBL_OBJ := 	${MY_OBJ_DIR}/lbl.o ${MY_OBJ_DIR}/rt_utl.o ${MY_OBJ_DIR}/hitran.o ${MY_OBJ_DIR}/disort.o ${MY_OBJ_DIR}/d1mach.o ${MY_OBJ_DIR}/r1mach.o	${MY_OBJ_DIR}/linpak.o
LBL_SRC_LITTLE_F := lbl.f 
lbl:	${MY_BIN_DIR}/lbl
${MY_BIN_DIR}/lbl:	$(LBL_OBJ) 
	${FC} -o $@ $(LBL_OBJ) -lcsz_f90 ${LDFLAGS}
	/bin/rm -f $(LBL_SRC_LITTLE_F)
	chmod 755 $@

HTRN_OBJ := 	${MY_OBJ_DIR}/htrn.o ${MY_OBJ_DIR}/getopt_bsd.o
htrn:	${MY_BIN_DIR}/htrn
${MY_BIN_DIR}/htrn:	$(HTRN_OBJ) 
	${CC} -o $@ $(HTRN_OBJ) ${LDFLAGS}
	chmod 755 $@
bin : ${MDL_BIN_TRG}
	rm -f ${MDL_BIN}
bin_cln: 
	rm -f ${MDL_BIN}
# Targets in bld
bld_cln :
	cd ${MY_BLD_DIR}; rm -f TAGS
clean : cln
cln : dat_cln dpn_cln obj_cln mod_cln
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
dst_cln : cln bld_cln lib_cln
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
inc : ${MDL_INC}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
lib : inc 
# Targets in mod
mod : ${MDL_OBJ}
	mv -f $(wildcard ${MY_OBJ_DIR}/*.mod) .
mod_cln :
	mv -f $(wildcard *.mod) ${MY_OBJ_DIR}
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
# Following incorporates dependency files into Makefile rules
#-include ${MDL_DPN}
endif

${MY_OBJ_DIR}/rt_utl.o:	${HOME}/aca/rt_utl.F
${MY_OBJ_DIR}/hitran.o:	${HOME}/aca/hitran.F
${MY_OBJ_DIR}/htrn_c++.o: ${HOME}/c++/htrn_c++.cc
