# $Id$

# Purpose: Makefile for index of refraction module idx_rfr

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

ifndef ${DBG}
# Debugging token N=No (default) Y=Yes
 DBG := N
endif # endif DBG
ifndef NETCDF_INC
 NETCDF_INC := /usr/local/include 
endif
ifndef NETCDF_LIB
 NETCDF_LIB := /usr/local/lib
endif
ifndef ${OMP} # OpenMP
 OMP := N
endif # endif OMP
ifndef OPTS
 OPTS := O
endif
ifndef ${PRC}
# Precision token D=Double (default) S=Single
 PRC := S
endif # endif PRC
ifndef PVM_ARCH
 PVM_ARCH := $(shell ${HOME}/sh/pvmgetarch)
endif
ifndef UNAMES
 UNAMES := $(shell uname -s)
endif

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
MDL_SRC := idx_rfr_H2O_ice_sgw.f idx_rfr_H2O_lqd_bpb.f
# Directories to search for source files
MDL_PTH := .
# Dependency list for executable
MDL_OBJ := $(addprefix ${MY_OBJ_DIR}/,$(addsuffix .o, $(basename ${MDL_SRC}))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix ${MY_DPN_DIR}/,$(addsuffix .d, $(basename ${MDL_SRC}))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.d file
VPATH := $(subst ${space},:,${MDL_PTH})
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,${MDL_PTH},-I${dir})

# Variables having to do with binary executables created by module
MDL_BIN_TRG := # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_DAT_TRG := idx_rfr_DKS91.nc idx_rfr_SAJ93.nc idx_rfr_Vol73.nc idx_rfr_shettle.nc idx_rfr_NNM98.nc idx_rfr_PaW75.nc idx_rfr_DHE02.nc idx_rfr_Ill_Kao_Mnt.nc idx_rfr_roush_illite.nc idx_rfr_roush_kaolinite.nc idx_rfr_limestone.nc idx_rfr_roush_montmorillonite.nc idx_rfr_shettle_Fe2O3_SiO2.nc idx_rfr_triaud_Fe2O3.nc
# Include targets
MDL_DAT_SYM_LNK := # Symbolic links
MDL_DAT_STB := $(MDL_DAT_TRG) $(MDL_DAT_SYM_LNK) # All targets in MY_DAT_DIR
MDL_DAT := $(addprefix ${MY_DAT_DIR}/,$(MDL_DAT_STB)) # dst_cln removes these files

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

# Rules for converting CDL files
${MY_DAT_DIR}/%.nc : %.cdl
	ncgen -b -o $@ $<

${MY_DAT_DIR}/idx_rfr_shettle_Fe2O3_SiO2.nc : idx_rfr_shettle_Fe2O3_SiO2.cdl
	ncgen -b -o $@ $<
	ncap -O -s 'idx_rfr_Fe2O3_avg_hitran96_rl=0.5*(idx_rfr_Fe2O3_oray_hitran96_rl+idx_rfr_Fe2O3_eray_hitran96_rl)' -s 'idx_rfr_Fe2O3_avg_hitran96_img=0.5*(idx_rfr_Fe2O3_oray_hitran96_img+idx_rfr_Fe2O3_eray_hitran96_img)' -s 'idx_rfr_SiO2_rl=0.5*(idx_rfr_SiO2_oray_hitran96_rl+idx_rfr_SiO2_eray_hitran96_rl)' -s 'idx_rfr_SiO2_img=0.5*(idx_rfr_SiO2_oray_hitran96_img+idx_rfr_SiO2_eray_hitran96_img)' ${DATA}/aca/idx_rfr_shettle_Fe2O3_SiO2.nc ${DATA}/aca/idx_rfr_shettle_Fe2O3_SiO2.nc

# Automatically generate a dependency file for each source file.
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.d
# Linux gcc may return an extra `.F' on Fortran names, e.g., `hello.F.o: hello.F'
# (.F)? gets rid of this extra `.F'
${MY_DPN_DIR}/%.d : %.F
# Following command makes, e.g., f.d begin "f.o f.d : f.F ..."
# Since f.o is not preceded by ${MY_OBJ_DIR}, objects are not recompiled when sources are touched.
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.F.o/$*.o $@/g;' > $@
# Following command makes, e.g., f.d begin "/home/zender/obj/LINUX/f.o f.d : f.F ..."
# This works fairly well, but is a hack.
# First pattern substitutes MY_OBJ_DIR_RX, which has placeholders for slashes.
# Second pattern substitutes slashes for the placeholders.
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*(.F)?.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@
# Following command makes, e.g., f.d begin "${MY_OBJ_DIR}/f.o f.d : f.F ..."
# This would be the ideal command but I can't get the dollar sign to show up.
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.F.o/\${dollar}MY_OBJ_DIR\/$*.o $@/g;' > $@

${MY_DPN_DIR}/%.d : %.f
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*(.F)?.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.c
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/$*.o $@/g;' > $@
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.cc
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/$*.o $@/g;' > $@
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

# First LDFLAGS is for typical C programs with netCDF, math, and networking
# Second LDFLAGS enables C/Fortran linking

ifeq (${PVM_ARCH},AIX)
ifeq (${OMP},Y)
 CC := xlc_r -qlanglvl=extc99
 FC := xlf90_r
 OMP_FLG := -qsmp=noauto 
else
 CC := xlc
 CPP_DFN += -U_OPENMP
 FC := f77
endif # endif OMP
C++ := xlC
CFLAGS := -qmaxmem=-1
CPP := /lib/cpp -P
# AIX xlf requires wrapping all CPP tokens of type -DFOO and -UFOO in -WF flag
CPP_DFN := -D${PVM_ARCH}
CPP_DFN_F := -WF,$(subst $(space),$(comma),$(CPP_DFN))
CPPINC := ${CPP_PTH} -I${NETCDF_INC}
FFLAGS := ${OMP_FLG} -qmaxmem=-1 -qarch=auto -qspillsize=2500
FIXEDFLAGS := -qsuffix=f=f:cpp=F -qfixed=132
FREEFLAGS := -qsuffix=f=f90:cpp=F90
LD := ld
# MATCH adds to LDFLAGS: -bmaxdata:2147483648 -bmaxstack:268435456
LDFLAGS := ${OMP_FLG} -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${PRC},D)
 FFLAGS += -qrealsize=8 -qdpc=e
else
 FFLAGS += -qrealsize=4
endif
ifeq (${OPTS},O)
 CFLAGS :=  -O2
 FFLAGS += -O3 -qstrict -Q
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=ov:zero:inv:en -qsigtrap -qextchk -C
endif
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif AIX

ifeq (${PVM_ARCH},ALPHA)
ifeq (${OMP},Y)
 OMP_FLG := -omp
endif # endif OMP
C++ := cxx -std strict_ansi
CC := cc
CFLAGS := ${OMP_FLG}
CPP := cpp
CPPFLAGS := -D${PVM_ARCH} -DOSF1 ${CPP_PTH} -I${NETCDF_INC}
FC := f90
FFLAGS := -c ${OMP_FLG} -automatic -fpe3
FIXEDFLAGS := -extend_source
FREEFLAGS := -DHIDE_SHR_MSG -free
LD := ld
LDFLAGS := ${OMP_FLG} -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
endif
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
 FFLAGS := -N 132
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS := -g -N 132
endif
ifeq (${OPTS},X)
 CFLAGS := -g -N 132
 FFLAGS := -g -N 132 -e i
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
 FC := g77
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

ifeq (${UNAMES},IRIX64)
C++ := CC
CC := cc
CFLAGS := -64 -mips4
CPP := /usr/lib/cpp
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90 -cpp
FFLAGS := -64 -mips4 -extend_source -mp -mpio
LD := ld
# 20000706: Added -mp -mpio for OpenMP compliance
LDFLAGS := -64 -mp -mips4 -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${PRC},D)
 FFLAGS += -r8 -i4 
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
 FFLAGS += -g -check_bounds -trapuv
endif
endif
# endif IRIX64=SGI64,SGIMP64

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
LDFLAGS := -L/opt/SUNWspro/lib -lf77compat -L${MY_LIB_DIR} -lcsz_f77 -L${NETCDF_LIB} -lnetcdf -lsunmath -lsocket -lnsl
LDFLAGS += -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
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
#FFLAGS += -r8const
#FFLAGS += -r8 -i4
# Using xtypemap appears to cause aer to fail
 FFLAGS += -xtypemap=real:64,double:64,integer:32
ifeq (${FC},f77)
 FFLAGS += -r8 -i4
endif
endif
endif
# endif SunOS=SUN4SOL2,SUNMP

# Define any remaining variables

# Default targets
all: ${MY_OBJ_DIR}/idx_rfr_H2O_ice_sgw.o ${MY_OBJ_DIR}/idx_rfr_H2O_lqd_bpb.o
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug
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
dst_cln : cln bld_cln
install : inc 
tags : 
	etags *.cc *.c *.com *.F *.hh *.h Makefile
uninstall : dst_cln bin_cln inc_cln lib_cln
# Targets in dat
data:	dat
dat : $(MDL_DAT)
dat_cln : 
	rm -f $(MDL_DAT)
# Targets in dpn
depend:	dpn
dpn : ${MDL_DPN}
dpn_cln :
	rm -f ${MDL_DPN}
# Targets in inc
include: inc
inc : ${MDL_INC}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
library : lib
lib :
	-cd ${HOME}/c++; gmake lib
lib_cln :
	-cd ${MY_LIB_DIR};
# Targets in obj
object : obj
obj : ${MDL_OBJ}
obj_cln :
	rm -f ${MDL_OBJ}

# include ${MDL_DPN}
