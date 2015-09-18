# $Id$

# Purpose: Makefile for HDF module hdf

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
ifndef HDF_INC
 HDF_INC := /usr/local/hdf5/include 
endif
ifndef HDF_LIB
 HDF_LIB := -L/usr/local/hdf5/lib -lhdf5 -ldf -lz
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
MDL_SRC := hdf.c
# Directories to search for source files
MDL_PTH := . ${HOME}/include /contrib/include /usr/include/gsl
# Dependency list for executable
MDL_OBJ := $(addprefix ${MY_OBJ_DIR}/,$(addsuffix .o, $(basename ${MDL_SRC}))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix ${MY_DPN_DIR}/,$(addsuffix .d, $(basename ${MDL_SRC}))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.d file
VPATH := $(subst ${space},:,${MDL_PTH})
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,${MDL_PTH},-I${dir})

# Variables having to do with binary executables created by module
MDL_BIN_TRG := hdf # Binary targets
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

# Default Fortran pattern rules: CRAY and AIX must override these rules
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<

# Automatically generate a dependency file for each source file.
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.d
# Linux gcc may return an extra `.F' on Fortran names, e.g., `hello.F.o: hello.F'
# (.F)? gets rid of this extran `.F'
${MY_DPN_DIR}/%.d : %.F
# Following command makes, e.g., f.d begin "f.o f.d : f.F ..."
# Since f.o is not preceded by ${MY_OBJ_DIR}, objects are not recompiled when sources are touched.
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.F\.o/$*.o $@/g;' > $@
# Following command makes, e.g., f.d begin "/home/zender/obj/LINUX/f.o f.d : f.F ..."
# This works fairly well, but is a hack
# First pattern substitutes MY_OBJ_DIR_RX, which has placeholders for slashes
# Second pattern substitutes slashes for the placeholders
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*(\.F)?\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@
# Following command makes, e.g., f.d begin "${MY_OBJ_DIR}/f.o f.d : f.F ..."
# This would be the ideal command but I can't get the dollar sign to show up.
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.F\.o/\${dollar}MY_OBJ_DIR\/$*.o $@/g;' > $@

${MY_DPN_DIR}/%.d : %.c
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.cc
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

# The first LDFLAGS is for typical C programs with hdf, math, and networking
# The second LDFLAGS enables C/Fortran linking
ifeq (${PVM_ARCH},CRAY)
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I$(HDF_INC)
FC := f90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lnco $(HDF_LIB)
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

# Works on LINUX and LINUXALPHA
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
C++ := g++
# 19980517 As of gcc-2.8.0 using -ansi switch breaks because resolv.h structure have changed
#CC := gcc -ansi
CC := gcc -Wall
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS := -D${PVM_ARCH} -Di386 ${CPP_PTH} -I$(HDF_INC)
FC := pgf90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lnco $(HDF_LIB) 
LDFLAGS += 
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS := -O
 FFLAGS := -O -Mextend
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS := -g -Mextend
endif
endif
# endif LINUX

ifeq (${PVM_ARCH},AIX)
#C++ := xlC
C++ := g++
CC := gcc -ansi
CPP := gcc
# NB: Hack to make AIX look like LINUX (g77 + gcc + g++)
CPPFLAGS := -DLINUX ${CPP_PTH} -I$(HDF_INC)
#FC := xlf
FC := g77
FFLAGS := -ffixed-line-length-132
#FFLAGS := -qfixed=132
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lnco $(HDF_LIB)
#LDFLAGS += -lxlf90 -lxlf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O
 CPP := ${CPP} ${CPPFLAGS}
 PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS += -O
# FFLAGS += -O -NS2000 -qfixed=132
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 CPP := ${CPP} ${CPPFLAGS}
 PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS += -g
# FFLAGS := -g -NS2000 -qfixed=132
endif
endif
# endif AIX

# SGI6, SGI64, SGIMP64
ifneq (,$(findstring SGI,${PVM_ARCH}))
#C++ := CC -LANG:std
#CC := cc
# 20000302: -w suppresses warnings which will swamp linker
C++ := g++ -w
CC := gcc
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I$(HDF_INC)
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
# endif SGI64, SGIMP64
endif
ifeq (${CC},gcc)
 LDFLAGS += $(GCC_LDFLAGS_SZ_SPC) -L${MY_LIB_DIR} -lcsz_c++ $(HDF_LIB)
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
 LDFLAGS += ${SGI_ABI_FLG} -L${MY_LIB_DIR} -lcsz_c++ $(HDF_LIB)
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
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I$(HDF_INC)
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lnco $(HDF_LIB)
LDFLAGS += -cg92 -L/opt/SUNWspro/SC3.0/lib/cg92 -lF77 -lM77 -lresolv -lm
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

ifeq (${PVM_ARCH},SUN4SOL2)
C++ := g++
CC := gcc -Wall
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I$(HDF_INC)
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lnco $(HDF_LIB) -lsunmath -lsocket -lnsl
LDFLAGS += -lF77 -lM77 -lresolv -lm
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
# NB: 98/06/01 -C (range-checking) is not supported by Sun f90
ifeq (${FC},f77)
 FFLAGS += ${FFLAGS} -C
endif
endif
endif
# endif SUN4SOL2

ifeq (${PVM_ARCH},SUNMP)
C++ := g++
#CC := gcc -ansi -pedantic
CC := gcc -Wall
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I$(HDF_INC)
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lnco $(HDF_LIB) -lsunmath -lthread -lsocket -lnsl
LDFLAGS += -lF77 -lM77 -lresolv -lm
LEX := flex
LINT := lint
YACC := bison
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
# NB: 98/06/01 -C (range-checking) is not supported by Sun f90
ifeq (${FC},f77)
 FFLAGS += ${FFLAGS} -C
endif
endif
endif
# endif SUNMP

ifeq (${PVM_ARCH},WIN32)
C_ONLY := Y
C++ := g++
CC := gcc -ansi
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS := -D${PVM_ARCH} -Di386 ${CPP_PTH} -I$(HDF_INC)
FC := g77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lnco $(HDF_LIB)
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
# Electric fence only works on malloc(), so is no good on most C++ programs
ifneq (${null},$(findstring LINUX,${PVM_ARCH}))
ifneq (${null},$(findstring ${OPTS},DX))
 LDFLAGS += -lefence
endif
endif
# endif LINUX or LINUXALPHA

# Define any remaining variables
libcsz_c := ${MY_LIB_DIR}/libcsz_c

# Default targets
all: ${MDL_BIN_TRG}
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes which matter
.SUFFIXES: .cc .c .o .F .d
# Delete targets which were not successfully made
.DELETE_ON_ERROR:
lib_dir :
	- if ! test -d ${MY_LIB_DIR}; then mkdir ${MY_LIB_DIR}; fi
bin_dir :
	- if ! test -d ${MY_BIN_DIR}; then mkdir ${MY_BIN_DIR}; fi
obj_dir :
	- if ! test -d ${MY_OBJ_DIR}; then mkdir ${MY_OBJ_DIR}; fi

# Targets in bin
HDF_OBJ := ${MY_OBJ_DIR}/hdf.o 
ifeq (${null},$(findstring LINUX,${PVM_ARCH}))
 HDF_OBJ += ${MY_OBJ_DIR}/alloca.o ${MY_OBJ_DIR}/getopt.o ${MY_OBJ_DIR}/getopt1.o
endif
hdf:	${MY_BIN_DIR}/hdf
${MY_BIN_DIR}/hdf:	$(HDF_OBJ) 
	${CC} -o $@ $(HDF_OBJ) ${LDFLAGS}
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
	@printf "FC = ${FC}\n"
	@printf "FFLAGS = ${FFLAGS}\n"
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
	@printf "HDF_INC = $(HDF_INC)\n"
	@printf "HDF_LIB = $(HDF_LIB)\n"
	@printf "OMP = ${OMP}\n"
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
	cp -f -p ${MDL_INC_TRG} ${MY_INC_DIR}
inc_cln :
	rm -f ${MDL_INC}
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
