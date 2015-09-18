# $Id$

# Purpose: Makefile for IDL module idl

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
ifdef OMP
 CPPFLAGS := -DOMP
 CFLAGS := -mp -mpio
 LDFLAGS := -mp
endif
ifndef PVM_ARCH
 PVM_ARCH := $(shell ${MY_BLD_DIR}/pvmgetarch)
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
# We replace `/' by `cszzsc', call perl, then replace `cszzxc' by `/'
# Unique character(s) to substitute for ${slash} before passing to perl regex
slash_rx := cszzsc
# Unique character(s) to replace by ${slash_rx} before passing to perl regex
slash := /
MY_OBJ_DIR_RX := $(subst ${slash},${slash_rx},${MY_OBJ_DIR})
MY_DPN_DIR_RX := $(subst ${slash},${slash_rx},${MY_DPN_DIR})

# Source file names with directories removed
MDL_SRC := ccc.cc cls.cc vec.cc utl.cc xtr.cc tdy.cc flx_sfc.cc blm.cc cln.cc gsl.cc
# Directories to search for source files
MDL_PTH := . ${HOME}/include /usr/local/include/gsl
# Dependency list for executable
MDL_OBJ := $(addprefix ${MY_OBJ_DIR}/,$(addsuffix .o, $(basename ${MDL_SRC}))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix ${MY_DPN_DIR}/,$(addsuffix .d, $(basename ${MDL_SRC}))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.d file
VPATH := $(subst ${space},:,${MDL_PTH})
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,${MDL_PTH},-I${dir})

# Variables having to do with binary executables created by module
MDL_BIN_TRG := ccc gsl # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := libcsz_c++.hh utl.hh xtr.hh vec.hh dbg.hh phys_cst.hh libcsm_c++.hh lsm.hh flx_sfc.hh blm.hh tdy.hh chm.hh # Include targets
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

# Rules for installing header files
${MY_INC_DIR}/%.hh : %.hh
	cp -f -p $(notdir $@) $@

# Automatically generate a dependency file for each source file
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
# This would be the ideal command but I can't get the dollar sign to show up
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.F\.o/\${dollar}MY_OBJ_DIR\/$*.o $@/g;' > $@

${MY_DPN_DIR}/%.d : %.f
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.c
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.cc
# NB: Use ${C++} rather than ${CPP} on C++ files for now because, e.g., SUNMP cpp does not 
# preprocess .cc files quite correctly
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/$*.o $@/g;' > $@
	${C++} -M ${CPPFLAGS} $< | perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

# The first LDFLAGS is for typical C programs with netCDF, math, and networking
# The second LDFLAGS enables C/Fortran linking
ifeq (${PVM_ARCH},CRAY)
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf
LDFLAGS += -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O
 FFLAGS = -N 132
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS = -g -N 132
endif
ifeq (${OPTS},X)
 CFLAGS += -g -N 132
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
C++ := g++
# 1998/05/17 As of gcc-2.8.0 using -ansi switch breaks because resolv.h structure have changed
#CC := gcc -ansi
CC := gcc
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := pgf90
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq (${CC},gcc)
ifeq (${OPTS},O)
 CFLAGS += -O -Wall
endif
ifeq (${OPTS},D)
 CFLAGS += -g -Wall
endif
ifeq (${OPTS},R)
 CFLAGS += -Wall
endif
ifeq (${OPTS},X)
 CFLAGS += -g -O -Wall
endif
 C++FLAGS := ${CFLAGS}
endif
ifeq (${FC},pgf90)
ifeq (${OPTS},O)
 FFLAGS := -fast -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign
endif
ifeq (${OPTS},D)
 FFLAGS := -g -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign
endif
ifeq (${OPTS},R)
 FFLAGS := -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign
endif
ifeq (${OPTS},X)
 FFLAGS := -g -Mbounds -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign
endif
endif
ifeq (${FC},g77)
ifeq (${OPTS},O)
 FFLAGS := -O -ffixed-line-length-132 -fno-second-underscore
endif
ifeq (${OPTS},D)
 FFLAGS := -g -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq (${OPTS},R)
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq (${OPTS},X)
 FFLAGS := -g -O -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge -fbounds-check
endif
endif
endif
# endif LINUX

ifeq (${PVM_ARCH},RS6K)
C++ := g++
CC := gcc -ansi
CPP := /lib/cpp -P
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := xlf
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lncaru -lnetcdf_c++ -lnetcdf
LDFLAGS += -lxlf90 -lxlf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O2
 CPP := ${CPP} ${CPPFLAGS}
 PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS := -O -NS2000 -qfixed=132
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 CPP := ${CPP} ${CPPFLAGS}
 PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS := -g -NS2000 -qfixed=132
endif
${MY_OBJ_DIR}/%.o : %.F
	${CPP} ${CPPFLAGS} $< ${MY_OBJ_DIR}/$(basename $<).f 
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) ${MY_OBJ_DIR}/$(basename $<).f
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif RS6K

ifeq (${PVM_ARCH},AIX)
#C++ := xlC
C++ := g++
CC := gcc -ansi
CPP := gcc
# NB: Hack to make AIX look like LINUX (g77 + gcc + g++)
CPPFLAGS += -DSUNMP ${CPP_PTH} -I${NETCDF_INC}
#FC := xlf
FC := g77
FFLAGS := -ffixed-line-length-132
#FFLAGS := -qfixed=132
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf
#LDFLAGS += -lxlf90 -lxlf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O
 CPP := ${CPP} ${CPPFLAGS}
 PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS += -O
# FFLAGS += -O -NS2000 -qfixed=132
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 CPP := ${CPP} ${CPPFLAGS}
 PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS += -g
# FFLAGS := -g -NS2000 -qfixed=132
endif
endif
# endif AIX

ifeq (${PVM_ARCH},SGI5)
C++ := g++
CC := gcc -ansi
CPP := cpp
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf
LDFLAGS += -lF77 -lI77 -lU77 -lftn
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS := -O2 -e -Nl99 -extend_source
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS := -g -e -Nl99 -extend_source
endif
endif
# endif SGI5

ifeq (${PVM_ARCH},SGI64)
C++ := g++
CC := cc
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS += -64 -mips4 -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS += -64 -mips4 -O2
 FFLAGS := -64 -mips4 -O2 -extend_source
endif
ifeq (${OPTS},R)
 CFLAGS += -64 -mips4
 FFLAGS := -64 -mips4 -extend_source
endif
ifeq (${OPTS},D)
 CFLAGS += -64 -mips4 -g
 FFLAGS := -64 -mips4 -g -extend_source
endif
ifeq (${OPTS},X)
 CFLAGS += -64 -mips4 -g -trapuv
 FFLAGS := -64 -mips4 -g -extend_source -check_bounds -trapuv
endif
endif
# endif SGI64

ifeq (${PVM_ARCH},SGIMP64)
#C++ := CC
C++ := /usr/local/egcs/bin/g++
CC := gcc
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
LD := ld
#LDFLAGS += -64 -mips4 -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++
LDFLAGS += -mabi=64 -mips4 -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
# CFLAGS += -64 -mips4 -O2
 C++FLAGS := -O2 -Wall -mabi=64
 CFLAGS += -O2 -Wall -mabi=64
 FFLAGS := -64 -mips4 -O2 -extend_source
endif
ifeq (${OPTS},R)
# CFLAGS += -64 -mips4
 C++FLAGS := -Wall -mabi=64
 CFLAGS += -Wall -mabi=64
 FFLAGS := -64 -mips4 -extend_source
endif
ifeq (${OPTS},D)
# CFLAGS += -64 -mips4 -g
 C++FLAGS := -g -Wall -mabi=64
 CFLAGS += -g -Wall -mabi=64
 FFLAGS := -64 -mips4 -g -extend_source
endif
ifeq (${OPTS},X)
# CFLAGS += -64 -mips4 -g -trapuv
 C++FLAGS := -g -O -Wall -mabi=64
 CFLAGS += -g -O -Wall -mabi=64
 FFLAGS := -64 -mips4 -g -extend_source -check_bounds -trapuv
endif
endif
# endif SGIMP64

ifeq (${PVM_ARCH},SUN4)
C++ := g++
CC := acc
CPP := cpp
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf
LDFLAGS += -cg92 -L/opt/SUNWspro/SC3.0/lib/cg92 -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS := -fast -e -Nl99
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS := -g -e -Nl99
endif
endif
# endif SUN4

ifeq (${PVM_ARCH},SUN4SOL2)
C++ := g++
CC := gcc -ansi
CPP := cpp
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf -lsunmath -lsocket -lnsl
LDFLAGS += -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O2
 FFLAGS := -O -e
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS := -g -e
endif
ifeq (${OPTS},X)
 CFLAGS += -g 
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
#C++ := /contrib/egcs-1.1.2/bin/g++
#CC := gcc -ansi -pedantic
CC := gcc 
CPP := cpp
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf -lsunmath -lthread -lsocket -lnsl
LDFLAGS += -lF77 -lM77 -lresolv
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 C++FLAGS := -O2 -Wall
 CFLAGS += -O2 -Wall
 FFLAGS := -O -e
endif
ifeq (${OPTS},D)
 C++FLAGS := -g -Wall
 CFLAGS += -g -Wall
 FFLAGS := -g -e
endif
ifeq (${OPTS},X)
 C++FLAGS := -g -O -Wall
 CFLAGS += -g -O -Wall
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
CPPFLAGS += -D${PVM_ARCH} -Di386 ${CPP_PTH} -I${NETCDF_INC}
FC := g77
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lnco -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf
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
# endif WIN32

ifdef OMP
 C++ := pgCC
 CC := pgcc
endif
# endif OMP

# Define any remaining variables
libcsz_c++ := ${MY_LIB_DIR}/libcsz_c++
libcsm_c++ := ${MY_LIB_DIR}/libcsm_c++

# Default targets
all: ${MDL_BIN_TRG} lib 
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug 
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
CCC_OBJ := ${MY_OBJ_DIR}/ccc.o ${MY_OBJ_DIR}/cls.o
ifneq (${PVM_ARCH},LINUX)
 CCC_OBJ += ${MY_OBJ_DIR}/getopt.o ${MY_OBJ_DIR}/getopt1.o
endif
ccc:	lib $(CCC_OBJ)
	${C++} -o ${MY_BIN_DIR}/$@ $(CCC_OBJ) ${LDFLAGS} -lcsm_c++
	chmod 755 ${MY_BIN_DIR}/$@
GSL_OBJ := ${MY_OBJ_DIR}/gsl.o ${MY_OBJ_DIR}/cls.o
ifneq (${PVM_ARCH},LINUX)
 GSL_OBJ += ${MY_OBJ_DIR}/getopt.o ${MY_OBJ_DIR}/getopt1.o
endif
gsl:	lib $(GSL_OBJ)
	${C++} -o ${MY_BIN_DIR}/$@ $(GSL_OBJ) ${LDFLAGS} -L/usr/local/lib/gsl -lgslspecfunc -lgslerr
	chmod 755 ${MY_BIN_DIR}/$@
bin : ${MDL_BIN_TRG}
bin_cln: 
	rm -f ${MDL_BIN}
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
lib : inc $(libcsz_c++).a $(libcsm_c++).a
$(libcsz_c++).a : $(libcsz_c++).a(${MY_OBJ_DIR}/utl.o) \
	$(libcsz_c++).a(${MY_OBJ_DIR}/vec.o) \
	$(libcsz_c++).a(${MY_OBJ_DIR}/xtr.o)
$(libcsm_c++).a : $(libcsm_c++).a(${MY_OBJ_DIR}/blm.o) \
	$(libcsm_c++).a(${MY_OBJ_DIR}/cln.o) \
	$(libcsm_c++).a(${MY_OBJ_DIR}/chm.o) \
	$(libcsm_c++).a(${MY_OBJ_DIR}/flx_sfc.o) \
	$(libcsm_c++).a(${MY_OBJ_DIR}/tdy.o)
lib_cln :
	rm -f $(libcsz_c++).a
	rm -f $(libcsm_c++).a
# Targets in obj
obj : ${MDL_OBJ}
obj_cln :
	rm -f ${MDL_OBJ}

# Housekeeping
clean : cln
cln : dat_cln dpn_cln obj_cln
debug: dbg
dbg:
	printf "MDL_SRC = ${MDL_SRC}\n"
	printf "MDL_OBJ = ${MDL_OBJ}\n"
	printf "MDL_DPN = ${MDL_DPN}\n"
	printf "VPATH = ${VPATH}\n"
	printf "OBJ_PAT = $(OBJ_PAT)\n"
	printf "MY_OBJ_DIR_RX = ${MY_OBJ_DIR_RX}\n"
	printf "MY_DPN_DIR_RX = ${MY_DPN_DIR_RX}\n"
distclean : dst_cln
dst_cln : cln bld_cln lib_cln
install : inc 
tags : 
	etags *.pro *.com *.sh *.txt *.tbl Makefile ../aca/*.pro ../arese/*.pro ../aeroce/*.pro ../mie/*.pro ../avhrr/*.pro ../dot/idl ../dot/idle
test : tst
tst :
	ccc --dbg=5
uninstall : dst_cln bin_cln inc_cln lib_cln

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


