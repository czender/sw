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
MDL_SRC := getopt.c getopt1.c alloca.c getopt_test.c
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
MDL_BIN_TRG := getopt_test # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All C++ files in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := getopt.h # Include targets
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
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.F.o/$*.o $@/g;' > $@
# Following command makes, e.g., f.d begin "/home/zender/obj/LINUX/f.o f.d : f.F ..."
# This works fairly well, but is a hack.
# First pattern substitutes MY_OBJ_DIR_RX, which has placeholders for slashes.
# Second pattern substitutes slashes for the placeholders.
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*(.F)?.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@
# Following command makes, e.g., f.d begin "${MY_OBJ_DIR}/f.o f.d : f.F ..."
# This would be the ideal command but I can't get the dollar sign to show up.
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.F.o/\${dollar}MY_OBJ_DIR\/$*.o $@/g;' > $@

${MY_DPN_DIR}/%.d : %.c
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/$*.o $@/g;' > $@
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

${MY_DPN_DIR}/%.d : %.cc
#	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/$*.o $@/g;' > $@
	${CPP} -M ${CPPFLAGS} $< | perl -p -e 's/$*.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' > $@

# The first LDFLAGS is for typical C programs with netCDF, math, and networking
# The second LDFLAGS enables C/Fortran linking

ifeq (${PVM_ARCH},AIX)
C++ := xlC
CC := cc
# Additional switch to fix compiler warnings on csz.c
# -qmaxmem=num Limit the amount of memory used by space intensive optimizations to <num> kilobytes
CFLAGS := -qmaxmem=8192
CPP := cc -c
# NB: Hack to make AIX look like LINUX (g77 + gcc + g++)
#CPPFLAGS := -DLINUX ${CPP_PTH} -I${NETCDF_INC}
CPPFLAGS := -D${PVM_ARCH} -DNO_NETCDF_2 ${CPP_PTH} -I${NETCDF_INC}
FC := xlf
FFLAGS := -qfixed=132
LD := ld
LDFLAGS += $(TMP_LDFLAGS) $(TMP_LIBS) -lm
ifeq ($(USE_FORTRAN_ARITHMETIC),Y)
#LDFLAGS += -lxlf90 -lxlf
endif
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O
# CPP := ${CPP} ${CPPFLAGS}
# PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS += -O -NS2000 -qfixed=132
endif
ifeq (${OPTS},D)
 CFLAGS += -g
# CPP := ${CPP} ${CPPFLAGS}
# PREPROCESS.F := ${CPP} ${CPPFLAGS}
 FFLAGS := -g -NS2000 -qfixed=132
endif
# Additional flags for AIX:
# -M Generate information to be included in a "make" description file; output goes to .u file
# -c Do not send object files to the linkage editor
# -P Preprocess but do not compile; output goes to .i file
# Using -P causes additional warning messages about lm 
# Not using -P causes *.o files to be created twice
${MY_DPN_DIR}/%.d : %.c
	${CPP} -M ${CPPFLAGS} $< ;perl -p -e 's/$*\.o/${MY_OBJ_DIR_RX}\/$*.o ${MY_DPN_DIR_RX}\/$(notdir $@)/g;s/${slash_rx}/\${slash}/g' $*.u > $@ ; \
	rm -f $*.i $*.o $*.u;
endif
# endif AIX

ifeq (${PVM_ARCH},ALPHA)
C++ := cc
CC := cc
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
LD := ld
LDFLAGS += -L${MY_LIB_DIR} -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++ -lnetcdf
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
# endif gcc
ifeq (${FC},f90)
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
# endif f90
endif
# endif ALPHA

ifeq (${PVM_ARCH},CRAY)
NETCDF_INC := /usr/local/include
NETCDF_LIB := /usr/local/lib
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf
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
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS := -D${PVM_ARCH} -Di386 ${CPP_PTH} -I${NETCDF_INC}
FC := g77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf 
LDFLAGS += 
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS := -O
 FFLAGS := -O -ffixed-line-length-132
endif
ifeq (${OPTS},D)
 CFLAGS := -g
 FFLAGS := -g -ffixed-line-length-132
endif
endif
# endif LINUX

ifeq (${PVM_ARCH},SGI6)
C++ := CC
CC := cc
CPP := cpp
CPPFLAGS := -D${CCM_TRG_FLG} ${CPP_PTH} -I${NETCDF_INC}
FC := f90 -cpp
LD := ld
LDFLAGS += -mabi=32 -mips3 -L${MY_LIB_DIR} -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++
LEX := flex
LINT := lint
YACC := bison
ifeq (${OPTS},O)
 CFLAGS := -n32 -mips3 -O2
 FFLAGS := -n32 -mips3 -O2 -extend_source
endif
ifeq (${OPTS},R)
 CFLAGS := -n32 -mips3
 FFLAGS := -n32 -mips3 -extend_source
endif
ifeq (${OPTS},D)
 CFLAGS := -n32 -mips3 -g
 FFLAGS := -n32 -mips3 -g -extend_source
endif
ifeq (${OPTS},X)
 CFLAGS := -n32 -mips3 -g -trapuv
 FFLAGS := -n32 -mips3 -g -extend_source -check_bounds -trapuv
endif
 FFLAGS += -r8 -i4
endif
# endif SGI6

# SGI64, SGIMP64
ifneq (,$(findstring SGI,${PVM_ARCH}))
ifneq (,$(findstring 64,${PVM_ARCH}))
C++ := CC -LANG:std
CC := cc
#C++ := g++
#CC := gcc
CPPFLAGS += -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f90 -cpp
LD := ld
LEX := flex
LINT := lint
YACC := bison
ifeq (${CC},gcc)
LDFLAGS += -mabi=64 -mips4 -L${MY_LIB_DIR} -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++
ifeq (${OPTS},O)
 C++FLAGS := -O2 -Wall -mabi=64
 CFLAGS += -O2 -Wall -mabi=64
 FFLAGS := -64 -mips4 -O2 -extend_source
endif
ifeq (${OPTS},R)
 C++FLAGS := -Wall -mabi=64
 CFLAGS += -Wall -mabi=64
 FFLAGS := -64 -mips4 -extend_source
endif
ifeq (${OPTS},D)
 C++FLAGS := -g -Wall -mabi=64
 CFLAGS += -g -Wall -mabi=64
 FFLAGS := -64 -mips4 -g -extend_source
endif
ifeq (${OPTS},X)
 C++FLAGS := -g -O -Wall -mabi=64
 CFLAGS += -g -O -Wall -mabi=64
 FFLAGS := -64 -mips4 -g -extend_source -check_bounds -trapuv
endif
endif
# end ifeq (${CC},gcc)
ifeq (${CC},cc)
LDFLAGS += -64 -mips4 -L${MY_LIB_DIR} -lcsz_c++ -L${NETCDF_LIB} -lnetcdf_c++
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
# 1998/08/20: Using range checking results in unresolved symbol errors unless linking to this library
endif
endif
endif
endif
# endif SGI64,SGIMP64

ifeq (${PVM_ARCH},SUN4)
NETCDF_INC := /contrib/include
NETCDF_LIB := /contrib/lib
C++ := g++
CC := acc
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf
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
NETCDF_INC := /contrib/include
NETCDF_LIB := /contrib/lib
C++ := g++
CC := gcc -ansi
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf -lsunmath -lsocket -lnsl
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
NETCDF_INC := /contrib/include
NETCDF_LIB := /contrib/lib
C++ := g++
#CC := gcc -ansi -pedantic
CC := gcc 
CPPFLAGS := -D${PVM_ARCH} ${CPP_PTH} -I${NETCDF_INC}
FC := f77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf -lsunmath -lthread -lsocket -lnsl
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
NETCDF_INC := /usr/local/include
NETCDF_LIB := /usr/local/lib
C_ONLY := Y
C++ := g++
CC := gcc -ansi
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS := -D${PVM_ARCH} -Di386 ${CPP_PTH} -I${NETCDF_INC}
FC := g77
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -L${NETCDF_LIB} -lnetcdf
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
GETOPT_TEST_OBJ=	${MY_OBJ_DIR}/getopt_test.o
ifeq (${null},$(findstring LINUX,${PVM_ARCH}))
 GETOPT_TEST_OBJ := $(GETOPT_TEST_OBJ) ${MY_OBJ_DIR}/getopt.o ${MY_OBJ_DIR}/getopt1.o ${MY_OBJ_DIR}/alloca.o
endif
getopt_test:	${MY_BIN_DIR}/getopt_test
${MY_BIN_DIR}/getopt_test:	$(GETOPT_TEST_OBJ) 
getopt_test:   $(GETOPT_TEST_OBJ)
	${CC} -o $@ $(GETOPT_TEST_OBJ) ${LDFLAGS} 
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
lib :
	-cd ${MY_LIB_DIR};
lib_cln :
	-cd ${MY_LIB_DIR};
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
