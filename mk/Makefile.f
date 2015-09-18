# $Id$ -*-Makefile-*-

# First LDFLAGS is for typical C programs with netCDF, math, and networking
# Second LDFLAGS enables C/Fortran linking

# Usage: 
# /bin/cp $HOME/mk/Makefile.f /fs/cgd/home0/zender/mk
# /bin/cp $HOME/include/* /fs/cgd/home0/zender/include
# rcp -p $HOME/mk/Makefile.f sanitas.cgd.ucar.edu:$HOME/mk

ifndef OPTS
 OPTS := O
endif
# NB: Be careful using the -r8 -i4 switches. 
# I rewrote most RT programs to compile under LINUX g77, which doesn't have -r8 -i4.
# So instead I use the preprocessor method of declaring tokens like "COMPUTATIONAL_PRECISION".
# This method is not transparently easy to use in conjunction with the -r8 -i4 method.
# The errors that result from using both together can be subtle and hard to trace.
# In general I think it safest to use a token, or "-r8 -i4", but not both.
ifndef $(precision)
 precision := single
endif
ifndef PVM_ARCH
 PVM_ARCH := $(shell pvmgetarch)
endif
ifndef MY_OBJ_DIR
 MY_OBJ_DIR := .
endif
ifndef MY_INC_DIR
 MY_INC_DIR := .
endif
ifndef MY_LIB_DIR
 MY_LIB_DIR := .
endif
ifndef MY_BIN_DIR
 MY_BIN_DIR := .
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
ifdef GSL_LIB
 GSL_LIB := -L$(GSL_LIB)
endif

# Source file names with directories removed
MDL_SRC := 
# Directories to search for source files
MDL_PTH := . $(HOME)/include
# Dependency list for executable
MDL_OBJ := $(addprefix $(MY_OBJ_DIR)/,$(addsuffix .o, $(basename $(MDL_SRC)))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix $(MY_DPN_DIR)/,$(addsuffix .d, $(basename $(MDL_SRC)))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.d file
VPATH := $(subst $(space),:,$(MDL_PTH))
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,$(MDL_PTH),-I$(dir))

# Redefine default C and C++ pattern rules
$(MY_OBJ_DIR)/%.o : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $(MY_OBJ_DIR)/$(notdir $@)
$(MY_OBJ_DIR)/%.o : %.cc
	$(C++) $(CPPFLAGS) $(C++FLAGS) -c $< -o $(MY_OBJ_DIR)/$(notdir $@)

# Default Fortran pattern rules: CRAY always overrides these rules, LINUX and SGIMP64 do for F90
$(MY_OBJ_DIR)/%.o : %.F90
#	$(CPP) $(CPPFLAGS) $< $(MY_OBJ_DIR)/$(patsubst %.F90,%.f90,$(notdir $<))
#	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $(MY_OBJ_DIR)/$(patsubst %.F90,%.f90,$(notdir $<))
#	$(FC) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f90
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<

ifeq ($(PVM_ARCH),CRAY)
C++ := g++
CC := cc
CPP := cpp
CPPFLAGS := -D$(PVM_ARCH) -I./ -I$(MY_INC_DIR) -I/usr/local/include
FC := f90
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lnco -L/usr/local/lib -lnetcdf -lm
LDFLAGS += -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O
 FFLAGS = -N 132
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS = -g -N 132
endif
ifeq ($(OPTS),X)
 CFLAGS := -g
 FFLAGS = -g -N 132 -e i
endif
ifeq ($(precision),double)
 FFLAGS += -DDOUBLE_PRECISION
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(CPP) -P $(CPPFLAGS) $< > $(patsubst %.F,%.f,$(notdir $<))
	$(FC) -c $(FFLAGS) $(patsubst %.F,%.f,$(notdir $<)) 
	-mv -f $(notdir $@) $(MY_OBJ_DIR)
	rm -f $(patsubst %.F,%.f,$(notdir $<)) 
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) $<
	mv -f $(notdir $@) $(MY_OBJ_DIR)
endif

ifeq ($(PVM_ARCH),LINUX)
C++ := g++
CC := gcc
CPPFLAGS := -D$(PVM_ARCH)  -I./ -I$(MY_INC_DIR) -I$(NETCDF_INC)
FC := pgf90
LD := ld
LDFLAGS := -mp -L$(MY_LIB_DIR) -lcsz_f77 -L$(NETCDF_LIB) -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq ($(CC),gcc)
ifeq ($(OPTS),O)
 CFLAGS := -O -Wall
endif
ifeq ($(OPTS),D)
 CFLAGS := -g -Wall
endif
ifeq ($(OPTS),R)
 CFLAGS := -Wall
endif
ifeq ($(OPTS),X)
 CFLAGS := -g -O -Wall
endif
 C++FLAGS := $(CFLAGS)
endif
ifeq ($(FC),pgf90)
ifeq ($(OPTS),O)
 FFLAGS := -fast -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign
endif
ifeq ($(OPTS),D)
 FFLAGS := -g -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign
endif
ifeq ($(OPTS),R)
 FFLAGS := -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign
endif
ifeq ($(OPTS),X)
 FFLAGS := -g -Mbounds -Mextend -Mnosecond_underscore -mp -byteswapio -Mrecursive -Mdalign
endif
endif
ifeq ($(FC),g77)
ifeq ($(OPTS),O)
 FFLAGS := -O -ffixed-line-length-132 -fno-second-underscore
endif
ifeq ($(OPTS),D)
 FFLAGS := -g -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq ($(OPTS),R)
 FFLAGS := -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge
endif
ifeq ($(OPTS),X)
 FFLAGS := -g -O -ffixed-line-length-132 -fno-second-underscore -fdebug-kludge -fbounds-check
endif
endif
endif

ifeq ($(PVM_ARCH),RS6K)
C++ := g++
CC := gcc -ansi
CPP := /lib/cpp -P
CPPFLAGS := -D$(PVM_ARCH) -I./ -I$(MY_INC_DIR) -I/usr/local/include
FC := xlf
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lnco -L/usr/local/lib -lncaru -lnetcdf -lm
LDFLAGS += -lxlf90 -lxlf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 CPP := $(CPP) $(CPPFLAGS)
 PREPROCESS.F := $(CPP) $(CPPFLAGS)
 FFLAGS := -O -NS2000 -qfixed=132
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 CPP := $(CPP) $(CPPFLAGS)
 PREPROCESS.F := $(CPP) $(CPPFLAGS)
 FFLAGS := -g -NS2000 -qfixed=132
endif
ifeq ($(precision),double)
 FFLAGS += -DDOUBLE_PRECISION -qREALSIZE=8 -qINTSIZE=4
endif
ifeq ($(precision),double)
 FFLAGS += -DDOUBLE_PRECISION
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(CPP) $(CPPFLAGS) $< $(MY_OBJ_DIR)/$(basename $<).f 
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $(MY_OBJ_DIR)/$(basename $<).f
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SGI5)
C++ := g++
CC := gcc -ansi
CPPFLAGS := -D$(PVM_ARCH) -I./ -I$(MY_INC_DIR) -I/opt/netcdf-2.3.2/include
FC := f77
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lnco -L/opt/netcdf-2.3.2/lib -lnetcdf -lm
LDFLAGS += -lF77 -lI77 -lU77 -lftn
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -O2 -e -Nl200 -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e -Nl200 -extend_source
endif
ifeq ($(precision),double)
 FFLAGS += -DDOUBLE_PRECISION
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SGI64)
C++ := c++
CC := cc
CPPFLAGS := -D$(PVM_ARCH) -I./ -I$(MY_INC_DIR) -I/fs/local/include -I/usr/local/include
FC := f77
LD := ld
LDFLAGS := -64 -L$(MY_LIB_DIR) -L/fs/local/lib64 -lnetcdf
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 CFLAGS := -64 -O2
 FFLAGS := -64 -O2 -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -64 -g
 FFLAGS := -64 -g -extend_source
endif
ifeq ($(OPTS),X)
 CFLAGS := -64 -g -trapuv
 FFLAGS := -64 -g -extend_source -check_bounds -trapuv
endif
ifeq ($(precision),double)
 FFLAGS += -DDOUBLE_PRECISION -r8 -i4
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SGIMP64)
C++ := CC
CC := cc
CPP := cpp
# spoof SGI64
CPPFLAGS := -DSGI64 $(CPP_PTH) -I$(NETCDF_INC)
FC := f90
LD := ld
# 20000706: Added -mp -mpio for OpenMP compliance
LDFLAGS := -64 -mp -mips4 -L$(MY_LIB_DIR) -lcsz_f77 -L/fs/cgd/home0/zender/lib/SGI64 -lrecipes_f -L$(NETCDF_LIB) -lnetcdf -L/contrib/lib -lspecfun
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 CFLAGS := -64 -mips4 -O2
 FFLAGS := -64 -mips4 -O2 -extend_source -mp -mpio 
endif
ifeq ($(OPTS),R)
 CFLAGS := -64 -mips4
 FFLAGS := -64 -mips4 -extend_source -mp -mpio 
endif
ifeq ($(OPTS),D)
 CFLAGS := -64 -mips4 -g
 FFLAGS := -64 -mips4 -g -extend_source -mp -mpio 
endif
ifeq ($(OPTS),X)
 CFLAGS := -64 -mips4 -g -trapuv
 FFLAGS := -64 -mips4 -g -extend_source -check_bounds -trapuv -mp -mpio 
endif
# 20000706: -DHIDE_SHR_MSG eliminates messages like "fff.F90", line 481: Warning: Referenced scalar variable THR_NBR is SHARED by default
$(MY_OBJ_DIR)/%.o : %.F90
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -DHIDE_SHR_MSG -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif
# endif SGIMP64

ifeq ($(PVM_ARCH),SUN4)
C++ := g++
CC := acc
CPPFLAGS := -D$(PVM_ARCH) -I./ -I$(MY_INC_DIR) -I/contrib/include
FC := f77
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lcsz_f77 -lnco -L/contrib/lib -lnetcdf -lm
LDFLAGS += -cg92 -L/opt/SUNWspro/SC3.0/lib/cg92 -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -fast -e -Nl200
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e -Nl200
endif
ifeq ($(precision),double)
 FFLAGS += -DDOUBLE_PRECISION
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SUN4SOL2)
C++ := g++
CC := gcc -ansi
CPPFLAGS := -D$(PVM_ARCH) -I./ -I$(MY_INC_DIR) -I/contrib/include
FC := f77
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lcsz_f77 -lnco -L/contrib/lib -lnetcdf -lsunmath -lsocket -lnsl -lm
LDFLAGS += -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -fast -e
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e
endif
ifeq ($(OPTS),X)
 CFLAGS := -g 
 FFLAGS := -g -e
# NB: 98/06/01 -C (range-checking) is not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += -C
endif
endif
ifeq ($(FC),f77)
 FFLAGS += -Nl200
endif
ifeq ($(precision),double)
 FFLAGS += -DDOUBLE_PRECISION -r8 -i4
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SUNMP)
C++ := g++
CC := gcc -ansi
CPPFLAGS := -D$(PVM_ARCH) -I./ -I$(MY_INC_DIR) -I/contrib/include
FC := f90
#FC := g77
#FC := g77
# libs in /pkg/craysoft/f90e/lib
#FC := /usr/local/bin/f90 
LD := ld
#LDFLAGS := -L$(MY_LIB_DIR) -lcsz_f77 -lnco -L/contrib/lib -lnetcdf -lsunmath -lthread -lsocket -lnsl -lm
#LDFLAGS += -lF77 -lM77 -lresolv
LDFLAGS := -L/contrib/lib -lnetcdf -lsunmath -lnsl
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -fast -e
endif
ifeq ($(OPTS),D)
 CFLAGS := -g 
 FFLAGS := -g -e
endif
ifeq ($(OPTS),X)
 CFLAGS := -g 
 FFLAGS := -g -e
# NB: 98/06/01 -C (range-checking) is not supported by Sun f90
ifeq ($(FC),f77)
 FFLAGS += -C
endif
endif
ifeq ($(FC),f77)
 FFLAGS += -Nl200
endif
ifeq ($(precision),double)
 FFLAGS += -DDOUBLE_PRECISION -r8 -i4
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif
