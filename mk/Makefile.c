# $Id$ -*-Makefile-*-

# Usage: 
# /bin/cp /home/zender/mk/Makefile.c /fs/cgd/home0/zender/mk
# /bin/cp /home/zender/include/* /fs/cgd/home0/zender/include

# First LDFLAGS is for typical C programs with netCDF, math, and networking
# Second LDFLAGS enables C/Fortran linking

ifndef OPTS
 OPTS := O
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

# Redefine the default patterns for C compiling
$(MY_OBJ_DIR)/%.o : %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $(MY_OBJ_DIR)/$(notdir $@)
$(MY_OBJ_DIR)/%.o : %.cc
	$(C++) $(C++FLAGS) $(CPPFLAGS) -c $< -o $(MY_OBJ_DIR)/$(notdir $@)

ifeq ($(PVM_ARCH),CRAY)
C++ := g++
CC := cc
CPPFLAGS := -D$(PVM_ARCH) -I$(MY_INC_DIR) -I/usr/local/include
FC := f90
FFLAGS := -f free -N 132
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lnco -L/usr/local/lib -lnetcdf -lm
LDFLAGS := $(LDFLAGS) -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $<
	mv -f $(notdir $@) $(MY_OBJ_DIR)
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) $<
	mv -f $(notdir $@) $(MY_OBJ_DIR)
endif

ifeq ($(PVM_ARCH),LINUX)
C++ := g++
CC := gcc -ansi
# NB: nameser.h needs -Di386, but gcc is sending -Di586 (on pentiums)
CPPFLAGS := -D$(PVM_ARCH) -I$(MY_INC_DIR) -Di386 -I/usr/local/include
FC := g77
LD := ld
#LDFLAGS := -lstdc++ -L$(MY_LIB_DIR) -lnco -L/usr/local/lib -lnetcdf_c++ -lnetcdf -lm
LDFLAGS := $(LDFLAGS)
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 C++FLAGS := -O2
 CFLAGS := -O
 FFLAGS := -O -ffixed-line-length-132
endif
ifeq ($(OPTS),D)
 C++FLAGS := -gstabs
 CFLAGS := -g
 FFLAGS := -g -ffixed-line-length-132
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),RS6K)
C++ := g++
CC := gcc -ansi
CPP := /lib/cpp -P
CPPFLAGS := -D$(PVM_ARCH) -I$(MY_INC_DIR) -I/usr/local/include
FC := xlf
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lnco -L/usr/local/lib -lncaru -lnetcdf -lm
LDFLAGS := $(LDFLAGS) -lxlf90 -lxlf
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
$(MY_OBJ_DIR)/%.o : %.F
	$(CPP) $(CPPFLAGS) $< $(MY_OBJ_DIR)/$(basename $<).f 
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $(MY_OBJ_DIR)/$(basename $<).f
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SGI5)
C++ := g++
CC := gcc -ansi
CPPFLAGS := -D$(PVM_ARCH) -I$(MY_INC_DIR) -I/opt/netcdf-2.3.2/include
FC := f77
LD := ld
LDFLAGS := -L$(MY_LIB_DIR) -lnco -L/opt/netcdf-2.3.2/lib -lnetcdf -lm
LDFLAGS := $(LDFLAGS) -lF77 -lI77 -lU77 -lftn
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 CFLAGS := -O2
 FFLAGS := -O2 -e -Nl99 -extend_source
endif
ifeq ($(OPTS),D)
 CFLAGS := -g
 FFLAGS := -g -e -Nl99 -extend_source
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SGI64)
C++ := c++
CC := cc
CPPFLAGS := -D$(PVM_ARCH) -I$(MY_INC_DIR) -I/fs/local/include -I/usr/local/include
FC := f77
LD := ld
LDFLAGS := -64 -L$(MY_LIB_DIR) -lnco -L/fs/local/lib64 -lnetcdf
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
# -mips4 does not seem to be necessary, at least in C. Could it be hurting?
 CFLAGS := -64 -mips4 -g -trapuv
 FFLAGS := -64 -mips4 -g -extend_source -check_bounds -trapuv
# 98/08/20: Using range checking results in unresolved symbol errors unless linking to this library
 LDFLAGS += -L/usr/lib64 -lftn
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SUN4)
C++ := g++
CC := acc
CPPFLAGS := -D$(PVM_ARCH) -I$(MY_INC_DIR) -I/contrib/include
FC := f77
LD := ld
LDFLAGS := -lstdc++ -L$(MY_LIB_DIR) -lnco -L/contrib/lib -lnetcdf_c++ -lnetcdf -lm
LDFLAGS := -lF77 -lM77 -lresolv $(LDFLAGS) 
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 C++FLAGS := -O2
 CFLAGS := -O2
 FFLAGS := -fast -e -Nl99
endif
ifeq ($(OPTS),D)
 C++FLAGS := -g
 CFLAGS := -g
 FFLAGS := -g -e -Nl99
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SUN4SOL2)
C++ := g++
CC := gcc -ansi
CPPFLAGS := -D$(PVM_ARCH) -I$(MY_INC_DIR) -I/contrib/include
FC := f77
LD := ld
LDFLAGS := -lstdc++ -L$(MY_LIB_DIR) -lnco -L/contrib/lib -lnetcdf_c++ -lnetcdf -lsunmath -lsocket -lnsl -lm
LDFLAGS := $(LDFLAGS) -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq ($(OPTS),O)
 C++FLAGS := -O2
 CFLAGS := -O2
 FFLAGS := -fast -e -Nl99
endif
ifeq ($(OPTS),D)
 C++FLAGS := -g
 CFLAGS := -g
 FFLAGS := -g -e -Nl99
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifeq ($(PVM_ARCH),SUNMP)
C++ := g++
#CC := gcc -ansi
CC := gcc 
CPPFLAGS := -D$(PVM_ARCH) -I$(MY_INC_DIR) -I/contrib/include
FC := f77
LD := ld
LDFLAGS := -lstdc++ -L$(MY_LIB_DIR) -lnco -L/contrib/lib -lnetcdf_c++ -lnetcdf -lsunmath -lthread -lsocket -lnsl -lm 
# NB: placing Fortran libs before -lsunmath ensure fortran IEEE functions are defined
LDFLAGS := -lF77 -lM77 -lresolv $(LDFLAGS) 
LEX := flex
LINT := lint
YACC := bison
ifeq ($(OPTS),O)
 C++FLAGS := -O2
 CFLAGS := -O2
 FFLAGS := -fast -e -Nl99
endif
ifeq ($(OPTS),D)
 C++FLAGS := -gstabs
 CFLAGS := -g
 FFLAGS := -g -e -Nl99
endif
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif

ifndef $(precision)
 precision := single
endif
ifeq ($(precision),double)
ifeq ($(PVM_ARCH),SUNMP)
 FFLAGS := $(FFLAGS) -r8 -i4
endif
ifeq ($(PVM_ARCH),RS6K)
 FFLAGS := $(FFLAGS) -qREALSIZE=8 -qINTSIZE=4
endif
ifeq ($(PVM_ARCH),IRIX5)
 FFLAGS := $(FFLAGS) -r8 -i4
endif
endif
