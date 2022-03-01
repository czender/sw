#########################################################################
# RCS Identification
#########################################################################
# $Author: zender $
# $Date$
# $Id$
# $Revision$
# $Locker:  $
# $RCSfile: Makefile,v $
# $Source: /home/zender/cvs/cld/Makefile,v $
# $Id$
# $State: Exp $
#
# $Log: not supported by cvs2svn $
# Revision 1.1.1.1  1998-09-15 02:06:40  zender
# Imported sources
#
# Revision 5.4  1995/06/19  04:54:58  zender
# have made many changes. got rid of fp_in,out,err. added flat updraft
# switch. changed -S -s to both refer to ice. got rid of while loop that
# became infinite when there was no cloud. changed vertical netCDF coord.
# to altitude. fixed bug in appending date and rev. to cmdline.
#
# Revision 5.3  1995/05/24  23:11:17  zender
# updated and validated the cloud model for ansi compatability,
# removed berkeley calls, streamlined architecture dependent
# fortran --> c calling methodology. validated results against
# comps II runs.
#
# Revision 5.2  1995/05/23  00:43:16  zender
# synchronization checkin prior to temperature sensitivity study
# for SPCP stable cloud mods.
#
# Revision 5.1  1994/02/14  21:28:54  zender
# Added DMC94 homogeneous ammonium sulphate solution nucleation,
# added -D 79 SWSCF/SWCF diagnostic. About to change default shape
# to hexagonal ice crystals, but will save bullet routines.
# NetCDF routines are not keeping up with the changes....
#
# Revision 5.0  1993/08/30  01:55:35  zender
# this version was used to make the paper (preprint at least).
# about to start working on an animation.
#
# Revision 4.12  1993/07/21  23:26:28  zender
# added Liou_IR_fudge() as a default routine whenever LIOU is true.
#
# Revision 4.11  1993/06/19  22:05:25  zender
# added ebert&curry parameterization to LW computations, and
# made this the default for the region outside of the cloud grid.
#
# Revision 4.10  1993/06/11  04:35:31  zender
# Played w/ Makefiles, made all CRAY2 refs. into plain CRAY,
# and put all shavano stuff on peace. changed surface type to
# land for predictable albedos, changed tropics SAZ to 30 deg.,
# removed gaussians from initial distributions. changed plev = 113.
# slingo option doesn't work for some reason (-D 53)
#
# Revision 4.9  1993/06/09  00:59:14  zender
# prepared the -Z initial_distribution_type switch to handle everyone's
# different clouds. Made Knollenberg fit his observations, prepared
# optical properties for 50 sizes and removed all num_sizes == 40 specific
# stuff. changed tau_CCN to 1. s.
#
#########################################################################  
ifndef $(OPTS)
 OPTS := O
endif
#
ifeq ($(PVM_ARCH),RS6K)
CPP := /lib/cpp -P
CC := gcc
CPPFLAGS := -I/usr/local/include
C++ := g++
#FC := xlf90
FC := xlf
LD := ld
LINT := lint
endif
#
ifeq ($(PVM_ARCH),SUNMP)
CC := acc
#CC := gcc
CPPFLAGS := -I/opt/netcdf2.3.2/include
C++ := g++
FC := f77
LD := ld
LINT := lint
endif
#
ifeq ($(PVM_ARCH),SGI5)
CC := gcc
#CC := gcc
CPPFLAGS := -I/opt/netcdf-2.3.2/include
C++ := g++
FC := f77
LD := ld
LINT := lint
endif
#
ifeq ($(PVM_ARCH),CRAY)
CC := cc
CPPFLAGS := -I/usr/local/include
C++ := g++
FC := cf77
LD := ld
LINT := lint
endif
#
ifeq ($(OPTS),O)
 STRING := Using Optimal compilation switches . . .
ifeq ($(PVM_ARCH),SUNMP)
 CFLAGS := -O2
 FFLAGS := -fast -e -Nl99
# FFLAGS := -O -e -Nl99
endif
ifeq ($(PVM_ARCH),SGI5)
 CFLAGS := -O2
 FFLAGS := -O2 
endif
ifeq ($(PVM_ARCH),RS6K)
 CFLAGS := -O2
 CPP := $(CPP) $(CPPFLAGS)
 PREPROCESS.F := $(CPP) $(CPPFLAGS)
 FFLAGS := -O -NS2000 -qfixed=132
endif
ifeq ($(PVM_ARCH),CRAY)
 CFLAGS := -O
endif
endif
#
ifeq ($(OPTS),D)
 STRING := Using Debugging compilation switches . . .
ifeq ($(PVM_ARCH),SUNMP)
 CFLAGS := -g
 FFLAGS := -g -e -Nl99
endif
ifeq ($(PVM_ARCH),SGI5)
 CFLAGS := -g
 FFLAGS := -g
endif
ifeq ($(PVM_ARCH),RS6K)
 CFLAGS := -g
 CPP := $(CPP) $(CPPFLAGS)
 PREPROCESS.F := $(CPP) $(CPPFLAGS)
 FFLAGS := -g -NS2000 -qfixed=132
endif
ifeq ($(PVM_ARCH),CRAY)
 CFLAGS := -g
endif
endif
#
ifeq ($(OPTS),R)
 STRING := Using Regular compilation switches . . .
 CFLAGS := 
 FFLAGS := -e -Nl99
endif
#
ifeq ($(PVM_ARCH),SUNMP)
 LDFLAGS := -lF77 -lM77
endif
ifeq ($(PVM_ARCH),SGI5)
 LDFLAGS := -lF77 -lI77 -lU77 -lftn
 FFLAGS := $(FFLAGS) -col120 -Nl99
endif
ifeq ($(CC),gcc)
 CFLAGS := $(CFLAGS)
endif
ifeq ($(PVM_ARCH),SUNMP)
# LDFLAGS := $(LDFLAGS) -lsunmath -lthread -lsocket -lnsl
 LDFLAGS := $(LDFLAGS) -lsunmath
 CPPFLAGS := $(CPPFLAGS) -I/opt/pvm3/include
endif
ifeq ($(PVM_ARCH),SUN4SOL2)
# LDFLAGS := $(LDFLAGS) -lsunmath -lsocket -lnsl
 LDFLAGS := $(LDFLAGS) -lsunmath
# CPPFLAGS := $(CPPFLAGS) -I/opt/pvm3/include
endif
ifeq ($(OS),SunOS)
# CPPFLAGS := $(CPPFLAGS) -I/opt/pvm3/include
endif
ifeq ($(OS),AIX)
# LDFLAGS := -L/usr/local/lib -L/usr/local/pvm3/lib/RS6K $(LDFLAGS) -lncaru
 LDFLAGS := -L/usr/local/lib $(LDFLAGS) -lncaru -lxlf90 -lxlf
# CPPFLAGS := $(CPPFLAGS) -I/usr/local/pvm3/include -I/usr/local/pvm3/pvmgs
 CPPFLAGS := $(CPPFLAGS)
endif
ifeq ($(OS),UNICOS)
 LDFLAGS := -L/usr/local/lib -L/lib -lf
 CPPFLAGS := $(CPPFLAGS)
endif
ifeq ($(PVM_ARCH),SUNMP)
 LDFLAGS := -cg92 -L/opt/SUNWspro/SC3.0/lib/cg92 -L/opt/netcdf2.3.2/lib $(LDFLAGS) 
 FFLAGS := $(FFLAGS) -cg92
endif
ifeq ($(PVM_ARCH),SGI5)
 LDFLAGS := -L/opt/netcdf-2.3.2/lib $(LDFLAGS) 
endif
ifeq ($(CC),acc)
 CFLAGS := $(CFLAGS) -cg92
endif

 LDFLAGS := $(LDFLAGS) -lnetcdf -lm
 CPPFLAGS := $(CPPFLAGS) -D$(PVM_ARCH) -D$(OS) -D$(PVM_ARCH)

########################################################################
# Old stuff...should be updated for hetero-cpu environment
########################################################################
ATHENA_LIBS= -lXaw -lXmu -lXt -lXext -lX11
MOTIF_CFLAGS= -D_NO_PROTO -I/usr/local/sun4/x11r5/include
########################################################################
# End old stuff
########################################################################
#
$(MY_OBJ_DIR)/%.o : %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
#
ifeq ($(PVM_ARCH),CRAY)
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $<
	mv -f $(notdir $@) $(MY_OBJ_DIR)
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) $<
	mv -f $(notdir $@) $(MY_OBJ_DIR)
endif
#
ifeq ($(PVM_ARCH),SUNMP)
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif
#
ifeq ($(PVM_ARCH),SGI5)
$(MY_OBJ_DIR)/%.o : %.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif
#
ifeq ($(PVM_ARCH),RS6K)
$(MY_OBJ_DIR)/%.o : %.F
	$(CPP) $(CPPFLAGS) $< $(MY_OBJ_DIR)/$(basename $<).f 
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $(MY_OBJ_DIR)/$(basename $<).f
$(MY_OBJ_DIR)/%.o : %.f
	$(FC) -c $(FFLAGS) -o $(MY_OBJ_DIR)/$(notdir $@) $<
endif
#
IDX_RFR_SRC=	idx_rfr.c idx_rfr_h2o_ice_sgw.f flx_slr_frc.f 
IDX_RFR_OBJ=	$(MY_OBJ_DIR)/idx_rfr.o \
		$(MY_OBJ_DIR)/idx_rfr_h2o_ice_sgw.o \
		$(MY_OBJ_DIR)/flx_slr_frc.o 
idx_rfr:	$(IDX_RFR_OBJ)
	$(CC) $(CPP_FLAGS) -o $(MY_BIN_DIR)/idx_rfr $(IDX_RFR_OBJ) $(LDFLAGS) 
	chmod 744 $(MY_BIN_DIR)/idx_rfr
#
CLD_SRC=	main.c utilities.c etbfct.f ice.c mie_bohren_c.f \
		idx_rfr_h2o_ice_sgw.f ccm2rad_c.f interp.c spline.f 
CLD_C_SRC=	main.c utilities.c ice.c interp.c 
CLD_F_SRC=	etbfct.f mie_bohren_c.f idx_rfr_h2o_ice_sgw.f ccm2rad_c.f spline.f
CLD_OBJ=	$(MY_OBJ_DIR)/main.o \
		$(MY_OBJ_DIR)/utilities.o \
		$(MY_OBJ_DIR)/etbfct.o \
		$(MY_OBJ_DIR)/ice.o \
		$(MY_OBJ_DIR)/mie_bohren_c.o \
		$(MY_OBJ_DIR)/idx_rfr_h2o_ice_sgw.o \
		$(MY_OBJ_DIR)/ccm2rad_c.o \
		$(MY_OBJ_DIR)/interp.o \
		$(MY_OBJ_DIR)/spline.o
cld:	$(CLD_OBJ) defs.h globals.h
#	$(LINT) $(LINTFLAGS) $(CLD_C_SRC)
	@echo $(STRING)
	$(CC) -o $(MY_BIN_DIR)/cld $(CLD_OBJ) $(LDFLAGS) 

MOTIF_LIBS= -lXm -lXt -lX11
XCFLAGS= -I/usr/local/sun4/x11r5/include -D_NO_PROTO
XCLD_SRC=	xcld.c xcldicon.h
XCLD_OBJ=	$(MY_OBJ_DIR)/xcld.o
xcld:	$(XCLD_OBJ) 
		$(CC) -o $(MY_BIN_DIR)/xcld $(XCLD_OBJ) $(LDFLAGS) -lXm -lXt -lX11

$(MY_OBJ_DIR)/xcld.o:	xcld.c
		$(CC) -I  -c xcld.c

CFC_TEST_SRC=	cfc_test.c cfc_test_f.f 
CFC_TEST_OBJ=	$(MY_OBJ_DIR)/cfc_test.o \
		$(MY_OBJ_DIR)/cfc_test_f.o 
cfc_test:	$(CFC_TEST_OBJ) 
#	$(CC) -o $(MY_BIN_DIR)/cfc_test $(CFC_TEST_OBJ) $(LDFLAGS) $(SUBSET_F77LIBS)
#	$(CC) -o $(MY_BIN_DIR)/cfc_test $(CFC_TEST_OBJ) $(LDFLAGS) 
#	$(NCARGCC) -o $(MY_BIN_DIR)/cfc_test $(CFC_TEST_OBJ) $(LDFLAGS) 
	$(CC) -o $(MY_BIN_DIR)/cfc_test $(CFC_TEST_OBJ) $(LDFLAGS) $(BIGGER_SUBSET_F77LIBS)

GREG_SRC=	etbfct.f fct.f $(HOME)/ncarg/ncargf.f
GREG_OBJ=	etbfct.o fct.o $(HOME)/ncarg/ncargf.o
greg:	$(GREG_OBJ) 
	$(FC) -o $(MY_BIN_DIR)/greg $(GREG_OBJ) $(LDFLAGS) 

INTERP_TEST_OBJ=	interp_test.o \
			utilities.c \
			/data2/zender/nr/ratint.c \
			/data2/zender/nr/polint.c \
			/data2/zender/nr/nrutil.c\
			spline.o 
interp_test:	$(INTERP_TEST_OBJ) 
	$(CC) $(INTERP_TEST_OBJ) -o $(MY_BIN_DIR)/interp_test $(BIGGER_SUBSET_F77LIBS) $(LDFLAGS) 

IOGRAPHS_SRC=	iographs.c utilities.c ncarg_colors.f
IOGRAPHS_OBJ=	iographs.o utilities.o ncarg_colors.o
iographs:	$(IOGRAPHS_OBJ) 
#		$(LINT) $(LINTFLAGS) $(IOGRAPHS_SRC)
		$(NCARGCC) $(IOGRAPHS_OBJ) $(NETCDFLIB) -o $(MY_BIN_DIR)/iographs $(LDFLAGS) 

IOGRAPHS.3.2_SRC=	iographs.c utilities.c ncarg_colors.f
IOGRAPHS.3.2_OBJ=	iographs.o utilities.o ncarg_colors.o
iographs.3.2:	$(IOGRAPHS.3.2_OBJ) 
#		$(LINT) $(LINTFLAGS) $(IOGRAPHS.3.2_SRC)
		ncargcc $(IOGRAPHS.3.2_OBJ) $(NETCDFLIB) -o $(MY_BIN_DIR)/iographs $(LDFLAGS) 

IO_DIFF_SRC=	io_diff.c utilities.c 
IO_DIFF_OBJ=	io_diff.o utilities.o 
io_diff:	$(IO_DIFF_OBJ) 
		$(NCARGCC) $(IO_DIFF_OBJ) $(NETCDFLIB) -o $(MY_BIN_DIR)/io_diff $(LDFLAGS) 

IO_MOVIE_SRC=	io_movie.c utilities.c ncarg_colors.f
IO_MOVIE_OBJ=	io_movie.o utilities.o ncarg_colors.o
io_movie:	$(IO_MOVIE_OBJ) 
#		$(LINT) $(LINTFLAGS) $(IO_MOVIE_SRC)
		ncargcc $(IO_MOVIE_OBJ) $(HOME)/ncarg/libncargC.a $(HOME)/ncarg/libncarg_gksC.a $(LDFLAGS) $(NETCDFLIB) -o $(MY_BIN_DIR)/io_movie

IO_MOVIE.3.2_SRC=	io_movie.3.2.c utilities.c ncarg_colors.f
IO_MOVIE.3.2_OBJ=	io_movie.3.2.o utilities.o ncarg_colors.o
io_movie.3.2:	$(IO_MOVIE.3.2_OBJ) 
#		$(LINT) $(LINTFLAGS) $(IO_MOVIE.3.2_SRC)
		ncargcc $(IO_MOVIE.3.2_OBJ) $(NETCDFLIB) -o $(MY_BIN_DIR)/io_movie.3.2 $(LDFLAGS) 

IO_TEST_SRC=	io_test.c utilities.c
IO_TEST_OBJ=	io_test.o utilities.o
io_test:	$(IO_TEST_OBJ) 
		$(NCARGCC) $(IO_TEST_OBJ) $(NETCDFLIB) -o $(MY_BIN_DIR)/io_test $(LDFLAGS) 

MIE_BOHREN_SRC=	mie_bohren.f
MIE_BOHREN_OBJ=	$(MY_OBJ_DIR)/mie_bohren.o
mie_bohren:	$(MIE_BOHREN_OBJ) 
	$(FC) $(MIE_BOHREN_OBJ) -o $(MY_BIN_DIR)/mie_bohren $(LDFLAGS)

F_TEMPLATE_SRC := f_template.F date_time.c
F_TEMPLATE_OBJ := $(MY_OBJ_DIR)/f_template.o \
		$(MY_OBJ_DIR)/date_time.o
f_template:     $(F_TEMPLATE_OBJ)
	$(FC) -o $(MY_BIN_DIR)/f_template $(F_TEMPLATE_OBJ) $(LDFLAGS) 
	chmod 744 $(MY_BIN_DIR)/f_template

clean:
	rm -f $(CLD_OBJ)
	rm -f $(IDX_RFR_OBJ)
	rm -f $(XCLD_OBJ)
	rm -f $(CFC_TEST_OBJ)
	rm -f $(GREG_OBJ)
	rm -f $(INTERP_TEST_OBJ)
	rm -f $(IOGRAPHS_OBJ)
	rm -f $(IOGRAPHS.3.2_OBJ)
	rm -f $(IO_DIFF_OBJ)
	rm -f $(IO_MOVIE_OBJ)
	rm -f $(IO_MOVIE.3.2_OBJ)
	rm -f $(IO_TEST_OBJ)
	rm -f $(MIE_BOHREN_OBJ)
# DO NOT DELETE THIS LINE -- make depend depends on it.








