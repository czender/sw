# $Id$ -*-makefile-*-

# Purpose: Makefile for radiation codes

# Copyright (C) 1994--present Charlie Zender
# You may copy, distribute, and/or modify this software under the terms of the GNU General Public License (GPL) Version 3
# See http://www.gnu.org/copyleft/gpl.html for full license text

# Quickie test copies:
# scp ~/sw/aca/Makefile givre.ess.uci.edu:aca/Makefile	
# scp ~/sw/aca/Makefile greenplanet.ps.uci.edu:aca/Makefile	

# Compilation:
# swnb2 is intimately linked to the fortran library libcsz_f90
# These two pieces must be compiled with the same intrinsic precisions
# Debugging compilation:
# cd ~/sw/f  ;make cln;make OMP=N OPTS=X PRC=D
# cd ~/sw/aca;make cln;make OMP=N OPTS=X PRC=D
# Production compilation:
# cd ~/sw/f  ;make cln;make OMP=N OPTS=O PRC=D
# cd ~/sw/aca;make cln;make OMP=N OPTS=O PRC=D
# cd ~/sw/f  ;make cln;make OMP=N OPTS=D MPI=Y PRC=D # greenplanet NB: Intel MPI IPO does not work with -ipo option set by ICC_RCH_ARGS when OPTS=O
# cd ~/sw/aca;make cln;make OMP=N OPTS=D MPI=Y PRC=D # greenplanet NB: Intel MPI IPO does not work with -ipo option set by ICC_RCH_ARGS when OPTS=O
# cd ~/sw/f  ;make cln;make OPTS=D MPI=Y OMP=N PRC=D # NB: Intel MPI IPO does not work with -ipo option set by ICC_RCH_ARGS when OPTS=O
# cd ~/sw/aca;make cln;make OPTS=D MPI=Y OMP=N PRC=D # NB: Intel MPI IPO does not work with -ipo option set by ICC_RCH_ARGS when OPTS=O

include ${HOME}/mk/Makefile.tex

ifndef MY_OBJ_DIR
 MY_OBJ_DIR := ../obj
endif
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

# Primary tokens which determine build options
# Specify non-default when invoking make, e.g. make DODS=Y
ifndef ${ABI}
# 32- vs. 64-bit ABI: 32=32-bit mode, 64=64-bit mode (default) if available
 ABI := 64
endif # endif ABI
ifndef DATA_RT
 DATA_RT := ${DATA}/aca
endif # endif DATA_RT
ifndef ${DBG}
# Debugging token N=No (default) Y=Yes
 DBG := N
endif # endif DBG
ifndef DPN_SFX
 DPN_SFX := d
endif # endif DPN_SFX
ifndef GCC_RCH_ARG
 GCC_RCH_ARG :=
endif # endif GCC_RCH_ARG
ifndef ICC_RCH_ARG
 ICC_RCH_ARG :=
endif # endif ICC_RCH_ARG
ifndef MDL_SFX
 MDL_SFX := mod
endif # endif MDL_SFX
ifndef ${MPI} # MPI
 MPI := N
endif # endif MPI
ifndef ${NETCDF4} # netCDF4 support
 NETCDF4 := Y
endif # endif NETCDF4
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
ifndef ${NETCDFF} # Link to libnetcdff.a
 NETCDFF := Y
endif # endif NETCDFF
ifndef ${OMP} # OpenMP
 OMP := N
endif # endif OMP
ifndef OPTS
 OPTS := O
endif
ifndef ${PRC}
# Precision token D=Double S=Single N=Native
 PRC := D
endif # endif PRC
ifndef PVM_ARCH
 PVM_ARCH := $(shell ${HOME}/sh/pvmgetarch)
endif
ifndef STC
# Created statically linked executable
 STC := N
endif
ifndef SZ
# Link to Szip library
 SZ := N
endif
ifndef UNAMES
 UNAMES := $(shell uname -s)
endif
ifndef USR_TKN
 USR_TKN := 
endif # endif USR_TKN

ifeq (${NETCDF4},Y)
 CPPFLAGS += -DENABLE_NETCDF4
 USR_TKN += -DENABLE_NETCDF4
endif # endif NETCDF4

# Derived-tokens based on primary tokens
# These tokens should not be altered by hand
ifeq (${PRC},D)
 USR_TKN := -DPRC_DBL ${USR_TKN}
endif # !PRC
ifeq (${PRC},L)
 USR_TKN := -DPRC_LDB ${USR_TKN}
endif # !PRC
ifeq (${PRC},S)
 USR_TKN := -DPRC_FLT ${USR_TKN}
endif # !PRC
ifeq (${PRC},N)
 USR_TKN := -DPRC_NTV ${USR_TKN}
endif # !PRC

# NB: CPP_TKN is recursively expanded variable, define only when components are ready
CPP_TKN = ${USR_TKN} -D${PVM_ARCH} -DHOSTNAME='${HOST}' -DUSER='${USER}'

# Default build. DAP, NETCDF4, and PNETCDF clauses may overwrite this
NC_LDFLAGS := -L${NETCDF_LIB}
NC_LIBS := -lnetcdf -lcurl

ifeq (${NETCDF4},Y)
# Enable netCDF4 functionality
ifndef HDF5_ROOT
 HDF5_ROOT := ${NETCDF4_ROOT}
endif
ifndef HDF5_INC
 HDF5_INC := ${HDF5_ROOT}/include # Directory containing hdf5.h
endif
ifndef HDF5_LIB
 HDF5_LIB := ${HDF5_ROOT}/lib # Directory containing libhdf5.a
endif
# fxm: Temporary kludge: Machines with MPI_ROOT set are assumed to have parallel filesystems for NCO builds
ifdef MPI_ROOT
 HDF5_INC += -I${MPI_ROOT}/include # Directory containing mpi.h
 HDF5_LIB += ${MPI_ROOT}/lib64 # Directory containing libmpi.a
endif # !MPI_ROOT
 NETCDF_INC := ${NETCDF4_ROOT}/include $(shell nc-config --fflags)
 NETCDF_LIB := ${NETCDF4_ROOT}/lib
 NC_LDFLAGS := -L${NETCDF_LIB} $(addprefix -L,${HDF5_LIB})
# NC_LIBS := -lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl
 NC_LIBS := $(shell nc-config --flibs)
ifdef MPI_ROOT
 NC_LIBS += -lmpi
endif # !MPI_ROOT
ifdef CURL_LIB
 NC_LDFLAGS += -L${CURL_LIB}
endif # end 
ifeq (${SZ},Y)
ifndef SZ_LIB
 SZ_LIB := /usr/lib # Directory containing libsz.a
endif
 NC_LDFLAGS += -L${SZ_LIB}
 NC_LIBS += -lsz
endif # !SZ
ifneq (${HDF5_INC},${NETCDF_INC})
 NETCDF_INC += -I${HDF5_INC}
endif # end 
endif # endif NETCDF4

# Link to libnetcdff.a
ifeq (${NETCDFF},Y)
 NC_LIBS := -lnetcdff ${NC_LIBS}
endif # endif NETCDFF

# Derived-tokens based on primary tokens
# These tokens should not be altered by hand
# Works on LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, FREEBSD, MACOSX, and WIN32
ifneq (${null},$(findstring ${PVM_ARCH},LINUXALPHALINUXAMD64LINUXARMFREEBSDMACOSXWIN32))
# Decide among the plethora of Linux compilers
ifndef LINUX_CXX
# C++ compiler for Linux
 LINUX_CXX := g++
 #LINUX_CXX := como
 #LINUX_CXX := icc
 #LINUX_CXX := insure
 #LINUX_CXX := pathCC
 #LINUX_CXX := pgCC
endif # endif LINUX_CXX
ifndef LINUX_CC
# C compiler for Linux
 LINUX_CC := gcc -std=c99 -pedantic -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := como --c99
 #LINUX_CC := icc -std=c99 -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := insure
 #LINUX_CC := pathcc -std=c99 -D_BSD_SOURCE -D_POSIX_SOURCE
 #LINUX_CC := pgcc
endif # endif LINUX_CC
ifndef LINUX_FC
# Fortran compiler for Linux
 #LINUX_FC := g95
 LINUX_FC := gfortran
 #LINUX_FC := gfortran -std=f95
 #LINUX_FC := lf95
 #LINUX_FC := ifort
 #LINUX_FC := pathf95
 #LINUX_FC := pgf95
endif # endif LINUX_FC
endif # endif LINUX

# NB: Do NOT add comment lines, e.g., # This is a comma, to character definitions
null :=
space := ${null} ${null}
comma := ,
newline := \n

# Source file names with directories removed
MDL_SRC := H2OH2O.F90 NO2.F90 O2X.F90 O2.F90 O3.F90 CFC11.F90 CFC12.F90 clm.F90 htrn2nb.F machine.F mie2nc.F nbm.F swnb2.F90 swnb3.F90 bdref.f d1mach.f disort.f disort2.f linpak.f r1mach.f rdhtrn.f
# disort3.f disobrdf3.f bdref3.f errpack3.f lapack3.f linpak3.f rdi1mach3.f disotest3.f disobrdf3.f
# Directories to search for source files
MDL_PTH := . ${HOME}/include ${HOME}/sw/f
# Dependency list for executable
MDL_OBJ := $(addprefix ${MY_OBJ_DIR}/,$(addsuffix .o, $(basename ${MDL_SRC}))) 
# Dependency (make) file for each object file
MDL_DPN := $(addprefix ${MY_DPN_DIR}/,$(addsuffix .${DPN_SFX}, $(basename ${MDL_SRC}))) 
# VPATH helps make find dependencies (which are not pathname qualified) in *.${DPN_SFX} file
VPATH := $(subst ${space},:,${MDL_PTH})
# Prepend -I to use for compiler argument
CPP_PTH := $(foreach dir,${MDL_PTH},-I${dir})

# Variables having to do with binary executables created by module
MDL_BIN_TRG := swnb2 swnb3 disotest3 htrn2nb O2 O3 CFC11 CFC12 NO2 H2OH2O O2X clm # Binary targets
MDL_BIN_SYM_LNK := # Symbolic links
MDL_BIN_STB := ${MDL_BIN_TRG} ${MDL_BIN_SYM_LNK} # All executables in MY_BIN_DIR
MDL_BIN := $(addprefix ${MY_BIN_DIR}/,${MDL_BIN_STB}) # dst_cln removes these files

# Variables having to do with header files created by module
MDL_INC_TRG := # Include targets
MDL_INC_SYM_LNK := # Symbolic links
MDL_INC_STB := ${MDL_INC_TRG} ${MDL_INC_SYM_LNK} # All header files in MY_INC_DIR
MDL_INC := $(addprefix ${MY_INC_DIR}/,${MDL_INC_STB}) # dst_cln removes these files

# Files to install in ${DATA}/aca
DATA_RT_TRG = abs_xsx_H2OH2O_CFT98.txt abs_xsx_H2OH2O_Chy98.txt abs_xsx_H2OH2O_TGC98.txt abs_xsx_O2O2.txt abs_xsx_WMO85.txt arese_tst.txt mls_afgl_73lvl.txt mls_cld.txt mls_clr.txt mls_icrccm_35lvl.txt mls_icrccm_92lvl.txt mls_snw.txt mls_snw_FZR07.txt mls_snw_JRI.txt mls_snw_lgge.txt mlw_afgl_73lvl.txt prf_sas_65N.txt prf_sas_75N.txt prf_sas_Smm.txt sas_afgl_73lvl.txt saw_afgl_73lvl.txt saw_icrccm_18lvl.txt snw_noatm.txt snw_lgge.txt snw_lgge_13mm_32lvl.txt snw_lgge_25mm_32lvl.txt snw_lgge_hmg.txt snw_lgge_nhm.txt snw_crocus_08lvl.txt snw_nfn_01lvl.txt spc_KFB84.txt spc_Kur95_01wvn.txt spc_Kur95_20wvn.txt spc_Kur97.txt spc_NeL84.txt trp_afgl_73lvl.txt trp_icrccm_18lvl.txt trp_icrccm_35lvl.txt trp_icrccm_92lvl.txt
# DATA_RT targets
DATA_RT_SYM_LNK := # Symbolic links
DATA_RT_STB := ${DATA_RT_TRG} ${DATA_RT_SYM_LNK} # All files in DATA_RT
DATA_RT_FL := $(addprefix ${DATA_RT}/,${DATA_RT_STB}) # dst_cln removes these files

# Rule to install RT data file is copy from source directory to RT data directory
${DATA_RT}/% : %
	cp -f -p $(notdir $<) $@

# Files to install on website
C_FL := date_time.c # C files
DAT_FL := abs_xsx_H2OH2O.nc abs_xsx_NO2.nc abs_xsx_O2O2.nc abs_xsx_O2.nc abs_xsx_O3_WMO85.nc abs_xsx_O3_JPL15.nc abs_xsx_O3_HTR16.nc abs_xsx_CFC11_HTR16.nc abs_xsx_CFC12_HTR16.nc aer_H2SO4_PaW75.nc aer_dust_like.nc aer_h2so4_215K.nc aer_h2so4_300K.nc aer_meteoric_dust.nc aer_mineral_dust.nc aer_saharan_dust.nc aer_sulfate.nc aer_volcanic_dust.nc ice_05.nc ice_20.nc lqd_05.nc lqd_07.nc lqd_10.nc mls_cld.nc mls_cld.txt mls_clr.nc mls_clr.txt mls_icrccm_92lvl.nc mls_icrccm_92lvl.txt nst_FSBR.nc spc_Kur95_01wvn.nc spc_Kur95_20wvn.nc spc_LaN68.nc spc_NeL84.nc spc_ThD71.nc swnb_CO2.nc swnb_H2O.nc swnb_O2.nc swnb_OH.nc trp_icrccm_92lvl.nc trp_icrccm_92lvl.txt # Data files
DST_FL := swnb2.F90 clm.F90 disort2.f d1mach.f r1mach.f linpak.f bdref.f # Distribution files
DST_MFL := ChangeLog INSTALL Makefile README TODO VERSION # Distribution metafiles
#DST_MFL += FAQ index.shtml MANIFEST
F_FL := dbg_mdl shr_kind_mod erf_mdl gmm_mdl sng_mdl mmr_mdl nf90_utl prn_mdl utl_mdl xtr_mdl vec_mdl # Fortran files in libcsz_f90
F_FL += phys_cst_mdl cape_mdl wnd_PDF_mdl clm_mdl dmr_mdl fnd_cst_mdl htrn_mdl drv_cst_mdl flx_slr_mdl slr_crd_mdl wvl_mdl rt_mdl tdy_mdl # Fortran files in libcsm_f90
SH_FL := pvmgetarch # Shell files

# C and C++ rules
${MY_OBJ_DIR}/%.o : %.c
	${CC} ${CPPFLAGS} ${CFLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)
${MY_OBJ_DIR}/%.o : %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o ${MY_OBJ_DIR}/$(notdir $@)

# Fortran rules
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
 %.${MDL_SFX} : %.F90
	- if test -f ${NETCDF_LIB}/netcdf.${MDL_SFX}; then ln -s -f ${NETCDF_LIB}/netcdf.${MDL_SFX} .; else if test -f ${NETCDF_INC}/netcdf.${MDL_SFX} ; then ln -s -f ${NETCDF_INC}/netcdf.${MDL_SFX} .; fi; fi
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f90
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FFLAGS} ${CPPFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FFLAGS} -o ${MY_OBJ_DIR}/$(notdir $@) $<
%.f : %.F
	${CPP} -P ${CPPFLAGS} -o $@ $<
%.f90 : %.F90
	${CPP} -P ${CPPFLAGS} -o $@ $<

# Dependency rules
# Automatically generate a dependency file for each source file
# $* is the stem, e.g., f
# $@ is the filename of the target, e.g., f.${DPN_SFX}
DPN_GNR := makdep
${MY_DPN_DIR}/%.${DPN_SFX} : %.c
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.cc
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.F
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.F90
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.f
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@
${MY_DPN_DIR}/%.${DPN_SFX} : %.f90
	@echo "Building dependency file $@"
	@${DPN_GNR} -f ${CPP_PTH} -D ${MY_DPN_DIR} -O ${MY_OBJ_DIR} $< > $@

# Rules for installing header files
${MY_INC_DIR}/%.com : %.com
	cp -f -p $(notdir $@) $@
${MY_INC_DIR}/%.F90 : %.F90
	cp -f -p $(notdir $@) $@

# First LDFLAGS is for typical C programs with netCDF, math, and networking
# Second LDFLAGS enables C/Fortran linking

# Works on AIX and AIX46K
ifneq (${null},$(findstring AIX,${PVM_ARCH}))
ifeq (${OMP},Y)
 CC := xlc_r -qlanglvl=extc99
 FC := xlf95_r
# -qsmp=noauto : Turn on SMP/OMP code generation but do no automatic parallelization
# -qsmp=omp : Use industry standard OMP without IBM extensions
 OMP_FLG := -qsmp=omp
else
 CC := xlc
 CPP_DFN += -U_OPENMP
 FC := xlf95
endif # endif OMP
CXX := xlC
CFLAGS := ${OMP_FLG} -qmaxmem=-1
CPP := /lib/cpp -P
# AIX xlf requires passing all CPP tokens of type -DFOO and -UFOO with -WF flag
CPP_DFN := ${CPP_TKN}
CPP_DFN_F := -WF,$(subst $(space),$(comma),$(CPP_DFN))
CPP_DFN_F := $(subst $(comma)$(comma),$(comma),$(CPP_DFN_F))
CPPINC := ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
# -qarch=auto : Automatically detect architecture of compiling machine and assume execution on same machine
# -qsrcmsg prints transgressing source line with finger
# -qsuppress=1501-245 : Suppress RLIM_INFINITY memory message due to ulimits
# -qtune=auto : Optimize executable for architecture detected during compilation
FFLAGS := ${OMP_FLG} -qarch=auto -qmaxmem=-1 -qspillsize=2500 -qsuppress=1501-245 -qtune=auto
FIXEDFLAGS := -qsuffix=f=f:cpp=F -qfixed=132
FREEFLAGS := -qsuffix=f=f90:cpp=F90
LD := ld
# Heap is ~eight times larger than stack
# Heap can grow to ~2 GB with 32-bit addressing
# Local arrays or arrays whose size is known at compile time will be allocated from the heap
# Arrays whose size is not known at compile time will always be allocated from the stack
# -bmaxdata: Overrides ulimits stack and heap settings
# MATCH adds to LDFLAGS: -bmaxdata:2147483648 -bmaxstack:268435456
# Obtain maximum allowable: -bmaxdata:0x80000000 -bmaxstack:0x10000000
#LDFLAGS := ${OMP_FLG} -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS} -bnoquiet
LDFLAGS := ${OMP_FLG} -L${MY_LIB_DIR} ${NC_LDFLAGS} ${NC_LIBS} -bnoquiet -qsuppress=1501-245 
LEX := lex
LINT := lint
YACC := yacc
ifeq (${PRC},D)
# -qdpc=e : Promote all floating-point constants to double precision
 FFLAGS += -qrealsize=8 -qintsize=4 -qdpc=e
else
 FFLAGS += -qrealsize=4 -qintsize=4
endif
ifeq (${OPTS},O)
# -qstrict: Ensure that -O3 optimizations do not alter program semantics
# -Q : Inline all appropriate subprograms
 CFLAGS += -O3 -g -qstrict -Q
 FFLAGS += -O3 -g -qstrict -Q
endif
ifeq (${OPTS},D)
 CFLAGS += -g
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=inv:ov:zero:en -qfullpath -qsigtrap -qextchk -qwarn64
endif
ifeq (${OPTS},X)
 CFLAGS += -g
# -qcheck=all generates run-time bounds, divide-by-zero, and null-dereference checks
# -qthreaded: Ensure all internally-generated temporaries have automatic storage
# -qwarn64: Detect truncation of 8 byte integers to 4 bytes
# -qflttrap=... # Detect and trap floating point arithmetic problems
# :und: underflow # Underflow NB: Underflow is often OK and too hard to get rid of
# :imp: imprecise # Check for specified exceptions on subprogram entry and exit only
# :inv: invalid # Invalid operations
# :ov: overflow # Overflow
# :zero: zero # Divide-by-zero
# :en: enable # Enable specified exception trapping
# -qfullpath # Show full path back to floating point exceptions
# -qhalt=e stop compiler if error severity equals or exceeds i, w, e, s, u
 FFLAGS += -g -qinitauto=7FF7FFFF -qflttrap=inv:ov:und:zero:en -qfullpath -qsigtrap -qextchk -qthreaded -qwarn64 -C
endif
# -q64: Select 64-bit compiler mode (required for accessing large files)
# -qwarn64: Warn on possible long-to-integer or pointer-to-integer truncation
ifeq (${ABI},64)
 AR := ar -X 64
 CFLAGS += -q64
 FFLAGS += -q64
 LDFLAGS += -q64 
else
 CPPFLAGS += -D_LARGE_FILES
endif # endif ABI
CXXFLAGS := ${CFLAGS}
${MY_OBJ_DIR}/%.o : %.F90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f90
	${FC} -c ${FREEFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.F
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
${MY_OBJ_DIR}/%.o : %.f
	${FC} -c ${FIXEDFLAGS} ${FFLAGS} ${CPP_DFN_F} ${CPPINC} -o ${MY_OBJ_DIR}/$(notdir $@) $<
endif
# endif AIX

ifeq (${PVM_ARCH},ALPHA)
ifeq (${OMP},Y)
 OMP_FLG := -omp
endif # endif OMP
CXX := cxx -std strict_ansi
CC := cc
CFLAGS := ${OMP_FLG}
CPP := cpp
CPPFLAGS := -DOSF1 ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f95
FFLAGS := -c ${OMP_FLG} -automatic -fpe3
FIXEDFLAGS := -extend_source
FREEFLAGS := -DHIDE_SHR_MSG -free
LD := ld
LDFLAGS := ${OMP_FLG} -L${MY_LIB_DIR} -lcsm_f90 -lcsz_f90 ${NC_LDFLAGS} ${NC_LIBS}
LEX := lex
LINT := lint
YACC := yacc
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
else
 FFLAGS += -r4 -i4
endif # PRC != D
ifeq (${OPTS},O)
 CFLAGS += -O2 -g
 FFLAGS += -O3 -g -inline speed
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
CXX := g++
CC := cc
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f90
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsm_f90 -lcsz_f90 ${NC_LDFLAGS} ${NC_LIBS}
LDFLAGS += -L/lib -lf
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS := -O -g
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
# 1999/05/24 Added -C for cpp to retain comments (else cpp thinks Fortran string concatenator // is a C++ comment and removes it)
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
CXX := ${LINUX_CXX}
CC := ${LINUX_CC}
CPP := ${CC}
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := ${LINUX_FC}
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsm_f90 -lcsz_f90 ${NC_LDFLAGS} ${NC_LIBS}
LEX := flex
LINT := lint
YACC := bison
# GNU Compiler Collection
ifeq (${CC},gcc)
 CFLAGS := -Wall
ifeq (${OPTS},O)
 CFLAGS += -O -g ${GCC_RCH_ARG}
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
 CXXFLAGS := ${CFLAGS}
endif # end GNU Compiler Collection
# G77 Fortran compiler
ifeq (${FC},g77)
 FFLAGS := -ffixed-line-length-132
ifeq (${OPTS},O)
 FFLAGS += -O -g ${GCC_RCH_ARG}
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
endif # end G77 Fortran compiler
# G95 Compiler
ifeq (${FC},g95)
# FFLAGS := -ffixed-line-length-132 -Wline-truncation -Waliasing -Wimplicit-none -Wunused-label -Wall
# FFLAGS := -ffixed-line-length-132 -Wimplicit-none -Wunused-label -Wall
 FFLAGS := -ffixed-line-length-132 -Wimplicit-none -Wall -Wuninitialized
ifeq (${PRC},D)
 FFLAGS += -r8
else
 FFLAGS += 
endif
ifeq (${OPTS},O)
 FFLAGS += -O -g
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 FFLAGS += -g -O -fbounds-check -pedantic
# -finteger : Initialize uninitialized scalar integer variables
# -flogical=none/true/false : Initialize uninitialized scalar logical variables
# -fpointer=none/null/invalid : Initialize uninitialized scalar pointers
# -freal=none/zero/nan/inf/+inf/-inf : Initialize uninitialized scalar logical variables
 FFLAGS += -finteger=-2147483648 -flogical=none -freal=nan -fpointer=invalid
endif
endif # end G95 Compiler
# GNU Fortran Compiler (called as gfortran-mp-4.9 on MACOSX)
ifneq (${null},$(findstring gfortran,${FC}))
 FFLAGS := -ffixed-line-length-132 -Wall
ifeq (${OMP},Y)
 OMP_FLG_C := -fopenmp
 OMP_FLG_F := -fopenmp
 LDFLAGS += -lgomp -lpthread
endif # endif OMP
ifeq (${PRC},D)
# https://gcc.gnu.org/onlinedocs/gfortran/Fortran-Dialect-Options.html
# FFLAGS += -fdefault-real-8
 FFLAGS += -fdefault-real-8 -fdefault-double-8
else
 FFLAGS += 
endif
ifeq (${OPTS},O)
 FFLAGS += -O -g
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += -g
endif
ifeq (${OPTS},X)
# FFLAGS += -g -O -fbounds-check -pedantic -Wconversion -Werror -Wimplicit-interface -Wunderflow
 FFLAGS += -g -O -fbounds-check -pedantic -Wconversion -Wimplicit-interface -Wunderflow
endif
endif # end GNU Fortran Compiler
# Intel (Kai) C Compiler
ifeq (${CC},icc)
 CFLAGS := 
 OMP_FLG_C := -openmp
ifeq (${OPTS},O)
 CFLAGS += -O -g
endif
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},R)
 CFLAGS +=
endif
ifeq (${OPTS},X)
 CFLAGS += -g
endif
 CXXFLAGS := ${CFLAGS}
endif # end Intel (Kai) C Compiler
# Intel (Kai) Fortran Compiler
ifeq (${FC},ifort)
# -132: specify 132 column lines for fixed form sources (same as -extend_source)
# -axW: generate generic IA32 code and specialized code exclusively for Pentium 4 processors (requires linking to libsvml.a)
# -diag-disable 5268: quiet warnings that comments exceed allowed length
# -e95: die on warnings for non-standard Fortran95 (which includes OpenMP)
# -fast: enable -xP -O3 -ipo -static
# -fpconstant: extend precision of single precision constants assigned to double precision variables to double precision
# -fpp2: necessary, but not sufficient, for  OpenMP
# -ftz: enable/disable flush denormal results to zero
# -ip: enable single-file IP optimizations
# -ipo: enable multi-file IP optimizations
# -module: specifies path to place .mod files
# -mp: maintain floating point precision (disables some optimizations) (deprecated 20090129)
# -nus: do not append an underscore to external subroutine names
# -openmp: enable compiler to generate multi-threaded OpenMP code
# -vec_report1: indicate vectorized loops (DEFAULT)
# -tune pn4: optimize for Pentium(R) 4 processor (DEFAULT)
# -us: append an underscore to external subroutine names
# -vms: enable VMS and DEC statement extensions
# -w95: supress warnings for non-standard fortran (which includes OpenMP) (deprecated 20140129)
 FFLAGS := -diag-disable 5268 -extend_source -fpconstant -fpp2 -implicitnone
# -lsvml needed with -[a]xW
 LDFLAGS += -lsvml
 OMP_FLG_F := -openmp
ifeq (${PRC},D)
 FFLAGS += -i4 -r8
else
 FFLAGS += -i4
endif
ifeq (${OPTS},O)
 FFLAGS += -O3 -g -vec_report1 -ftz ${ICC_RCH_ARG}
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
# -C: extensive runtime diagnostics
# -e95 issues warnings for non-standard fortran (which includes OpenMP)
 FFLAGS += -g -C -e95
endif
ifeq (${MPI},Y)
 CC := $(subst icc,mpicc,${CC})
 CXX := $(subst icpc,mpiCC,${CXX})
 FC := $(subst ifort,mpif90,${FC})
endif # endif MPI
endif # end Intel (Kai) Fortran Compiler
# Lahey-Fujitsu Fortran Compiler
ifeq (${FC},lf95)
# --wide: Allow wide source format for fixed form
# NB: Convert to/from big-endian with runtime option: a.out -Wl,-T[unit_nbr]
 FFLAGS := -DLinux -DHIDE_SHR_MSG --mod .:${NETCDF_INC} --wide
 OMP_FLG_F := --openmp
ifeq (${PRC},D)
# Passing --dbl confuses resolution of overloaded float/double functions
# Ensure such functions are enclosed in #ifndef PRC_DBL
 FFLAGS += --dbl
else
 FFLAGS += 
endif
ifeq ($(STC),Y)
 LDFLAGS := --staticlink ${LDFLAGS}
endif # endif STC
ifeq (${OPTS},O)
# --tpp: Pentium Pro/Pentium II code generation
# --[n]sav: [Do not] save local variables
 FFLAGS += -O3 -g --tpp --nchk --ntrace --nsav --nwarn
endif
ifeq (${OPTS},D)
 FFLAGS += -g --chk a,e,s --pca --trace --trap
endif
ifeq (${OPTS},R)
 FFLAGS += 
endif
ifeq (${OPTS},X)
 FFLAGS += -g --chk a,e,s,u,x --chkglobal --pca --trace --trap
endif
endif # end Lahey-Fujitsu Fortran Compiler
# Pathscale (QLogic) Fortran Compiler
ifeq (pathf95,$(firstword ${FC}))
# man -k pathscale for full listing
# -col120 : Allow F77 source to be 120 columns
# -O2 = -O : Default optimization
# -Ofast = -O3:
# -woff801 : Turn off unknown compiler directive warnings triggered by !cdir
 FFLAGS := -extend_source -byteswapio -woff801
 LDFLAGS += 
 OMP_FLG_F := -apo -mp
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
else
 FFLAGS += -r4 -i4
endif
ifeq (${OPTS},O)
 FFLAGS += -O3 -g ${PSC_RCH_ARG}
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += -O2
endif
ifeq (${OPTS},X)
# -C: Bounds checking
 FFLAGS += -g -C -trapuv
endif
endif # end Pathscale (QLogic) Fortran Compiler
# Portland Group C++ Compiler
ifeq (${CXX},pgCC)
# Enable Large File Support (LFS) by default
 CFLAGS := -Mlfs
 LDFLAGS += -Mlfs
 OMP_FLG := -mp
ifeq (${OPTS},D)
 CFLAGS += -g
endif
ifeq (${OPTS},O)
 CFLAGS += -fast
endif
ifeq (${OPTS},R)
 CFLAGS += 
endif
ifeq (${OPTS},X)
 CFLAGS += -g -Mbounds
endif
 CXXFLAGS := ${CFLAGS}
endif # endif Portland Group C++ Compiler
# Portland Group Fortran Compiler
ifneq (${null},$(findstring pgf,${FC}))
# -Kieee : Used by CCSM
# -Ktrap=align|denorm|divz|fp|inexact|inv|none|ovf|unf : IEEE Trap conditions
# -Ktrap=fp : Equivalent to -Ktrap=inv,divz,ovf
# -M[no]daz : Treat denormalized numbers as zero
 FFLAGS := -Mextend -byteswapio -Mrecursive -Mdalign -Ktrap=fp -Mlfs
 OMP_FLG_F := -mp
ifeq (${PRC},D)
 FFLAGS += -r8 -i4
else
 FFLAGS += -r4 -i4
endif
ifeq ($(STC),Y)
 LDFLAGS := -Wl,-Bstatic ${LDFLAGS}
endif # endif STC
ifeq (${OPTS},O)
# -fast    Common optimizations: -O2 -Munroll=c:1 -Mnoframe -Mlre
# -fastsse == -fast -Mvect=sse -Mscalarsse -Mcache_align -Mflushzendif
# -mvect=nosse : Used by CCSM
# 20061105: -fastsse breaks CAM unexpectedly and is not worth it anyway
 FFLAGS += -g -fast -mvect=nosse
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 FFLAGS += -g -Mbounds
endif
endif # end Portland Group Fortran Compiler
# End compiler-specific switches for Linux
# Linux-generic processing for after compiler-specific argument have been set:
ifeq (${OMP},Y)
 CFLAGS += ${OMP_FLG_C}
 FFLAGS += ${OMP_FLG_F}
 LDFLAGS := ${OMP_FLG_F} ${LDFLAGS}
endif # endif OMP
ifeq (${MPI},Y)
 CC := $(subst icc,mpicc,${CC})
 CXX := $(subst icpc,mpiCC,${CXX})
 FC := $(subst ifort,mpif90,${FC})
endif # endif MPI
endif
# endif LINUX, LINUXALPHA, LINUXAMD64, LINUXARM, and FREEBSD and MACOSX

ifeq (${UNAMES},IRIX64)
CXX := CC
CC := cc -c99
CFLAGS := -64 -mips4
CPP := /usr/lib/cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} -I${NETCDF_LIB}
FC := f90 -cpp
FFLAGS := -64 -mips4 -extend_source -mp -mpio
LD := ld
# 20000706: Added -mp -mpio for OpenMP compliance
LDFLAGS := -64 -mp -mips4 -L${MY_LIB_DIR} -lrecipes_f ${NC_LDFLAGS} ${NC_LIBS}
# 20021220: Pat Worley on cam-dev says when using MPI on SGI, link library order is important 
#LDFLAGS += -mp -lmp -lmpi
LEX := flex
LINT := lint
YACC := bison
ifeq (${PRC},D)
 FFLAGS += -r8 -i4 
else
 FFLAGS += -r4 -i4
endif
ifeq (${OPTS},O)
 CFLAGS += -O2 -g3
 FFLAGS += -O2 -g3
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
 FFLAGS += -g -check_bounds -ansi -uvar
endif
endif
# endif IRIX64=SGI64,SGIMP64

ifeq (${UNAMES},SunOS) 
CXX := g++
CC := gcc -ansi
CFLAGS := -Wall
CPP := cpp
CPPFLAGS := ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC} 
FC := f90 -DHIDE_SHR_MSG
#FFLAGS := -xs -stackvar -e -Qoption f90comp -r8const
# Adding -dalign may be required if netCDF library was compiled with it
# -M: Module directory
FFLAGS := -xs -stackvar -e -M${NETCDF_INC}
LD := ld
LDFLAGS := -L/opt/SUNWspro/lib -lf77compat -L${MY_LIB_DIR} -lcsm_f90 -lcsz_f90 ${NC_LDFLAGS} ${NC_LIBS} -lsunmath -lsocket -lnsl
LDFLAGS += -lF77 -lM77 -lresolv
LEX := lex
LINT := lint
YACC := yacc
ifeq (${OPTS},O)
 CFLAGS += -O2 -g
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
else
 FFLAGS += -xtypemap=real:32,double:64,integer:32
ifeq (${FC},f77)
 FFLAGS += -r4 -i4
endif
endif
endif
# endif SunOS=SUN4SOL2,SUNMP

ifeq (${PVM_ARCH},WIN32)
CXX := ${LINUX_CXX}
CC := ${LINUX_CC}
# NB: nameser.h needs -Di386, but gcc sends -Di586 (on pentiums)
CPPFLAGS := -Di386 ${CPP_TKN} ${CPP_PTH} -I${NETCDF_INC}
FC := ${LINUX_FC}
LD := ld
LDFLAGS := -L${MY_LIB_DIR} -lcsm_f90 -lcsz_f90 ${NC_LDFLAGS} ${NC_LIBS}
LDFLAGS += 
LEX := flex
LINT := lint
YACC := bison
# G95 Compiler
ifeq (${FC},g95)
 FFLAGS := -ffixed-line-length-132 -Wimplicit-none -Wall -Wuninitialized
ifeq (${PRC},D)
 FFLAGS += -r8
else
 FFLAGS += 
endif
ifeq (${OPTS},O)
 FFLAGS += -O -g
endif
ifeq (${OPTS},D)
 FFLAGS += -g
endif
ifeq (${OPTS},R)
 FFLAGS += -g
endif
ifeq (${OPTS},X)
 FFLAGS += -g -O -fbounds-check -pedantic
# -finteger : Initialize uninitialized scalar integer variables
# -flogical=none/true/false : Initialize uninitialized scalar logical variables
# -fpointer=none/null/invalid : Initialize uninitialized scalar pointers
# -freal=none/zero/nan/inf/+inf/-inf : Initialize uninitialized scalar logical variables
# -ftrace=frame or -ftrace=full
# FFLAGS += -finteger=-2147483648 -flogical=none -freal=nan -fpointer=invalid
 FFLAGS += -finteger=-2147483648 -flogical=none -freal=nan -fpointer=invalid -ftrace=frame
endif
endif # end G95 Compiler
endif
# endif WIN32

# Link to Electric Fence library, see man efence
ifeq (${PVM_ARCH},LINUX)
ifneq (${null},$(findstring ${OPTS},DX))
# LDFLAGS += -lefence
endif
endif
# endif LINUX

# Disable OpenMP on platforms that automatically support it
ifeq (${OMP},N)
 FFLAGS := $(filter-out -omp,${FFLAGS})
 FIXEDFLAGS := $(filter-out -omp,${FIXEDFLAGS})
 CFLAGS := $(filter-out -mpio,${CFLAGS})
 FFLAGS := $(filter-out -mpio,${FFLAGS})
 LDFLAGS := $(filter-out -mpio,${LDFLAGS})
 CFLAGS := $(filter-out -mp,${CFLAGS})
 FFLAGS := $(filter-out -mp,${FFLAGS})
 LDFLAGS := $(filter-out -mp,${LDFLAGS})
 ifneq (${UNAMES},SunOS)
  CPPFLAGS += -U_OPENMP
 endif # endif SunOS
endif # endif OMP

# Define any remaining variables

# Default targets
all: ${MDL_BIN_TRG} lib
# .PHONY tells make to remake the following non-file targets 
.PHONY: all cln dst_cln debug ${MDL_BIN_TRG}
# Delete default suffixes---this should increase speed
.SUFFIXES: 
# Define suffixes that matter
.SUFFIXES: .cc .c .${DPN_SFX} .F .F90 .f90 .${MDL_SFX} .o
# Delete targets that were not successfully made
.DELETE_ON_ERROR:
lib_dir :
	- if ! test -d ${MY_LIB_DIR}; then mkdir -p ${MY_LIB_DIR}; fi
bin_dir :
	- if ! test -d ${MY_BIN_DIR}; then mkdir -p ${MY_BIN_DIR}; fi
obj_dir :
	- if ! test -d ${MY_OBJ_DIR}; then mkdir -p ${MY_OBJ_DIR}; fi

# Targets in bin
CLM_OBJ := ${MY_OBJ_DIR}/clm.o
CLM_OBJ += ${MY_OBJ_DIR}/date_time.o
clm:	${MY_BIN_DIR}/clm
${MY_BIN_DIR}/clm:	lib $(CLM_OBJ) 
	${FC} -o $@ $(CLM_OBJ) ${LDFLAGS}
	chmod 755 $@

HTRN2NB_OBJ := 	${MY_OBJ_DIR}/htrn2nb.o
htrn2nb: ${MY_BIN_DIR}/htrn2nb
${MY_BIN_DIR}/htrn2nb: lib $(HTRN2NB_OBJ) 
	${FC} -o $@ $(HTRN2NB_OBJ) ${LDFLAGS}
	chmod 755 $@

O2_OBJ := 	${MY_OBJ_DIR}/O2.o
O2:	${MY_BIN_DIR}/O2
${MY_BIN_DIR}/O2:	lib $(O2_OBJ) 
	${FC} -o $@ $(O2_OBJ) ${LDFLAGS} 
	chmod 755 $@

O3_OBJ := 	${MY_OBJ_DIR}/O3.o
O3:	${MY_BIN_DIR}/O3
${MY_BIN_DIR}/O3:	lib $(O3_OBJ) 
	${FC} -o $@ $(O3_OBJ) ${LDFLAGS} 
	chmod 755 $@

CFC11_OBJ := 	${MY_OBJ_DIR}/CFC11.o
CFC11:	${MY_BIN_DIR}/CFC11
${MY_BIN_DIR}/CFC11:	lib $(CFC11_OBJ) 
	${FC} -o $@ $(CFC11_OBJ) ${LDFLAGS} 
	chmod 755 $@

CFC12_OBJ := 	${MY_OBJ_DIR}/CFC12.o
CFC12:	${MY_BIN_DIR}/CFC12
${MY_BIN_DIR}/CFC12:	lib $(CFC12_OBJ) 
	${FC} -o $@ $(CFC12_OBJ) ${LDFLAGS} 
	chmod 755 $@

H2OH2O_OBJ := 	${MY_OBJ_DIR}/H2OH2O.o
H2OH2O:	${MY_BIN_DIR}/H2OH2O
${MY_BIN_DIR}/H2OH2O:	lib $(H2OH2O_OBJ) 
	${FC} -o $@ $(H2OH2O_OBJ) ${LDFLAGS}
	chmod 755 $@

MIE2NC_OBJ := ${MY_OBJ_DIR}/mie2nc.o
mie2nc:	${MY_BIN_DIR}/mie2nc
${MY_BIN_DIR}/mie2nc:	lib $(MIE2NC_OBJ) 
	${FC} -o $@ $(MIE2NC_OBJ) ${LDFLAGS} 
	chmod 755 $@

NBM_OBJ := ${MY_OBJ_DIR}/nbm.o ${MY_OBJ_DIR}/disort.o ${MY_OBJ_DIR}/d1mach.o ${MY_OBJ_DIR}/r1mach.o ${MY_OBJ_DIR}/linpak.o
nbm:	${MY_BIN_DIR}/nbm
${MY_BIN_DIR}/nbm:	lib $(NBM_OBJ) 
	${FC} -o $@ $(NBM_OBJ) ${LDFLAGS} 
	chmod 755 ${MY_BIN_DIR}/nbm

NO2_OBJ := 	${MY_OBJ_DIR}/NO2.o
NO2:	${MY_BIN_DIR}/NO2
${MY_BIN_DIR}/NO2:	lib $(NO2_OBJ) 
	${FC} -o $@ $(NO2_OBJ) ${LDFLAGS}
	chmod 755 $@

O2X_OBJ := 	${MY_OBJ_DIR}/O2X.o
O2X:	${MY_BIN_DIR}/O2X
${MY_BIN_DIR}/O2X:	lib $(O2X_OBJ) 
	${FC} -o $@ $(O2X_OBJ) ${LDFLAGS} 
	chmod 755 $@

SWNB_OBJ := ${MY_OBJ_DIR}/swnb.o ${MY_OBJ_DIR}/rt_utl.o ${MY_OBJ_DIR}/disort.o ${MY_OBJ_DIR}/d1mach.o ${MY_OBJ_DIR}/r1mach.o ${MY_OBJ_DIR}/linpak.o
swnb:	${MY_BIN_DIR}/swnb
${MY_BIN_DIR}/swnb:	lib $(SWNB_OBJ) 
	${FC} -o $@ $(SWNB_OBJ) ${LDFLAGS} 
	chmod 755 $@

SWNB2_OBJ := ${MY_OBJ_DIR}/swnb2.o ${MY_OBJ_DIR}/disort2.o ${MY_OBJ_DIR}/d1mach.o ${MY_OBJ_DIR}/r1mach.o ${MY_OBJ_DIR}/linpak.o ${MY_OBJ_DIR}/bdref.o
swnb2:	${MY_BIN_DIR}/swnb2
${MY_BIN_DIR}/swnb2:	lib $(SWNB2_OBJ) 
	${FC} -o $@ $(SWNB2_OBJ) ${LDFLAGS}
	chmod 755 $@

#SWNB3_OBJ := ${MY_OBJ_DIR}/swnb3.o ${MY_OBJ_DIR}/disort3.o ${MY_OBJ_DIR}/rdi1mach3.o ${MY_OBJ_DIR}/lapack3.o ${MY_OBJ_DIR}/linpak3.o ${MY_OBJ_DIR}/bdref3.o ${MY_OBJ_DIR}/errpack3.o
#swnb3:	${MY_BIN_DIR}/swnb3
#${MY_BIN_DIR}/swnb3:	lib $(SWNB3_OBJ) 
#	${FC} -o $@ $(SWNB3_OBJ) ${LDFLAGS}
#	chmod 755 $@

#DISOTEST3_OBJ := ${MY_OBJ_DIR}/disotest3.o ${MY_OBJ_DIR}/disort3.o ${MY_OBJ_DIR}/rdi1mach3.o ${MY_OBJ_DIR}/lapack3.o ${MY_OBJ_DIR}/linpak3.o ${MY_OBJ_DIR}/bdref3.o ${MY_OBJ_DIR}/errpack3.o ${MY_OBJ_DIR}/disobrdf3.o
#disotest3:	${MY_BIN_DIR}/disotest3
#${MY_BIN_DIR}/disotest3:	lib $(DISOTEST3_OBJ) 
#	${FC} -o $@ $(DISOTEST3_OBJ) ${LDFLAGS}
#	chmod 755 $@

# Target-specific variable values syntax TARGET ... : VARIABLE-ASSIGNMENT
# Rules begin in leftmost column else interpreted as commands

bin : ${MDL_BIN_TRG}
bin_cln: 
	rm -f ${MDL_BIN}
# Targets in bld
buildclean: bld_cln
bld_cln :
	cd ${MY_BLD_DIR}; rm -f TAGS
clean : cln
cln : dat_cln dpn_cln lib_cln mod_cln obj_cln
debug: dbg
dbg:
	@printf "ABI = ${ABI}\n"
	@printf "AR = ${AR}\n"
	@printf "CC = ${CC}\n"
	@printf "CFLAGS = ${CFLAGS}\n"
	@printf "CPP = ${CPP}\n"
	@printf "CPPFLAGS = ${CPPFLAGS}\n"
	@printf "CPPINC = ${CPPINC}\n"
	@printf "CPP_DFN = ${CPP_DFN}\n"
	@printf "CPP_DFN_F = ${CPP_DFN_F}\n"
	@printf "CPP_PTH = ${CPP_PTH}\n"
	@printf "CPP_TKN = ${CPP_TKN}\n"
	@printf "CXX = ${CXX}\n"
	@printf "CXXFLAGS = ${CXXFLAGS}\n"
	@printf "DBG = ${DBG}\n"
	@printf "FC = ${FC}\n"
	@printf "FFLAGS = ${FFLAGS}\n"
	@printf "FIXEDFLAGS = ${FIXEDFLAGS}\n"
	@printf "FREEFLAGS = ${FREEFLAGS}\n"
	@printf "GCC_RCH_ARG = ${GCC_RCH_ARG}\n"
	@printf "GSL_LIB = ${GSL_LIB}\n"
	@printf "HOST = ${HOST}\n"
	@printf "ICC_RCH_ARG = ${ICC_RCH_ARG}\n"
	@printf "HOSTNAME = ${HOSTNAME}\n"
	@printf "LDFLAGS = ${LDFLAGS}\n"
	@printf "LIBNETCDFF = ${LIBNETCDFF}\n"
	@printf "LINUX_CC = ${LINUX_CC}\n"
	@printf "LINUX_CXX = ${LINUX_CXX}\n"
	@printf "LINUX_FC = ${LINUX_FC}\n"
	@printf "MDL_DPN = ${MDL_DPN}\n"
	@printf "MDL_OBJ = ${MDL_OBJ}\n"
	@printf "MDL_PTH = ${MDL_PTH}\n"
	@printf "MDL_SRC = ${MDL_SRC}\n"
	@printf "MPI = ${MPI}\n"
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
	@printf "NC_LDFLAGS = ${NC_LDFLAGS}\n"
	@printf "NC_LIBS = ${NC_LIBS}\n"
	@printf "NETCDF4 = ${NETCDF4}\n"
	@printf "NETCDF4_ROOT = ${NETCDF4_ROOT}\n"
	@printf "NETCDF_INC = ${NETCDF_INC}\n"
	@printf "NETCDF_LIB = ${NETCDF_LIB}\n"
	@printf "NETCDF_ROOT = ${NETCDF_ROOT}\n"
	@printf "NETCDFF = ${NETCDFF}\n"
	@printf "OMP = ${OMP}\n"
	@printf "OMP_FLG_C = ${OMP_FLG_C}\n"
	@printf "OMP_FLG_F = ${OMP_FLG_F}\n"
	@printf "OPTS = ${OPTS}\n"
	@printf "PGI_RCH_ARG = ${PGI_RCH_ARG}\n"
	@printf "PRC = ${PRC}\n"
	@printf "PSC_RCH_ARG = ${PSC_RCH_ARG}\n"
	@printf "PVM_ARCH = ${PVM_ARCH}\n"
	@printf "UNAMES = ${UNAMES}\n"
	@printf "VPATH = ${VPATH}\n"
distclean : dst_cln
dst_cln : bld_cln cln lib_cln
install : inc 
tags : 
	etags *.cc *.c *.com *.F *.F90 *.f90 *.hh *.h Makefile *.txt
uninstall : bin_cln dst_cln inc_cln lib_cln
# Targets in dat
dat_cln :
	-cd ${MY_DAT_DIR}; rm -f foo.nc
data_rt : ${DATA_RT_FL}
data_rt_cln :
	rm -f ${DATA_RT_FL}

# Targets in dpn
dpn : ${MDL_DPN}
dpn_cln :
	rm -f ${MDL_DPN}
# Targets in dst
dst : ${C_FL} ${DAT_FL} ${DST_FL} ${DST_MFL} ${F_FL} ${SH_FL}
	tar cvzf ${DATA}/tmp/${MDL_NM}-${VERSION}.tar.gz ${C_FL} ${DAT_FL} ${DST_FL} ${DST_MFL} ${F_FL} ${SH_FL}
	scp -p ${DST_FL} ${DATA}/tmp/${MDL_NM}-${VERSION}.tar.gz ${MY_HTTP_MCH}:${MY_HTTP_DIR}
	-ssh ${MY_HTTP_MCH} "cd ${MY_HTTP_DIR};ln -s -f ${MDL_NM}-${VERSION}.tar.gz ${MDL_NM}.tar.gz"
# Targets in inc
inc : ${MDL_INC}
inc_cln :
	rm -f ${MDL_INC}
# Targets in lib
lib :
	-cd ${HOME}/sw/f; ${MAKE} lib
lib_cln :
	-cd ${MY_LIB_DIR};
# Targets in mod
mod : ${MDL_OBJ}
	mv -f $(wildcard ${MY_OBJ_DIR}/*.${MDL_SFX}) .
mod_cln :
	rm -f $(wildcard *.${MDL_SFX}) $(wildcard ${MY_OBJ_DIR}/*.${MDL_SFX})
#	-mv -f $(wildcard *.${MDL_SFX}) ${MY_OBJ_DIR}
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
# Incorporate dependency files into Makefile rules
-include ${MDL_DPN}
endif
