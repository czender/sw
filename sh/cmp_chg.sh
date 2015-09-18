#!/bin/sh

# $Id$

# Purpose: Change ABI of default system C++/Fortran libraries and modules 
# Currently supported ABI's are: g++, gfortran, icc, ifort, lf95, pgcc, pgf95
# User priveleges required to change C++ compiler (currently only user libraries depend on CXX)
# Root priveleges required to change Fortran compiler (system libraries like libnetcdf.a depend on FC)

# Usage: 
# Recommend changing fortran libraries with two commands:
# First changes environment variables in user's shell
# Second moves system libraries
# source ${HOME}/sh/cmp_chg.sh gcc g77;sudo ${HOME}/sh/cmp_chg.sh gcc g77
# source ${HOME}/sh/cmp_chg.sh gcc g95;sudo ${HOME}/sh/cmp_chg.sh gcc g95
# source ${HOME}/sh/cmp_chg.sh gcc gfortran;sudo ${HOME}/sh/cmp_chg.sh gcc gfortran
# source ${HOME}/sh/cmp_chg.sh gcc ifort;sudo ${HOME}/sh/cmp_chg.sh gcc ifort
# source ${HOME}/sh/cmp_chg.sh gcc lf95;sudo ${HOME}/sh/cmp_chg.sh gcc lf95
# source ${HOME}/sh/cmp_chg.sh gcc pgf95;sudo ${HOME}/sh/cmp_chg.sh gcc pgf95
# source ${HOME}/sh/cmp_chg.sh icc ifort;sudo ${HOME}/sh/cmp_chg.sh icc ifort
# source ${HOME}/sh/cmp_chg.sh icc g95;sudo ${HOME}/sh/cmp_chg.sh icc g95
# source ${HOME}/sh/cmp_chg.sh pgcc pgf95;sudo ${HOME}/sh/cmp_chg.sh pgcc pgf95

# cmp_c='gcc';cmp_f='gfortran';cmp_hyb="${cmp_c}-${cmp_f}";export LINUX_FC=${cmp_f}
# cmp_c='gcc';cmp_f='g95';cmp_hyb="${cmp_c}-${cmp_f}";export LINUX_FC=${cmp_f}
# cmp_c='gcc';cmp_f='pgf95';cmp_hyb="${cmp_c}-${cmp_f}";export LINUX_FC=${cmp_f}

# For simplicity, C compiler name determines C++ compiler name, e.g., gcc implies gcc for C code and g++ for C++ code, icc implies icc and icpc, etc.
cmp_c=$1; # C compiler name
cmp_f=$2; # Fortran compiler name
cmp_hyb="${cmp_c}-${cmp_f}";

# Change Fortran compiler
if [ "${cmp_f}" = 'g77' -o "${cmp_f}" = 'g95' -o "${cmp_f}" = 'gfortran' -o "${cmp_f}" = 'ifort' -o "${cmp_f}" = 'lf95' -o "${cmp_f}" = 'pgf95' ] ; then
# Set environment variable used by Makefiles
    case "${cmp_f}" in 
	g77* ) LINUX_FC='g77'; ;; # endif g77
	g95* ) LINUX_FC='g95'; ;; # endif g95
	gfortran* ) LINUX_FC='gfortran'; ;; # endif gfortran
	ifort* ) LINUX_FC='ifort'; ;; # endif ifort
	lf95* ) LINUX_FC='lf95'; ;; # endif lf95
	pgf95* ) LINUX_FC='pgf95'; ;; # endif pgf95
	* ) printf "ERROR: Compiler ${cmp_f} is unknown in cmp_chg.sh\n"; ;; # endif default
    esac # endcase ${cmp_f}
    export LINUX_FC # [cmd] Fortran compiler used by Makefiles
# Redirect personal libraries
    cd ${MY_LIB_DIR}
    printf "Changing pure Fortran libraries to ${cmp_f} ABI...\n"
    for lib in librecipes_f.a librecipes_f90.a libspecfun.a; do
	/bin/rm -f ${lib}
	ln -s -f ${lib}.${cmp_f} ${lib}
    done # end loop over lib
# Redirect system libraries (requires sudo priveleges)
    cd /usr/local/lib
# libnetcdf.a contains fortran object code too
    printf "Changing hybrid libraries to ${cmp_hyb} ABI...\n"
    sudo /bin/rm -f libnetcdf.a
    sudo ln -s -f libnetcdf.a.${cmp_hyb} libnetcdf.a
    cd /usr/local/include
    sudo /bin/rm -f netcdf.mod typesizes.mod
    sudo ln -s -f netcdf.mod.${cmp_f} netcdf.mod
    sudo ln -s -f typesizes.mod.${cmp_f} typesizes.mod
fi # !(g77,g95,gfortran,ifort,lf95,pgf95)

# C++ compiler
if [ "${cmp_c}" = 'gcc' -o "${cmp_c}" = 'icc' ] ; then
# Set environment variable used by Makefiles
    case "${cmp_c}" in 
	como* ) LINUX_CXX='como'; LINUX_CC='como --c99'; ;; # endif como
	gcc* ) LINUX_CXX='g++'; LINUX_CC='gcc -std=c99 -pedantic -D_BSD_SOURCE -D_POSIX_SOURCE'; ;; # endif g++
	icc* ) LINUX_CXX='icpc'; LINUX_CC='icc -std=c99 -D_BSD_SOURCE -D_POSIX_SOURCE'; ;; # endif icc
	insure* ) LINUX_CXX='insure'; LINUX_CC='insure'; ;; # endif insure
	* ) printf "ERROR: Compiler ${cmp_c} is unknown in cmp_chg.sh\n"; ;; # endif default
    esac # endcase ${cmp_c}
    export LINUX_CC # [cmd] C compiler used by Makefiles
    export LINUX_CXX # [cmd] C++ compiler used by Makefiles
    export MPICH_CC="${LINUX_CC} -I/usr/lib/mpich/include"
    export MPICH_CCC="${LINUX_CXX} -I/usr/lib/mpich/include"
# Redirect personal libraries
    printf "Changing C++ libraries to ${cmp_c} ABI..."
    cd ${HOME}/c++; make obj_cln lib_cln
    cd ${MY_LIB_DIR}
    for lib in libcsm_c++.a libcsz_c++.a libnco_c++.a; do
	/bin/rm -f ${lib}
	ln -s -f ${lib}.${cmp_c} ${lib}
    done # end loop over lib
fi # !(g++,icc)

# IITAF
printf "done\n"
