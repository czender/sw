! Companion module to tst.F90 for debugging Fortran90 code/compilers
! /bin/rm tst_mdl.o tst.o
t ! export FFLAGS='-Wimplicit-none -Wall -O -g'
! cd ~/f;f90 ${FFLAGS} -c -o tst_mdl.o tst_mdl.F90;f90 ${FFLAGS} -c -o tst.o tst.F90;f90 -o tst tst.o tst_mdl.o
! cd ~/f;g95 ${FFLAGS} -c -o tst_mdl.o tst_mdl.F90;g95 ${FFLAGS} -c -o tst.o tst.F90;g95 -o tst tst.o tst_mdl.o
! cd ~/f;gfortran ${FFLAGS} -c -o tst_mdl.o tst_mdl.F90;gfortran ${FFLAGS} -c -o tst.o tst.F90;gfortran -o tst tst.o tst_mdl.o
! cd ~/f;ifc ${FFLAGS} -c -o tst_mdl.o tst_mdl.F90;ifc ${FFLAGS} -c -o tst.o tst.F90;ifc -o tst tst.o tst_mdl.o
! cd ~/f;lf95 ${FFLAGS} -c -o tst_mdl.o tst_mdl.F90;lf95 ${FFLAGS} -c -o tst.o tst.F90;lf95 -o tst tst.o tst_mdl.o
! cd ~/f;pgf90 ${FFLAGS} -c -o tst_mdl.o tst_mdl.F90;pgf90 ${FFLAGS} -c -o tst.o tst.F90;pgf90 -o tst tst.o tst_mdl.o
! cd ~/f;xlf95_r ${FFLAGS} -c -qsuffix=f=f90:cpp=F90 -o tst_mdl.o tst_mdl.F90
! scp ~/f/tst.F90 ~/f/tst_mdl.F90 sand.ess.uci.edu:f
module tst_mdl
  implicit none
  integer,parameter::r8=selected_real_kind(p=12)

contains

  real(r8) function fnc_foo(flt_foo)
    real(r8) flt_foo
    fnc_foo=log(flt_foo)
  end function fnc_foo

  integer function fnc_1(arg)
    integer arg
    fnc_1=arg**2
  end function fnc_1
  
  integer function fnc_2(arg)
    integer arg
    fnc_2=2*arg
  end function fnc_2
  
end module tst_mdl
