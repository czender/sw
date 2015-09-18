! cd ~/f;f90 -o tst tst.F90
! cd ~/f;g95 -o tst tst.F90
! cd ~/f;gfortran -o tst tst.F90
! cd ~/f;ifc -o tst tst.F90
! cd ~/f;lf95 -o tst tst.F90
! cd ~/f;pgf90 -o tst tst.F90
! cd ~/f;xlf95_r -c -qsuffix=f=f90:cpp=F90 -o tst tst.F90
! scp ~/f/tst.F90 dust.ess.uci.edu:f/tst.F90
program tst
  use tst_mdl
  integer::arg_idx
  integer::int_foo
  character(len=200)::arg_val
  real(r8) flt_foo,val
  arg_idx=1
  call getarg(arg_idx,arg_val)
  flt_foo=exp(1.0)
  int_foo=flt_foo
  val=fnc_foo(flt_foo)
  print *,'int_foo = ',int_foo
  print *,'arg_idx = ',arg_idx,'arg_val = ',arg_val
  print *,'val = ',val
  write (6,*) "Hello, World!"
  ! csz 20061125 following line causes ICE on gfortran 4.1.2:
  write (6,"(a,es15.8)") "2.0**(-0.0) = ",2.0**(-0.0)
end program tst
