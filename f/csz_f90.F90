! $Id$ -*-f90-*-

! Purpose: Module to contain entire csz_f90 library

! Usage:
!use csz_f90 ! [mdl] Fortran9X utility library

module csz_f90 ! [mdl] Fortran9X utility library
  ! Commented-out files are not modules yet
!  use csz_mdl
  use nf90_utl
!  use prn_mdl
!  use sng_mdl
!  use tdy_mdl
!  use utl_mdl
!  use vec_mdl
  use xtr_mdl
  implicit none
contains
  ! Use C-preprocessor #include, not Fortran include, if any included files contain C-preprocessor commands
  ! Warning: Including modules in "contains" section is a syntax error
  ! Both these rules makes sense if you think about them...
!#include "csz_mdl.F90"
!#include "prn_mdl.F90"
!#include "sng_mdl.F90"
!#include "tdy_mdl.F90"
!#include "utl_mdl.F90"
!#include "vec_mdl.F90"
end module csz_f90
