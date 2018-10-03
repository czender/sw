! $Id$

! Purpose: Print formatting & utilities

! Copyright (C) 1994--2018 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
! use prn_mdl ! [mdl] Print formatting & utilities

module prn_mdl ! [mdl] Print formatting & utilities
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public dat_prn_f77 ! [sbr] Print block data Fortran77-style fixed format
  public dat_prn_f90 ! [sbr] Print block data Fortran90-style free format
  
contains 
  
  subroutine dat_prn_f77( & ! [sbr] Print block data Fortran77-style fixed format
       dat_nbr, & ! I [nbr] Number of elements in array
       dat_val, & ! I [frc] Values to print
       dat_nbr_per_ln, & ! I [nbr] Number of data per line
       dat_sng) ! I [sng] Variable name
    ! Purpose: Print block data to stderr in Fortran77-style fixed format
    ! Resulting code may then be included in Fortran77 code
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Input
    character(len=*),intent(in)::dat_sng ! I [sng] Variable name
    integer,intent(in)::dat_nbr ! I [nbr] Number of elements in array
    integer,intent(in)::dat_nbr_per_ln ! I [nbr] Number of data per line
    real,intent(in)::dat_val(dat_nbr) ! I [frc] Values to print
    ! Local
    integer dat_idx ! [idx] Counting index for data
    write (unit=0,fmt='(a11,a,a2)',advance="no") "      data ",dat_sng(1:ftn_strlen(dat_sng))," /"
    do dat_idx=1,dat_nbr
       if (mod(dat_idx,dat_nbr_per_ln) == 0) write (unit=0,fmt='(/,a11)',advance="no") "     $     "
       write (unit=0,fmt='(es14.7)',advance="no") dat_val(dat_idx)
       if (dat_idx /= dat_nbr) write (unit=0,fmt='(a2)',advance="no") ", "
    enddo ! end loop over dat_idx
    write (unit=0,fmt='(a2,/)') " /"
  end subroutine dat_prn_f77
  
  subroutine dat_prn_f90( & ! [sbr] Print block data Fortran90-style free format
       dat_nbr, & ! I [nbr] Number of elements in array
       dat_val, & ! I [frc] Values to print
       dat_nbr_per_ln, & ! I [nbr] Number of data per line
       dat_sng) ! I [sng] Variable name
    ! Purpose: Print input data to stderr in Fortran90-style fixed format
    ! Resulting code may then be included in Fortran90 code
    use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
    implicit none
    ! Input
    character(len=*),intent(in)::dat_sng ! I [sng] Variable name
    integer,intent(in)::dat_nbr ! I [nbr] Number of elements in array
    integer,intent(in)::dat_nbr_per_ln ! I [nbr] Number of data per line
    real,intent(in)::dat_val(dat_nbr) ! I [frc] Values to print
    ! Local
    character(8)::fmt_sng ! [sng] Format string
    integer dat_idx ! [idx] Counting index for data
    integer dgt_nbr ! [nbr] Number of digits
    dgt_nbr=ceiling(log10(real(dat_nbr)))
    ! Fix special cases where rounding up with ceiling() does not help
    if (mod(dat_nbr,10)==0) dgt_nbr=dgt_nbr+1
    fmt_sng='(a,iX,a)' ! [sng] Format string
    write (fmt_sng(5:5),'(i1)') dgt_nbr ! [sng] Format string
    write (unit=0,fmt=fmt_sng,advance="yes") "      integer,parameter::foo=",dat_nbr," ! [nbr] Number of foo"
    write (unit=0,fmt='(a,a,a5,/,a11)',advance="no") "      real(r8),dimension(foo),parameter::", &
         dat_sng(1:ftn_strlen(dat_sng)),"=(/ &","           "
    do dat_idx=1,dat_nbr
       write (unit=0,fmt='(es14.7)',advance="no") dat_val(dat_idx)
       if (dat_idx /= dat_nbr) write (unit=0,fmt='(a2)',advance="no") ", "
       if (mod(dat_idx,dat_nbr_per_ln) == 0) write (unit=0,fmt='(a1,/,a11)',advance="no") "&","           "
    enddo ! end loop over dat_idx
    write (unit=0,fmt='(a3,/)') " /)"
  end subroutine dat_prn_f90
  
end module prn_mdl ! [mdl] Print formatting & utilities
