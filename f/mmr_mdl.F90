! $Id$

! Purpose: Memory management

! Copyright (C) 1994--2014 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! Usage:
! scp ~/f/mmr_mdl.F90 dataproc.ucar.edu:f
! pgf90 -c -Mextend -Mnosecond_underscore -byteswapio -Mrecursive -Mdalign -Ktrap=fp -Mi4 -fast -mp -DPRC_FLT  -DLINUX -I. -I/home/zender/include -I/usr/local/include -o /home/zender/obj/LINUX/mmr_mdl.o mmr_mdl.F90
! xlf95_r -c -qsuffix=f=f90:cpp=F90 -o mmr_mdl mmr_mdl.F90

! Notes on module
! Pointers may not be declared with intent attribute, MeR96 p. 90, NyL97 p. 92
! Functions which return pointers can use result() clause NyL97 p. 92

! Notes on PRC_DBL, auto-promotion, overloading:
! Assume PRC_DBL is defined when compiler is passed autopromotion flags ("-r8")
! This means all floats are promoted to double precision
! As a result, ftn_arg_get_flt() is identical in precision to ftn_arg_get_dbl()
! Thus compiler cannot disambiguate overloading these two functions
! Result is OK on all known compilers except Lahey lf95 which complains
! Solution adopted is do not define ftn_arg_get_flt() if auto-promotion is enabled
! Alternate solution is never use auto-promotion with lf95

! Usage:
! use mmr_mdl ! [mdl] Memory management

module mmr_mdl ! [mdl] Memory management
  use sng_mdl,only:ftn_strlen ! [mdl] String manipulation
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public wrp_allocate ! [sbr] Overloaded memory allocation with error checking
  public wrp_deallocate ! [sbr] Overloaded memory de-allocation with error checking
  public mmr_tst ! [sbr] Test memory management
  
  ! Overloaded memory allocation routines
  interface wrp_allocate
     module procedure wrp_allocate_scl_FourByteInt,wrp_allocate_1D_FourByteInt, &
          wrp_allocate_scl_EightByteInt,wrp_allocate_1D_EightByteInt, &
#ifndef PRC_DBL 
          wrp_allocate_scl_FourByteReal,wrp_allocate_1D_FourByteReal, &
#endif /* PRC_DBL */
          wrp_allocate_scl_EightByteReal,wrp_allocate_1D_EightByteReal
  end interface ! wrp_allocate()

  ! Overloaded memory deallocation routines
  interface wrp_deallocate
     module procedure wrp_deallocate_FourByteInt,wrp_deallocate_EightByteInt, &
#ifndef PRC_DBL 
          wrp_deallocate_FourByteReal, &
#endif /* PRC_DBL */
          wrp_deallocate_EightByteReal
  end interface ! wrp_deallocate()

contains 

  subroutine mmr_tst ! [sbr] Test memory management
    integer::lmn_nbr
    integer::idx
    integer::rcd
    integer,pointer::var_ptr(:)
    integer,dimension(1)::dmn_cnt=(/10/)
    character(7)::var_nm="var_ptr"
    lmn_nbr=10
!    var_ptr=>wrp_allocate_fnc(lmn_nbr)
    call wrp_allocate(dmn_cnt,var_ptr,var_nm)
    do idx=1,lmn_nbr
       var_ptr(idx)=idx
       write (6,'(a7,i3,a2,i3)') 'var_ptr(',idx,')=',var_ptr(idx)
    enddo
    call wrp_deallocate(var_ptr,rcd,var_nm)
  end subroutine mmr_tst
  
  ! Lahey lf95 prints a bogus warning here
  ! Sun f90 prints:
  ! "function wrp_allocate_fnc(lmn_nbr) result(var_ptr) ! [sbr] Allocate 1D array"
  ! "mmr_mdl.F90", Line = 77, Column = 45: WARNING: The result of function name "VAR_PTR" in the function subprogram is not defined.
  function wrp_allocate_fnc(lmn_nbr) result(var_ptr) ! [sbr] Allocate 1D array
    integer,dimension(:),pointer::var_ptr
    integer,intent(in)::lmn_nbr
    integer::rcd
    allocate(var_ptr(lmn_nbr),stat=rcd)
    if(rcd /= 0) stop "allocate() failed in wrp_allocate_fnc"
  end function wrp_allocate_fnc

  subroutine wrp_allocate_1D_FourByteInt( & ! [sbr] Allocate 1D array
       dmn_cnt, & ! I [nbr] Array of dimension sizes
       var_ptr, & ! O [ptr] Pointer to values
       var_nm) ! I [sng] Variable name
    ! Full template for memory allocation
    ! All overloaded subroutines based on this structure but without comments
    ! Usage:
    implicit none
    ! Parameters
    character(len=*),parameter::sbr_nm='wrp_allocate_1D_FourByteInt' ! [sng] Subroutine name
    ! Input
    integer,intent(in)::dmn_cnt(:) ! I [nbr] Array of dimension sizes
    character(len=*),optional,intent(in)::var_nm ! I [sng] Variable name
    ! Output
    integer(selected_int_kind(6)),pointer::var_ptr(:) ! O [ptr] Pointer to values
    ! Local
    integer::rcd ! [rcd] Return success code
    ! Main code
    ! Check if variable is already associated
    if (associated(var_ptr)) then
       if(present(var_nm)) then 
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable ',var_nm,' already associated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable already associated'
       endif ! endif
       stop
    endif ! endif associated
    ! Allocate memory
    allocate(var_ptr(dmn_cnt(1)),stat=rcd) ! Coordinate variable
    ! Handle errors
    if (rcd /= 0) then  
       if(present(var_nm)) then 
          write (6,'(3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating ',var_nm
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating memory'
       endif ! endif
       write (6,'(2a,i3)') sbr_nm(1:ftn_strlen(sbr_nm)),': rcd = ',rcd
       stop  
    endif ! endif err  
    return
  end subroutine wrp_allocate_1D_FourByteInt

  subroutine wrp_allocate_1D_EightByteInt(dmn_cnt,var_ptr,var_nm)
    ! Skeleton template based on wrp_allocate_1D_FourByteInt()
    implicit none
    character(len=*),parameter::sbr_nm='wrp_allocate_1D_EightByteInt'
    integer,intent(in)::dmn_cnt(:)
    character(len=*),optional,intent(in)::var_nm
    integer(selected_int_kind(12)),pointer::var_ptr(:)
    integer::rcd
    if (associated(var_ptr)) then
       if(present(var_nm)) then 
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable ',var_nm,' already associated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable already associated'
       endif
       stop
    endif
    allocate(var_ptr(dmn_cnt(1)),stat=rcd)
    if (rcd /= 0) then  
       if(present(var_nm)) then 
          write (6,'(3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating ',var_nm
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating memory'
       endif
       write (6,'(2a,i3)') sbr_nm(1:ftn_strlen(sbr_nm)),': rcd = ',rcd
       stop  
    endif
    return
  end subroutine wrp_allocate_1D_EightByteInt

#ifndef PRC_DBL 
  subroutine wrp_allocate_1D_FourByteReal(dmn_cnt,var_ptr,var_nm)
    implicit none
    character(len=*),parameter::sbr_nm='wrp_allocate_1D_FourByteReal'
    integer,intent(in)::dmn_cnt(:)
    character(len=*),optional,intent(in)::var_nm
    real(selected_real_kind(p=6)),pointer::var_ptr(:)
    integer::rcd
    if (associated(var_ptr)) then
       if(present(var_nm)) then 
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable ',var_nm,' already associated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable already associated'
       endif
       stop
    endif
    allocate(var_ptr(dmn_cnt(1)),stat=rcd)
    if (rcd /= 0) then  
       if(present(var_nm)) then 
          write (6,'(3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating ',var_nm
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating memory'
       endif
       write (6,'(2a,i3)') sbr_nm(1:ftn_strlen(sbr_nm)),': rcd = ',rcd
       stop  
    endif
    return
  end subroutine wrp_allocate_1D_FourByteReal
#endif /* PRC_DBL */

  subroutine wrp_allocate_1D_EightByteReal(dmn_cnt,var_ptr,var_nm)
    implicit none
    character(len=*),parameter::sbr_nm='wrp_allocate_1D_EightByteReal'
    integer,intent(in)::dmn_cnt(:)
    character(len=*),optional,intent(in)::var_nm
    real(selected_real_kind(p=12)),pointer::var_ptr(:)
    integer::rcd
    if (associated(var_ptr)) then
       if(present(var_nm)) then 
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable ',var_nm,' already associated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable already associated'
       endif
       stop
    endif
    allocate(var_ptr(dmn_cnt(1)),stat=rcd) ! Coordinate variable
    if (rcd /= 0) then  
       if(present(var_nm)) then 
          write (6,'(3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating ',var_nm
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating memory'
       endif
       write (6,'(2a,i3)') sbr_nm(1:ftn_strlen(sbr_nm)),': rcd = ',rcd
       stop  
    endif
    return
  end subroutine wrp_allocate_1D_EightByteReal

  subroutine wrp_allocate_scl_FourByteInt(dmn_cnt,var_ptr,var_nm)
    ! Skeleton template based on wrp_allocate_1D_FourByteInt()
    implicit none
    character(len=*),parameter::sbr_nm='wrp_allocate_scl_FourByteInt'
    integer,intent(in)::dmn_cnt
    character(len=*),optional,intent(in)::var_nm
    integer(selected_int_kind(6)),pointer::var_ptr(:)
    integer::rcd
    if (associated(var_ptr)) then
       if(present(var_nm)) then 
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable ',var_nm,' already associated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable already associated'
       endif
       stop
    endif
    allocate(var_ptr(dmn_cnt),stat=rcd)
    if (rcd /= 0) then  
       if(present(var_nm)) then 
          write (6,'(3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating ',var_nm
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating memory'
       endif
       write (6,'(2a,i3)') sbr_nm(1:ftn_strlen(sbr_nm)),': rcd = ',rcd
       stop  
    endif
    return
  end subroutine wrp_allocate_scl_FourByteInt

  subroutine wrp_allocate_scl_EightByteInt(dmn_cnt,var_ptr,var_nm)
    ! Skeleton template based on wrp_allocate_1D_FourByteInt()
    implicit none
    character(len=*),parameter::sbr_nm='wrp_allocate_scl_EightByteInt'
    integer,intent(in)::dmn_cnt
    character(len=*),optional,intent(in)::var_nm
    integer(selected_int_kind(12)),pointer::var_ptr(:)
    integer::rcd
    if (associated(var_ptr)) then
       if(present(var_nm)) then 
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable ',var_nm,' already associated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable already associated'
       endif
       stop
    endif
    allocate(var_ptr(dmn_cnt),stat=rcd)
    if (rcd /= 0) then  
       if(present(var_nm)) then 
          write (6,'(3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating ',var_nm
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating memory'
       endif
       write (6,'(2a,i3)') sbr_nm(1:ftn_strlen(sbr_nm)),': rcd = ',rcd
       stop  
    endif
    return
  end subroutine wrp_allocate_scl_EightByteInt

#ifndef PRC_DBL 
  subroutine wrp_allocate_scl_FourByteReal(dmn_cnt,var_ptr,var_nm)
    ! Skeleton template based on wrp_allocate_1D_FourByteInt()
    implicit none
    character(len=*),parameter::sbr_nm='wrp_allocate_scl_FourByteReal'
    integer,intent(in)::dmn_cnt
    character(len=*),optional,intent(in)::var_nm
    real(selected_real_kind(p=6)),pointer::var_ptr(:)
    integer::rcd
    if (associated(var_ptr)) then
       if(present(var_nm)) then 
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable ',var_nm,' already associated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable already associated'
       endif
       stop
    endif
    allocate(var_ptr(dmn_cnt),stat=rcd)
    if (rcd /= 0) then  
       if(present(var_nm)) then 
          write (6,'(3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating ',var_nm
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating memory'
       endif
       write (6,'(2a,i3)') sbr_nm(1:ftn_strlen(sbr_nm)),': rcd = ',rcd
       stop  
    endif
    return
  end subroutine wrp_allocate_scl_FourByteReal
#endif /* PRC_DBL */

  subroutine wrp_allocate_scl_EightByteReal(dmn_cnt,var_ptr,var_nm)
    ! Skeleton template based on wrp_allocate_1D_EightByteInt()
    implicit none
    character(len=*),parameter::sbr_nm='wrp_allocate_scl_EightByteReal'
    integer,intent(in)::dmn_cnt
    character(len=*),optional,intent(in)::var_nm
    real(selected_real_kind(p=12)),pointer::var_ptr(:)
    integer::rcd
    if (associated(var_ptr)) then
       if(present(var_nm)) then 
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable ',var_nm,' already associated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR variable already associated'
       endif
       stop
    endif
    allocate(var_ptr(dmn_cnt),stat=rcd)
    if (rcd /= 0) then  
       if(present(var_nm)) then 
          write (6,'(3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating ',var_nm
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR allocating memory'
       endif
       write (6,'(2a,i3)') sbr_nm(1:ftn_strlen(sbr_nm)),': rcd = ',rcd
       stop  
    endif
    return
  end subroutine wrp_allocate_scl_EightByteReal

#ifndef PRC_DBL 
  subroutine wrp_deallocate_FourByteReal(var_ptr,rcd,var_nm_opt)
    ! Purpose: Master wrapper for deallocate(), all other instances same
    ! as this one but have comments stripped out
    ! If no error is indicated, return silently
    ! If error is indicated, print corresponding error message, return code, optional variable name, then exit
    implicit none
    character(len=*),parameter::sbr_nm='wrp_deallocate_FourByteReal'
    ! Input
    integer,intent(out)::rcd ! [rcd] Return code
    character(len=*),optional,intent(in)::var_nm_opt ! [rcd] Optional variable name
    real(selected_real_kind(p=6)),dimension(:),pointer::var_ptr ! [frc] Array to deallocate
    ! Local
    ! Main Code
    ! Initialize return code
    rcd=0 ! [rcd] Return code
    ! If no error or error is non-fatal, return value of rcd
    if (associated(var_ptr)) then
       deallocate(var_ptr,stat=rcd)
       nullify(var_ptr)
    else
       if(present(var_nm_opt)) then
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)), &
               ': WARNING variable ',var_nm_opt, &
               ' is not associated and will not be deallocated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)), &
               ': WARNING variable is not associated and will not be deallocated'
       endif
    endif ! endif associated
    ! Return if no error or if error matches exception rcd
    if (rcd == 0) return
    ! Houston, we have a problem...
    if(present(var_nm_opt)) then
       write (6,'(2a,i6,3a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR deallocate failed rcd = ',rcd, &
            ' for variable ',var_nm_opt,' in '//__FILE__

    else
       write (6,'(2a,i6,a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR deallocate failed rcd = ',rcd, &
            ' in __FILE__'//__FILE__
    endif
    ! Die, monster, die!
    stop
  end subroutine wrp_deallocate_FourByteReal ! end wrp_deallocate()
#endif /* PRC_DBL */
  
  subroutine wrp_deallocate_EightByteReal(var_ptr,rcd,var_nm_opt)
    implicit none
    character(len=*),parameter::sbr_nm='wrp_deallocate_EightByteReal'
    integer,intent(out)::rcd
    character(len=*),optional,intent(in)::var_nm_opt
    real(selected_real_kind(p=12)),dimension(:),pointer::var_ptr
    rcd=0
    if (associated(var_ptr)) then
       deallocate(var_ptr,stat=rcd)
       nullify(var_ptr)
    else
       if(present(var_nm_opt)) then
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)), &
               ': WARNING variable ',var_nm_opt, &
               ' is not associated and will not be deallocated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)), &
               ': WARNING variable is not associated and will not be deallocated'
       endif
    endif
    if (rcd == 0) return
    if(present(var_nm_opt)) then
       write (6,'(2a,i6,2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR deallocate failed rcd = ',rcd,' for variable ',var_nm_opt
    else
       write (6,'(2a,i6)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR deallocate failed rcd = ',rcd
    endif
    stop
  end subroutine wrp_deallocate_EightByteReal
  
  subroutine wrp_deallocate_FourByteInt(var_ptr,rcd,var_nm_opt)
    implicit none
    character(len=*),parameter::sbr_nm='wrp_deallocate_FourByteInt'
    integer,intent(out)::rcd
    character(len=*),optional,intent(in)::var_nm_opt
    integer(selected_int_kind(6)),dimension(:),pointer::var_ptr
    rcd=0
    if (associated(var_ptr)) then
       deallocate(var_ptr,stat=rcd)
       nullify(var_ptr)
    else
       if(present(var_nm_opt)) then
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)), &
               ': WARNING variable ',var_nm_opt, &
               ' is not associated and will not be deallocated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)), &
               ': WARNING variable is not associated and will not be deallocated'
       endif
    endif
    if (rcd == 0) return
    if(present(var_nm_opt)) then
       write (6,'(2a,i6,2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR deallocate failed rcd = ',rcd,' for variable ',var_nm_opt
    else
       write (6,'(2a,i6)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR deallocate failed rcd = ',rcd
    endif
    stop
  end subroutine wrp_deallocate_FourByteInt
  
  subroutine wrp_deallocate_EightByteInt(var_ptr,rcd,var_nm_opt)
    implicit none
    character(len=*),parameter::sbr_nm='wrp_deallocate_EightByteInt'
    integer,intent(out)::rcd
    character(len=*),optional,intent(in)::var_nm_opt
    integer(selected_int_kind(12)),dimension(:),pointer::var_ptr
    rcd=0
    if (associated(var_ptr)) then
       deallocate(var_ptr,stat=rcd)
       nullify(var_ptr)
    else
       if(present(var_nm_opt)) then
          write (6,'(4a)') sbr_nm(1:ftn_strlen(sbr_nm)), &
               ': WARNING variable ',var_nm_opt, &
               ' is not associated and will not be deallocated'
       else
          write (6,'(2a)') sbr_nm(1:ftn_strlen(sbr_nm)), &
               ': WARNING variable is not associated and will not be deallocated'
       endif
    endif
    if (rcd == 0) return
    if(present(var_nm_opt)) then
       write (6,'(2a,i6,2a)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR deallocate failed rcd = ',rcd,' for variable ',var_nm_opt
    else
       write (6,'(2a,i6)') sbr_nm(1:ftn_strlen(sbr_nm)),': ERROR deallocate failed rcd = ',rcd
    endif
    stop
  end subroutine wrp_deallocate_EightByteInt
  
end module mmr_mdl ! [mdl] Memory management
