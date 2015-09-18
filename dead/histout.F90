! $Id$ -*-f90-*-

! Purpose: Supply history/archiving routines

! Usage: 
! use histout,only:outfld ! [mdl] History/archiving

#include <params.h> /* Preprocessor tokens */

module histout ! [mdl] History/archiving
  use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
  implicit none
  private ! [stt] Symbols are private unless individually qualified as public
  public::outfld ! [sbr] Dummy routine for history tape write
  
contains
  
#ifdef CCM
  subroutine outfld(name,data,pcols,lchnk)
    ! Purpose: Dummy routine for history tape write (CCM)
    use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
    character(len=*),intent(in)::name ! [sng] Variable name
    integer,intent(in)::pcols ! [nbr] Number of columns
    integer,intent(in)::lchnk ! [idx] Chunk index
    real(r8),intent(in)::data(pcols,*) ! [frc] Data
    real(r8) CEWI_flt ! CEWI
    CEWI_flt=real(lchnk)+sum(data(:,1))+real(ichar(name(1:1))) ! CEWI
    return 
  end subroutine outfld            ! end outfld()
#else /* not CCM */
  subroutine outfld(name,data,plond,lat_idx,obuf)
    ! Purpose: Dummy routine for history tape write (MATCH)
    use shr_kind_mod,only:r8=>shr_kind_r8 ! [mdl] Precision r8, i8, ...
    character(len=*),intent(in)::name ! [sng] Variable name
    integer,intent(in)::plond ! [nbr] Number of longitudes
    integer,intent(in)::lat_idx ! [idx] Latitude index
    real(r8),intent(in)::data(plond,*) ! [frc] Data
    real(r8),intent(inout):obuf(*) ! [bfr] Output buffer
    return 
  end subroutine outfld            ! end outfld()
#endif /* not CCM */
  
end module histout
