#include <misc.h>
      subroutine freemem(address)
C-----------------------------------------------------------------------
C
C Interface routine to free dynamically allocated memory.  
C Replaces automatic arrays on Cray.  Isolates machine-specific code 
C for portable model.
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Rosinski April, 1993
C                    T. Acker, March 1996
C
C-----------------------------------------------------------------------
c
c $Id: freemem.F,v 1.4 1998/03/04 23:56:06 jet Exp $
c $Author: jet $
c
C-----------------------------------------------------------------------
      implicit none
C-----------------------------Arguments---------------------------------
C
C Input
C
      pointer (address,space)   ! Pointer to dynamically allocated space
      real space
C
C--------------------------Local Variables------------------------------
C
#if ( defined CRAY ) || ( defined T3D ) || ( defined T3E )
      integer ier               ! Error return code
#endif
C
c++ csz 20100707
c Intended to fix link problem on AIX xlf95
c#if ( defined SUN ) || ( defined RS6K )
#if ( defined RS6K )
      integer free
      external free
#endif
c-- csz 20100707
C-----------------------------------------------------------------------
#if ( defined CRAY ) || ( defined T3D ) || ( defined T3E )
CMIC$ GUARD
      call hpdeallc(address,ier,0)
CMIC$ END GUARD
C
C write error message to std. output
C
      if (ier.ne.0) then
         write(6,*)'FREEMEM: Bad return code = ',ier,
     $             ' from hpdeallc: Abort'
         if (ier.eq.-3) then
            write(6,*) address,' is outside heap bounds'
         else if (ier.eq.-4) then
            write(6,*)'Block is already free'
         else if (ier.eq.-5) then
            write(6,*) address,' is not at the beginning of a block'
         else if (ier.eq.-7) then
            write(6,*)'Control word for the next block overwritten'
         end if
         call endrun
      end if
#endif
c++ csz 20100707
c Intended to fix link problem on AIX xlf95
c#if ( defined RS6K )
c      address = free(%val(address))
c#endif
c-- csz 20100707
c#if ( defined SUN )
c      address = free(address)
c#endif
#if ( defined SGI ) || ( defined HP )
      call free(address)
#endif
C
      return
C
      end


 
