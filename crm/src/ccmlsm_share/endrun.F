#include <misc.h>
      subroutine endrun
C-----------------------------------------------------------------------
C
C Abort the model for abnormal termination
C
C---------------------------Code history--------------------------------
C
C Original version:  L. Bath, Apr 1992
C Standardized:      L. Bath, Jun 1992
C                    L. Buja, Feb 1996
C
C-----------------------------------------------------------------------
c
c $Id$
c $Author$
c
C-----------------------------------------------------------------------
      implicit none
C-----------------------------------------------------------------------
C
#if ( defined SPMD )
      integer info

      call MPI_FINALIZE(info)
#endif
      call abort
C
      end
 
