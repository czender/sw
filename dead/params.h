! $Id$ -*-f90-*- 

! Purpose: params.h sets all tokens used in dust parameterization

! Usage: 
! #include <params.h> /* Preprocessor tokens */ 

! DST_NBR is number of dust tracers
! DST_NBR = 4 selects default CCM/MATCH/CAM grid (ZBN03)
! DST_NBR = 8 selects default WRF grid
! DST_IDX_SRT is array index of first dust component in constituent array
! Provide defaults for DST_IDX_SRT and DST_NBR even when DST is not defined
! This permits files in dst to compile for non-dust control runs
! CCM and MATCH differ significantly in their enumeration of advected species
! CCM guarantees that water vapor is species 1, thus DST_IDX_SRT may be 2
! MATCH (with prognostic cloudwater) sets condensed water to species PCNST 
! MATCH carries water vapor in a separate array (shadv), thus DST_IDX_SRT is 1

#ifndef PARAMS_H
#define PARAMS_H

! Define Box Model parameters if appropriate, else default to global models
#ifdef BXM

! Begin Box Model section
#define PLAT 1 
#define PLON 1 
#define PLEV 1
#ifndef DST_NBR
# define DST_NBR 4
!# define DST_NBR 8
#endif /* DST_NBR */
#define DST_IDX_SRT 2
#define PCNST 1+DST_NBR

! End Box Model section
#else /* not BXM */ 

! Define CCM parameters if appropriate, else default to MATCH/MOZART
#ifdef CCM 

! Begin CCM section
#ifdef DST
#ifndef DST_NBR
# define DST_NBR 4
#endif /* DST_NBR */
#define DST_IDX_SRT 2
#define PCNST 1+DST_NBR
#else /* not DST */
#define DST_IDX_SRT 1
#define DST_NBR 1 
#define PCNST 1
#endif /* not DST */

#define PNATS 0
#define PLEV 18
#define PLEVR 18
#define POZLEV 23

! End CCM section
#else /* not CCM */ 

! Begin MATCH/MOZART section
#ifdef DST
#ifndef DST_NBR
# define DST_NBR 4
#endif /* DST_NBR */
#define DST_IDX_SRT 1
#define PCNST 1+DST_NBR
#else /* not DST */
#define DST_IDX_SRT 1
#define DST_NBR 1 
#define PCNST 2
#endif /* not DST */

#define PNATS  0
#define PLEV 28

#define CALC_ETADOT
#define DI_VDIFF
#define DI_CONV_CCM
#define DI_CLOUD_PHYS

! End MATCH/MOZART section
#endif /* not CCM */ 

! Define Resolution parameters for global models
! NB: MATCH/MOZART do not need PTRM, PTRN, or PTRK
#ifdef T5
#define PLON 16
#define PLAT 8
#define PTRM 5
#define PTRN 5
#define PTRK 5
#endif /* not T5 */

#ifdef T21
#define PLON 64
#define PLAT 32
#define PTRM 21
#define PTRN 21
#define PTRK 21
#endif /* not T21 */

#ifdef T31
#define PLON 96
#define PLAT 48
#define PTRM 31
#define PTRN 31
#define PTRK 31
#endif /* not T31 */

#ifdef T42
#define PLON 128
#define PLAT 64
#define PTRM 42
#define PTRN 42
#define PTRK 42
#endif /* not T42 */

#ifdef T62
#define PLON 192
#define PLAT 94
#define PTRM 62
#define PTRN 62
#define PTRK 62
#endif /* not T62 */

! End Global Model section
#endif /* not BXM */ 

#endif /* not PARAMS_H */ 

