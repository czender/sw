! $Id$

! Purpose: Compute solar spectral radiative fluxes

! Copyright (C) 1994--2017 Charlie Zender
! This software is distributed under the terms of the GNU General Public License
! See http://www.gnu.org/copyleft/gpl.html for full license text

! This file contains four routines to return the fraction
! of the solar flux contained within a give wavelength interval.
! Source code is all from Bruce P. Briegleb. There are 4 functions,
! the first two give data from Thekeakara and Drummond (ThD71), 
! the last two give data from Labs and Neckel (LaN68). 

! slfftd(lmin,lmax): lmin,lmax are wavelengths in microns.
! slfwtd(wavhgh,wavlow): wavhgh,wavlow are higest,lowest wavenumbers
! slffln(lmin,lmax): lmin,lmax are wavelengths in microns.
! slfwln(wavhgh,wavlow): wavhgh,wavlow are higest,lowest wavenumbers

! ThD71 binned their data every 0.01 microns from 0.2 up to 1.0 micron 
! then by 0.2 microns up to 5.0 microns for a total of 120 intervals. 
! Valid wavelengths are between 0.15 and 4.95 microns inclusive.

! LaN68 binned their data every 0.01 microns from 0.2 up to 0.6 microns 
! then by 0.1 microns up to 5 microns for a total of 120 intervals. 
! Valid wavelengths are between 0.2 and 5.0 microns inclusive.

! Users of ThD71: TaL89, Ste78, Sli89
! Users of LaN68: CCM2, 
! Users of NeL84: EbC92

! Notes on PRC_DBL, auto-promotion, overloading:
! Assume PRC_DBL is defined when compiler is passed autopromotion flags ("-r8")
! This means all floats are promoted to double precision
! As a result, ftn_arg_get_flt() is identical to ftn_arg_get_dbl() in precision
! Thus compiler cannot disambiguate overloading these two functions
! Result is OK on all known compilers except Lahey lf95 which complains
! Solution is to define ftn_arg_get_flt() only when auto-promotion is not enabled
! Alternate solution is never use auto-promotion with lf95

! Usage: 
! use flx_slr_mdl ! [mdl] Solar spectral fluxes

! fxm: ifc produces Comment 1 : In procedure SLFWLN_DBL label 201 is not referenced
! fxm: ifc produces Comment 1 : In procedure SLFWLN_DBL label 101 is not referenced

module flx_slr_mdl ! [mdl] Solar spectral fluxes
  implicit none
  public::slfftd,slffln,slfwtd,slfwln
#ifndef PRC_DBL 
  private::slfftd_flt,slffln_flt,slfwtd_flt,slfwln_flt
#endif /* PRC_DBL */
  private::slfftd_dbl,slffln_dbl,slfwtd_dbl,slfwln_dbl
  private::r4 ! r4: 4B (C float) default, 8B (C double) possible
  private::r8 ! r8: 8B (C double) default, 4B (C float) possible

  integer,parameter::r4=selected_real_kind(p=6) ! r4: 4B (C float) default, 8B (C double) possible
  integer,parameter::r8=selected_real_kind(p=12) ! r8: 8B (C double) default, 4B (C float) possible
  
  ! Overloaded spectral solar flux retrieval functions
  interface slfftd
#ifndef PRC_DBL 
     module procedure slfftd_flt
#endif /* PRC_DBL */
     module procedure slfftd_dbl
  end interface ! slfftd
  interface slffln
#ifndef PRC_DBL 
     module procedure slffln_flt
#endif /* PRC_DBL */
     module procedure slffln_dbl
  end interface ! slffln
  interface slfwtd
#ifndef PRC_DBL 
     module procedure slfwtd_flt
#endif /* PRC_DBL */
     module procedure slfwtd_dbl
  end interface ! slfwtd
  interface slfwln
#ifndef PRC_DBL 
     module procedure slfwln_flt
#endif /* PRC_DBL */
     module procedure slfwln_dbl
  end interface ! slfwln

contains
  
#ifndef PRC_DBL 
  real(r4) function slfftd_flt(lmin,lmax)
    real(r4),intent(in)::lmin
    real(r4),intent(in)::lmax
    slfftd_flt=slfftd_dbl(real(lmin,r8),real(lmax,r8))
  end function slfftd_flt
  real(r4) function slffln_flt(lmin,lmax)
    real(r4),intent(in)::lmin
    real(r4),intent(in)::lmax
    slffln_flt=slffln_dbl(real(lmin,r8),real(lmax,r8))
  end function slffln_flt
  real(r4) function slfwtd_flt(wavhgh,wavlow)
    real(r4),intent(in)::wavhgh
    real(r4),intent(in)::wavlow
    slfwtd_flt=slfwtd_dbl(real(wavhgh,r8),real(wavlow,r8))
  end function slfwtd_flt
  real(r4) function slfwln_flt(wavhgh,wavlow)
    real(r4),intent(in)::wavhgh
    real(r4),intent(in)::wavlow
    slfwln_flt=slfwln_dbl(real(wavhgh,r8),real(wavlow,r8))
  end function slfwln_flt
#endif /* PRC_DBL */

  real(r8) function slfftd_dbl(lmin,lmax)
    ! Purpose: Fractional solar flux between wavelengths lmin,lmax
    ! Solar flux fractions according to Thekeakara and Drummond (1971)
    ! Integrating values from 0.200 microns to 4.950 microns gives solar constant = 1351.3 W m-2
    ! Values below are actually fraction of solar flux shortwave of particular wavelength
    implicit none
    integer,parameter::maxwav=67
    real(r8),intent(in)::lmin
    real(r8),intent(in)::lmax
    real(r8),dimension(maxwav),parameter::lambda=(/ &
         .15 , &
         .20 , .22 , .23 , .24 , .25 , .26  ,  .27 , &
         .28 , &
         .29 , .30 , .31 , .32 , .33 , .34  ,  .35 , &
         .36 , .37 , .38 , .39 , .40 , .41  ,  .42 , &
         .43 , .44 , .45 , .46 , .47 , .48  ,  .49 , &
         .50 , .51 , .52 , .53 , .54 , .55  ,  .56 , &
         .57 , &
         .58 , .59 , .60 , .62 , .64 , .66  ,  .68 , &
         .70 , .72 , .75 , .80 , .90 ,1.00  , 1.20 , &
         1.40 ,1.60 ,1.80 ,2.00 ,2.20 ,2.40  , 2.60 , &
         2.80 ,3.00 ,3.20 ,3.40 ,3.60 ,3.80  , 4.00 , &
         5.00 /)
    real(r8),dimension(maxwav),parameter::dlambd=(/ &
         .0000, &
         .0001,.0005,.0010,.0014,.0019,.0027 ,.0041 , &
         .0056, &
         .0081,.0121,.0165,.0222,.0293,.0372 ,.0452 , &
         .0532,.0615,.0700,.0782,.0873,.0992 ,.1122 , &
         .1247,.1373,.1514,.1665,.1817,.1968 ,.2115 , &
         .2260,.2401,.2538,.2674,.2808,.2938 ,.3065 , &
         .3191, &
         .3318,.3444,.3568,.3810,.4042,.4266 ,.4481 , &
         .4688,.4886,.5169,.5602,.6336,.6946 ,.7839 , &
         .8434,.8861,.9159,.9349,.9483,.9589 ,.9667 , &
         .9731,.9783,.9822,.9850,.9872,.9891 ,.9906 , &
         .9951 /)
    integer::n
    ! Locals with simple initialization and no command-line override
    real(r8)::dmin=0.0 ! CEWI
    real(r8)::dmax=0.0 ! CEWI
    
    if(lmin >= 0.150 .and. lmin < 4.950) then
       
       ! find lower limit
       n=1
100    if(lmin >= lambda(n)) then
          n=n+1
          go to 100
       endif
       dmin=dlambd(n-1)+((dlambd(n)-dlambd(n-1))/(lambda(n) &
           -lambda(n-1)))*(lmin-lambda(n-1))
    else if(lmin < 0.150) then
       dmin=0.0
    else if(lmin > 4.950) then
       slfftd_dbl=0.0
       return
    endif
    
    ! find upper limit
    if(lmax > 0.200 .and. lmax <= 4.950) then
       n=1
200    if(lmax >= lambda(n)) then
          n=n+1
          go to 200
       endif
       dmax=dlambd(n-1)+((dlambd(n)-dlambd(n-1))/(lambda(n) &
           -lambda(n-1)))*(lmax-lambda(n-1))
    else if(lmax > 4.950) then
       dmax=1.0000
       ! csz++
       ! else if(lmax < 0.150) then
    else if(lmax <= 0.2) then
       ! csz--
       slfftd_dbl=0.0
       return
    endif
    
    ! solar flux fraction
    slfftd_dbl=dmax-dmin
    return
  end function slfftd_dbl
  
  real(r8) function slfwtd_dbl(wavhgh,wavlow)
    ! Purpose: Fractional solar flux between wavelnumbers wavhgh,wavlow
    ! Purpose: Compute solar flux fraction between highest (wavhgh) and lowest (wavlow) wavenumber
    ! Solar flux fractions according to Thekeakara and Drummond (1971)
    ! Integrating values from 0.200 microns to 4.950 microns gives solar constant = 1351.3 W m-2
    ! Values below are actually fraction of solar flux shortwave of particular wavelength
    implicit none
    integer,parameter::maxwav=67
    real(r8),intent(in)::wavhgh
    real(r8),intent(in)::wavlow
    real(r8) wavnum(maxwav)
    integer::iter=0
    real(r8),dimension(maxwav),parameter::lamdat=(/ &
         .15 , &
         .20 , .22 , .23 , .24 , .25 , .26  ,  .27 , &
         .28 , &
         .29 , .30 , .31 , .32 , .33 , .34  ,  .35 , &
         .36 , .37 , .38 , .39 , .40 , .41  ,  .42 , &
         .43 , .44 , .45 , .46 , .47 , .48  ,  .49 , &
         .50 , .51 , .52 , .53 , .54 , .55  ,  .56 , &
         .57 , &
         .58 , .59 , .60 , .62 , .64 , .66  ,  .68 , &
         .70 , .72 , .75 , .80 , .90 ,1.00  , 1.20 , &
         1.40 ,1.60 ,1.80 ,2.00 ,2.20 ,2.40  , 2.60 , &
         2.80 ,3.00 ,3.20 ,3.40 ,3.60 ,3.80  , 4.00 , &
         5.00 /)
    real(r8),dimension(maxwav),parameter::solwav=(/ &
         .0000, &
         .0001,.0005,.0010,.0014,.0019,.0027 ,.0041 , &
         .0056, &
         .0081,.0121,.0165,.0222,.0293,.0372 ,.0452 , &
         .0532,.0615,.0700,.0782,.0873,.0992 ,.1122 , &
         .1247,.1373,.1514,.1665,.1817,.1968 ,.2115 , &
         .2260,.2401,.2538,.2674,.2808,.2938 ,.3065 , &
         .3191, &
         .3318,.3444,.3568,.3810,.4042,.4266 ,.4481 , &
         .4688,.4886,.5169,.5602,.6336,.6946 ,.7839 , &
         .8434,.8861,.9159,.9349,.9483,.9589 ,.9667 , &
         .9731,.9783,.9822,.9850,.9872,.9891 ,.9906 , &
         .9951 /)
    integer::n
    integer::m
    real(r8)::dmin
    real(r8)::dmax
    
    if(iter == 0) then
       do n=1,maxwav
          wavnum(n)=10000.0/lamdat(n)
       end do
    endif
    iter=iter+1
    
    if(wavlow > 50000.0) then
       slfwtd_dbl=0.0
       return
    endif
    
    if(wavlow >= 2000.0) then
       m=1
       do n=2,maxwav
          
          ! Wavenumbers ordered in array 'wavnum' from large to small:
          if(wavnum(n-1) >= wavlow .and. wavlow >= wavnum(n)) then
             m=n
             goto 101
          endif
       end do
       
101    continue
       dmax=solwav(m-1)+((solwav(m)-solwav(m-1))/(wavnum(m) &
           -wavnum(m-1)))*(wavlow-wavnum(m-1))
    else
       dmax=1.0
    endif
    
    if(wavhgh < 2000.0) then
       slfwtd_dbl=0.0
       return
    endif
    
    if(wavhgh <= 50000.0) then
       m=1
       do n=2,maxwav
          
          ! Wavenumbers ordered in array 'wavnum' from large to small:
          if(wavnum(n-1) >= wavhgh .and. wavhgh >= wavnum(n)) then
             m=n
             goto 201
          endif
       end do
       
201    continue
       dmin=solwav(m-1)+((solwav(m)-solwav(m-1))/(wavnum(m) &
           -wavnum(m-1)))*(wavhgh-wavnum(m-1))
    else
       dmin=0.0
    endif
    
    ! Solar flux fraction
    slfwtd_dbl=dmax-dmin
    return
  end function slfwtd_dbl
  
  real(r8) function slffln_dbl(lmin,lmax)
    ! Purpose: Fractional solar flux between wavelengths lmin,lmax
    ! Solar flux fractions according to Labs and Neckel (1968)
    ! Integrating values from 0.200 microns to 4.950 microns gives solar constant = 1360.39 W m-2
    ! Values below are actually fraction of solar flux shortwave of particular wavelength
    implicit none
    integer,parameter::maxwav=122
    real(r8),intent(in)::lmin
    real(r8),intent(in)::lmax
    real(r8),dimension(maxwav),parameter::lambda=(/ &
         .200,.205,.215,.225,.235,.245,.255,.265,.275, &
         .285,.295, &
         .305,.315,.325,.335,.345,.355,.365,.375,.385,.395, &
         .405,.415,.425,.435,.445,.455,.465,.475,.485,.495, &
         .505,.515,.525,.535,.545,.555,.565,.575,.585,.595, &
         .605,.615,.625,.635,.645,.655,.665,.675,.685,.695, &
         .705,.715,.725,.735,.745,.755,.765,.775,.785,.795, &
         .805,.815,.825,.835,.845,.855,.865,.875,.885,.895, &
         .905,.915,.925,.935,.945,.955,.965,.975,.985,.995, &
         1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95, &
         2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95, &
         3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95, &
         4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95, &
         5.00 /)
    real(r8),dimension(maxwav),parameter::dlambd=(/ &
         .0000, &
         .0001,.0003,.0006,.0010,.0015,.0020,.0029,.0042,.0059,.0088, &
         .0127,.0172,.0227,.0291,.0358,.0427,.0502,.0580,.0654,.0732, &
         .0835,.0959,.1085,.1209,.1343,.1489,.1638,.1786,.1930,.2073, &
         .2217,.2356,.2493,.2633,.2774,.2911,.3047,.3183,.3318,.3451, &
         .3581,.3709,.3834,.3956,.4076,.4191,.4304,.4417,.4528,.4636, &
         .4741,.4844,.4945,.5043,.5139,.5232,.5324,.5414,.5502,.5588, &
         .5672,.5755,.5835,.5913,.5989,.6062,.6134,.6204,.6273,.6340, &
         .6407,.6472,.6536,.6598,.6660,.6719,.6778,.6835,.6892,.6947, &
         .7230,.7671,.8034,.8336,.8590,.8808,.8992,.9144,.9269,.9372, &
         .9457,.9529,.9590,.9642,.9686,.9724,.9757,.9786,.9811,.9834, &
         .9853,.9870,.9885,.9899,.9911,.9922,.9931,.9940,.9948,.9955, &
         .9962,.9968,.9973,.9978,.9982,.9987,.9990,.9994,.9997,.9999, &
         1.0 /)
    integer::n
    ! Locals with simple initialization and no command-line override
    real(r8)::dmin=0.0 ! CEWI
    real(r8)::dmax=0.0 ! CEWI

    if(lmin >= 0.200 .and. lmin < 5.000) then
       
       ! find lower limit
       n=1
100    if(lmin >= lambda(n)) then
          n=n+1
          go to 100
       endif
       dmin=dlambd(n-1)+((dlambd(n)-dlambd(n-1))/(lambda(n) &
           -lambda(n-1)))*(lmin-lambda(n-1))
    else if(lmin < 0.200) then
       dmin=0.0
       ! csz++
       ! else if(lmin > 5.000) then
    else if(lmin >= 5.000) then
       ! csz--
       slffln_dbl=0.0
       return
    endif
    
    ! find upper limit
    
    ! csz++
    ! if(lmax <= 5.000) then
    if(lmax > 0.200 .and. lmax < 5.0) then
       ! csz--
       n=1
200    if(lmax >= lambda(n)) then
          n=n+1
          go to 200
       endif
       dmax=dlambd(n-1)+((dlambd(n)-dlambd(n-1))/(lambda(n) &
           -lambda(n-1)))*(lmax-lambda(n-1))
    else if(lmax >= 5.000) then
       dmax=1.0000
       ! csz++
       ! else if(lmax < 0.200) then
    else if(lmax <= 0.200) then
       ! csz--
       slffln_dbl=0.0
       return
    endif
    
    ! solar flux fraction
    slffln_dbl=dmax-dmin
    return
  end function slffln_dbl
  
  real(r8) function slfwln_dbl(wavhgh,wavlow)
    ! Purpose: Fractional solar flux between wavenumbers wavhgh,wavlow
    ! Solar flux fractions according to Labs and Neckel (1968)
    ! Integrating values from 0.200 microns to 4.950 microns gives solar constant = 1360.39 W m-2
    ! Values below are actually fraction of solar flux shortwave of particular wavelength
    implicit none
    integer,parameter::maxwav=122
    real(r8),intent(in)::wavhgh
    real(r8),intent(in)::wavlow
    real(r8) wavnum(maxwav)
    integer::iter=0
    real(r8),dimension(maxwav),parameter::lamdat=(/ &
         .200,.205,.215,.225,.235,.245,.255,.265,.275, &
         .285,.295, &
         .305,.315,.325,.335,.345,.355,.365,.375,.385,.395, &
         .405,.415,.425,.435,.445,.455,.465,.475,.485,.495, &
         .505,.515,.525,.535,.545,.555,.565,.575,.585,.595, &
         .605,.615,.625,.635,.645,.655,.665,.675,.685,.695, &
         .705,.715,.725,.735,.745,.755,.765,.775,.785,.795, &
         .805,.815,.825,.835,.845,.855,.865,.875,.885,.895, &
         .905,.915,.925,.935,.945,.955,.965,.975,.985,.995, &
         1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95, &
         2.05,2.15,2.25,2.35,2.45,2.55,2.65,2.75,2.85,2.95, &
         3.05,3.15,3.25,3.35,3.45,3.55,3.65,3.75,3.85,3.95, &
         4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85,4.95, &
         5.00 /)
    real(r8),dimension(maxwav),parameter::solwav=(/ &
         .0000, &
         .0001,.0003,.0006,.0010,.0015,.0020,.0029,.0042,.0059,.0088, &
         .0127,.0172,.0227,.0291,.0358,.0427,.0502,.0580,.0654,.0732, &
         .0835,.0959,.1085,.1209,.1343,.1489,.1638,.1786,.1930,.2073, &
         .2217,.2356,.2493,.2633,.2774,.2911,.3047,.3183,.3318,.3451, &
         .3581,.3709,.3834,.3956,.4076,.4191,.4304,.4417,.4528,.4636, &
         .4741,.4844,.4945,.5043,.5139,.5232,.5324,.5414,.5502,.5588, &
         .5672,.5755,.5835,.5913,.5989,.6062,.6134,.6204,.6273,.6340, &
         .6407,.6472,.6536,.6598,.6660,.6719,.6778,.6835,.6892,.6947, &
         .7230,.7671,.8034,.8336,.8590,.8808,.8992,.9144,.9269,.9372, &
         .9457,.9529,.9590,.9642,.9686,.9724,.9757,.9786,.9811,.9834, &
         .9853,.9870,.9885,.9899,.9911,.9922,.9931,.9940,.9948,.9955, &
         .9962,.9968,.9973,.9978,.9982,.9987,.9990,.9994,.9997,.9999, &
         1.0000 /)
    integer::n
    integer::m
    real(r8)::dmin
    real(r8)::dmax
    
    if(iter == 0) then
       do n=1,maxwav
          wavnum(n)=10000.0/lamdat(n)
       end do
    endif
    iter=iter+1
    
    if(wavlow > 50000.0) then
       slfwln_dbl=0.0
       return
    endif
    
    if(wavlow >= 2000.0) then
       m=1
       do n=2,maxwav
          
          ! Wavenumbers ordered in array 'wavnum' from large to small:
          if(wavnum(n-1) >= wavlow .and. wavlow >= wavnum(n)) then
             m=n
          endif
       end do
       
       ! 101    continue
       dmax=solwav(m-1)+((solwav(m)-solwav(m-1))/(wavnum(m) &
           -wavnum(m-1)))*(wavlow-wavnum(m-1))
    else
       dmax=1.0
    endif
    
    if(wavhgh < 2000.0) then
       slfwln_dbl=0.0
       return
    endif
    
    if(wavhgh <= 50000.0) then
       m=1
       do n=2,maxwav
          
          ! wavenumbers ordered in array 'wavnum' from large to small:
          if(wavnum(n-1) >= wavhgh .and. wavhgh >= wavnum(n)) then
             m=n
          endif
       end do
       
       ! 201    continue
       dmin=solwav(m-1)+((solwav(m)-solwav(m-1))/(wavnum(m) &
           -wavnum(m-1)))*(wavhgh-wavnum(m-1))
    else
       dmin=0.0
    endif
    
    ! solar flux fraction
    slfwln_dbl=dmax-dmin
    return
  end function slfwln_dbl
  
end module flx_slr_mdl ! [mdl] Solar spectral fluxes
