c     $Id$
      
c     Purpose: Calculate solar position given year, day, time, latitude, and longitude. 
c     Outputs solar azimuth, elevation, hour angle, declination, and air mass.
c     This drives the slr_crd_jjm routines from terminal input.
      
c     NB: Original file ftp://hazy.asrc.albany.edu/pub/asunpos.f supplied 970414 by JJM
      
c     Reference: Joseph J. Michalsky, 1988, "The Astronomical Almanac's Algorithm for
c     Approximate Solar Position (1950--2050)", Solar Energy, v. 40, n. 3, pp. 227--235.

c     For future reference, output of the original program (before any modifications) for
c     Jan. 1, 2000, noontime on the equatator on the Greenwich meridian:
c     lat=0, lon=0, year=2000, time zone=0, day of year=1, local standard time=12.00.00.
c     solar azimuth=178.0599, elevation=66.9594, hour angle=-0.8252, declination=-23.0334,
c     air mass=1.09, solar distance=0.9833
      
      program slr_crd_jjm
c     Purpose:
c     Terminal driver for slr_crd_Mic88 routines.
      implicit none
c     Parameters
#include <parameter.com>
c     Common
c     Input
      character*20 foo_sng
      real gmt_hr_dcm 
      real gmt_hr 
      real gmt_mnt 
      real gmt_sec 
      real gmt_ydy 
      real gmt_yr 
      real lat_dgr 
      real lon_dgr
c-Input/Output Arguments
c-Output Arguments
      real air_mass
      real eqn_time_sec         ! Equation of time (second)
      real slr_azi_dgr          ! Solar azimuth (degree)
      real slr_dcl_dgr          ! Solar declination (degree)
      real slr_dmt_dgr          ! Diameter of solar disc (degree)
      real slr_dst_au           ! Solar distance (AU)
      real slr_elv_dgr          ! Solar elevation (degree)
      real slr_hr_ngl_dgr       ! Solar hour angle (degree)
      real slr_rgt_asc_dgr      ! Solar right ascension (degree)
      real slr_rfr_ngl_dgr      ! Solar refraction angle (degree)
c-Local workspace-
      double precision pi       ! [frc] 3
      real slr_cst
      real slr_flx_TOA 
      real slr_flx_nrm_TOA 
      real slr_zen_ngl 
      real slr_zen_ngl_cos 
      real slr_zen_ngl_dgr 
      real xnt_fac
c     External
      real airmass
      external airmass
c     Main code

c     Initialize defaults
      pi=4.*atan(1.)

c     ARM SGP CART site
      lat_dgr=36.6048
      lon_dgr=-97.4841

c     Local solar noon, Oct 11, 1995
      gmt_yr=1995
      gmt_ydy=284
      gmt_hr=18
      gmt_mnt=16
      gmt_sec=17

      write (6,'("Ready for coordinate input...")')

      write (6,'("Latitude (degrees North, default is ",f15.12,"):")') lat_dgr
      read (5,'(a)') foo_sng
      if (foo_sng.ne.'') read (foo_sng,'(f15.12)') lat_dgr
      write (6,'(f15.12)') lat_dgr
      foo_sng=''

      write (6,'("Longitude (degrees East, default is ",f15.11,"):")') lon_dgr
      read (5,'(a)') foo_sng
      if (foo_sng.ne.'') read (foo_sng,'(f15.11)') lon_dgr
      write (6,'(f15.11)') lon_dgr
      foo_sng=''

      write (6,'("Year (Christian, default is ",f5.0,"):")') gmt_yr
      read (5,'(a)') foo_sng
      if (foo_sng.ne.'') read (foo_sng,'(f5.0)') gmt_yr
      write (6,'(f5.0)') gmt_yr
      foo_sng=''

 60   continue
      write (6,'("Julian day of year, e.g., Feb 1=32 (default is ",f15.11,")")') gmt_ydy
      read (5,'(a)') foo_sng
      if (foo_sng.ne.'') read (foo_sng,'(f15.11)') gmt_ydy
      write (6,'(f15.11)') gmt_ydy
      foo_sng=''

 50   continue
      write (6,'("GMT hour (default is ",f15.11,"):")') gmt_hr
      read (5,'(a)') foo_sng
      if (foo_sng.ne.'') read (foo_sng,'(f15.11)') gmt_hr
      write (6,'(f15.11)') gmt_hr
      foo_sng=''

      write (6,'("GMT minute (default is ",f15.11,"):")') gmt_mnt
      read (5,'(a)') foo_sng
      if (foo_sng.ne.'') read (foo_sng,'(f15.11)') gmt_mnt
      write (6,'(f15.11)') gmt_mnt
      foo_sng=''

      write (6,'("GMT second (default is ",f15.11,"):")') gmt_sec
      read (5,'(a)') foo_sng
      if (foo_sng.ne.'') read (foo_sng,'(f15.11)') gmt_sec
      write (6,'(f15.11)') gmt_sec
      foo_sng=''

      gmt_hr_dcm=(gmt_hr*3600.0+gmt_mnt*60.0+gmt_sec)/3600.0
      call slr_crd_Mic88(
c     Input
     $     real(gmt_hr_dcm),
     $     real(gmt_ydy),
     $     real(gmt_yr),
     $     real(lat_dgr),
     $     real(lon_dgr),
c     Output
     $     eqn_time_sec,        ! Equation of time (second)
     $     slr_azi_dgr,         ! Solar azimuth (degree)
     $     slr_dcl_dgr,         ! Solar declination (degree)
     $     slr_dmt_dgr,         ! Diameter of solar disc (degree)
     $     slr_dst_au,          ! Solar distance (AU)
     $     slr_elv_dgr,         ! Solar elevation (degree)
     $     slr_hr_ngl_dgr,      ! Solar hour angle (degree)
     $     slr_rgt_asc_dgr,     ! Solar right ascension (degree)
     $     slr_rfr_ngl_dgr      ! Solar refraction angle (degree)
     $     )

      air_mass=airmass(slr_elv_dgr)

      slr_zen_ngl_dgr=90.-slr_elv_dgr
      slr_zen_ngl=pi*slr_zen_ngl_dgr/180.
      slr_zen_ngl_cos=cos(slr_zen_ngl)
      xnt_fac=1./(slr_dst_au*slr_dst_au)
      slr_cst=slr_cst_CCM
      slr_flx_nrm_TOA=slr_cst*xnt_fac
      slr_flx_TOA=slr_cst*xnt_fac*slr_zen_ngl_cos

      write (6,'("Solar zenith angle (degrees): ",f8.4)') slr_zen_ngl_dgr
      write (6,'("Solar zenith angle cosine: ",f8.7)') slr_zen_ngl_cos
      write (6,'("Solar constant (W/m2): ",f8.3)') slr_cst
      write (6,'("Solar distance (AU): ",f8.6)') slr_dst_au
      write (6,'("Eccentricity factor: ",f8.6)') xnt_fac
      write (6,'("Solar flux adjusted for eccentricity (W/m2): ",f8.3)') slr_flx_nrm_TOA
      write (6,'("Solar flux TOA (W/m2): ",f8.3)') slr_flx_TOA

      write (6,'("Solar azimuth (degrees): ",f8.4)') slr_azi_dgr
      write (6,'("Solar elevation (degrees): ",f8.4)') slr_elv_dgr
      write (6,'("Solar right ascension (degrees): ",f8.4)') slr_rgt_asc_dgr
      write (6,'("Solar declination (degrees): ",f8.4)') slr_dcl_dgr
      write (6,'("Solar hour angle (degrees): ",f8.4)') slr_hr_ngl_dgr
      write (6,'("Equation of time (s): ",f8.3)') eqn_time_sec
      write (6,'("Solar refraction angle (degrees): ",f8.4)') slr_rfr_ngl_dgr
      write (6,'("Diameter of solar disc (degrees): ",f8.4)') slr_dmt_dgr
      write (6,'("Air mass : ",f5.3)') air_mass

c      write (6,'("Type 9 for another time on same day, same location:")')
c      read (5,'(i1)') int_foo
c      if (int_foo.eq.9) go to 50
c      write (6,'("Type 99 for a different day and time, same location:")')
c      read (5,'(i2)') int_foo
c      if (int_foo.eq.99) go to 60

      end                       ! end slr_crd_jjm()
      
