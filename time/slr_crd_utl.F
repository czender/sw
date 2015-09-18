c     $Id$
      
c     Purpose: Utilities to compute solar geometry
      
c     NB: Original file ftp://hazy.asrc.albany.edu/pub/asunpos.f supplied 970414 by JJM
      
c     Reference: Joseph J. Michalsky, 1988, "The Astronomical Almanac's Algorithm for
c     Approximate Solar Position (1950--2050)", Solar Energy, v. 40, n. 3, pp. 227--235.
      
      subroutine slr_crd_Mic88_wrp(gmt_yr,gmt_ydy,gmt_hr_dcm,lat_dgr,lon_dgr,
     $     slr_azi_dgr,slr_elv_dgr,slr_hr_ngl_dgr,slr_dcl_dgr,slr_dst_au)
c     Purpose:
c     Wrapper for slr_crd_Mic88() routine.
c     Takes date and position and returns solar coordinates and TOA flux.
      implicit none
c------------------------------Parameters-------------------------------
c------------------------------Commons----------------------------------
c---------------------------Input Arguments-----------------------------
      real
     $     gmt_hr_dcm,
     $     gmt_ydy,
     $     gmt_yr,
     $     lat_dgr,
     $     lon_dgr
c---------------------------Input/Output Arguments----------------------
c---------------------------Output Arguments----------------------------
      real
     $     slr_azi_dgr,
     $     slr_dcl_dgr,
     $     slr_dst_au,          ! Solar distance in AUs
     $     slr_elv_dgr,
     $     slr_hr_ngl_dgr
c---------------------------Local workspace-----------------------------
      
c------------------------------Externals--------------------------------
c------------------------------Main code--------------------------------
      
      return
      end                       ! end slr_crd_Mic88_wrp()
      
      subroutine slr_crd_Mic88(
c     Inputs:
     $     gmt_hr_dcm,          ! GMT decimal hour of day
     $     gmt_ydy,             ! GMT integer day of year
     $     gmt_yr,              ! GMT year
     $     lat_dgr,             ! Latitude (degree)
     $     lon_dgr,             ! Longitude (degree)
c     Outputs:
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

c     Purpose:
c     Calculate local azimuth and elevation of the sun at a specific location 
c     and time using an approximation to equations used to generate tables in 
c     The Astronomical Almanac, U.S. Gov't Printing Office, Washington, D.C. (1985).
c     Refraction correction is added so sun position is apparent one. 
      
c     Limitations:
c     These formulae were constructed from formulae accurate for the epoch of the year 2000.
c     Thus, the highest quality output is for years between 1950--2050.

c     Note on input:
c     Currently, all time input is assumed to be of type real even though two of the
c     parameters, gmt_yr and gmt_ydy, are representable as natural numbers.
      
c     Input parameters:
c     gmt_yr=year, e.g., 1986                                           
c     gmt_ydy=day of year, e.g., feb 1=32                                 
c     gmt_hr_dcm=hours plus fraction in UT, e.g., 
c     8:30 am eastern daylight time is equal to 
c     8.5 + 5 (5 hours west of Greenwich) -1 (for daylight savings time correction)
c     lat_dgr=latitude in degrees (north is positive)                     
c     lon_dgr=longitude in degrees (east is positive)
      
      implicit none
c------------------------------Parameters-------------------------------
c------------------------------Commons----------------------------------
c---------------------------Input Arguments-----------------------------
      real
     $     gmt_hr_dcm,          ! GMT decimal hour of day
     $     gmt_ydy,             ! GMT integer day of year
     $     gmt_yr,              ! GMT year
     $     lat_dgr,             ! Latitude (degree)
     $     lon_dgr              ! Longitude (degree)
c---------------------------Input/Output Arguments----------------------
c---------------------------Output Arguments----------------------------
      real
     $     eqn_time_sec,        ! Equation of time (second)
     $     slr_azi_dgr,         ! Solar azimuth (degree)
     $     slr_dcl_dgr,         ! Solar declination (degree)
     $     slr_dmt_dgr,         ! Diameter of solar disc (degree)
     $     slr_dst_au,          ! Solar distance (AU)
     $     slr_elv_dgr,         ! Solar elevation (degree)
     $     slr_hr_ngl_dgr,      ! Solar hour angle (degree)
     $     slr_rgt_asc_dgr,     ! Solar right ascension (degree)
     $     slr_rfr_ngl_dgr      ! Solar refraction angle (degree)
c---------------------------Local workspace-----------------------------
      real
     $     day_leap_dlt_1949,   ! Leap days since 1949/01/01
     $     denominator,
     $     ecl_lon,             ! Ecliptic longitude (radian)
     $     ecl_lon_dgr,         ! Ecliptic longitude (degree)
     $     eqn_time_dgr,        ! Equation of time (degree)
     $     gmst_hr,             ! Greenwich mean sidereal time (hour)
     $     gmt_yr_dlt_1949,     ! GMT in years since 1949/01/01
     $     jln_day_dlt_1949_dcm, ! Julian days since 1949/01/01
     $     jln_day_dlt_2000_dcm, ! Julian days since 2000/01/01
     $     lat,                 ! Latitude (radian)
     $     lmst_hr,             ! Local mean sidereal time (hour)
     $     lmst_rdn,            ! Local mean sidereal time (radian)
     $     mean_anm,            ! Mean anomaly (radian)
     $     mean_anm_dgr,        ! Mean anomaly (degree)
     $     mean_lon_dgr,        ! Mean longitude (degree)
     $     numerator,
     $     obl_ecl,             ! Obliquity of the ecliptic (radian)
     $     obl_ecl_dgr,         ! Obliquity of the ecliptic (degree)
     $     slr_azi,             ! Solar azimuth (radian)
     $     slr_dcl,             ! Solar declination (radian)
     $     slr_elv,             ! Solar elevation (radian)
     $     slr_hr_ngl,           ! Solar hour angle (radian)
     $     slr_rgt_asc,         ! Solar right ascension (radian)
     $     slr_zen_ngl_dgr      ! Solar zenith angle (degree)
      
c     Define some constants, including one to change between degrees and radians
      double precision
     $     pi,
     $     dgr2rdn,
     $     rdn2dgr,
     $     twopi
c------------------------------Externals--------------------------------
c------------------------------Main code--------------------------------
      
c     Initialize some constants
      pi=4.*atan(1.)
      twopi=2.*pi
      dgr2rdn=pi/180.
      rdn2dgr=180./pi

c     Get current julian date (actually add 2,400,000 for jln_day)
      gmt_yr_dlt_1949=gmt_yr-1949.
      day_leap_dlt_1949=aint(gmt_yr_dlt_1949/4.)
c     JD 2,432,916.5 is midnight, 0 Jan 1949 UT
      jln_day_dlt_1949_dcm=32916.5+gmt_yr_dlt_1949*365.+day_leap_dlt_1949+gmt_ydy+gmt_hr_dcm/24.
c     32916.5 is midnight 0 Jan 1949 minus 2.4e6; day_leap_dlt_1949=leap days since 1949
c     The last year of a century is not a leap year unless year is evenly divisible by 400
      if (amod(gmt_yr,100.).eq.0.0.and.amod(gmt_yr,400.).ne.0.0) jln_day_dlt_1949_dcm=jln_day_dlt_1949_dcm-1.
      
c     Calculate ecliptic coordinates
      jln_day_dlt_2000_dcm=jln_day_dlt_1949_dcm-51545.0
c     JD 51545.0 + 2.4e6 = noon 1 Jan 2000
      
c     Force mean longitude between 0 and 360 degrees
      mean_lon_dgr=280.460+.9856474*jln_day_dlt_2000_dcm
      mean_lon_dgr=mod(mean_lon_dgr,360.)
      if (mean_lon_dgr.lt.0.) mean_lon_dgr=mean_lon_dgr+360.
      
c     Mean anomaly in radians between 0 and 2*pi
      mean_anm_dgr=357.528+.9856003*jln_day_dlt_2000_dcm
      mean_anm_dgr=mod(mean_anm_dgr,360.)
      if (mean_anm_dgr.lt.0.) mean_anm_dgr=mean_anm_dgr+360.
      mean_anm=mean_anm_dgr*dgr2rdn
      
c     Compute ecliptic longitude and obliquity of ecliptic in radians
      ecl_lon_dgr=mean_lon_dgr+1.915*sin(mean_anm)+.020*sin(2.*mean_anm)
      ecl_lon_dgr=mod(ecl_lon_dgr,360.)
      if (ecl_lon_dgr.lt.0.) ecl_lon_dgr=ecl_lon_dgr+360.
      obl_ecl_dgr=23.439-.0000004*jln_day_dlt_2000_dcm
      ecl_lon=ecl_lon_dgr*dgr2rdn
      obl_ecl=obl_ecl_dgr*dgr2rdn
      
c     Calculate right ascension and declination
      numerator=cos(obl_ecl)*sin(ecl_lon)
      denominator=cos(ecl_lon)
      slr_rgt_asc=atan(numerator/denominator)
c     Force slr_rgt_asc between 0 and 2*pi
      if (denominator.lt.0) then
         slr_rgt_asc=slr_rgt_asc+pi
      elseif (numerator.lt.0) then
         slr_rgt_asc=slr_rgt_asc+twopi
      endif
      slr_rgt_asc_dgr=slr_rgt_asc*rdn2dgr
      
c     slr_dcl in radians
      slr_dcl=asin(sin(obl_ecl)*sin(ecl_lon))
      
c     Calculate Greenwich mean sidereal time in hours
      gmst_hr=6.697375+.0657098242*jln_day_dlt_2000_dcm+gmt_hr_dcm 
c     Hour not changed to sidereal time since 'jln_day_dlt_2000_dcm' includes the fractional day 
      gmst_hr=mod(gmst_hr,24.)
      if (gmst_hr.lt.0.) gmst_hr=gmst_hr+24.
      
c     Calculate local mean sidereal time in radians 
      lmst_hr=gmst_hr+lon_dgr/15.
      lmst_hr=mod(lmst_hr,24.)
      if (lmst_hr.lt.0.) lmst_hr=lmst_hr+24.
      lmst_rdn=lmst_hr*15.*dgr2rdn
      
c     Calculate hour angle in radians between -pi and pi
      slr_hr_ngl=lmst_rdn-slr_rgt_asc
      if (slr_hr_ngl.lt.-pi) slr_hr_ngl=slr_hr_ngl+twopi
      if (slr_hr_ngl.gt.pi) slr_hr_ngl=slr_hr_ngl-twopi
      
c     Change latitude to radians
      lat=lat_dgr*dgr2rdn
      
c     Calculate azimuth and elevation
      slr_elv=asin(sin(slr_dcl)*sin(lat)+cos(slr_dcl)*cos(lat)*cos(slr_hr_ngl))
      slr_azi=asin(-cos(slr_dcl)*sin(slr_hr_ngl)/cos(slr_elv))
      
c     This puts azimuth between 0 and 2*pi radians
c     Following 5 lines are J. W. Spencer's correction:
      if (sin(slr_dcl)-sin(slr_elv)*sin(lat).ge.0.) then
         if (sin(slr_azi).lt.0.) slr_azi=slr_azi+twopi
      else
         slr_azi=pi-slr_azi
      endif
      
c     Following 4 lines are JJM's original method which can fail in Southern Hemisphere.
c     c   if slr_azi_dgr=90 degrees, elcritical=asin(sin(slr_dcl)/sin(lat))
c     c    elc=asin(sin(slr_dcl)/sin(lat))
c     c    if (slr_elv.ge.elc) slr_azi=pi-slr_azi
c     c    if (slr_elv.le.elc.and.slr_hr_ngl.gt.0.) slr_azi=twopi+slr_azi
      
c     Calculate refraction correction for US standard atmosphere
c     Need to have slr_elv in degrees before calculating correction
      slr_elv_dgr=slr_elv*rdn2dgr
c     
      if (slr_elv_dgr.ge.19.225) then 
         slr_rfr_ngl_dgr=.00452*3.51823/tan(slr_elv)
      else if (slr_elv_dgr.gt.-.766.and.slr_elv_dgr.lt.19.225) then
         slr_rfr_ngl_dgr=3.51823*(.1594+.0196*slr_elv_dgr+.00002*slr_elv_dgr**2)/
     $        (1.+.505*slr_elv_dgr+.0845*slr_elv_dgr**2)
      else if (slr_elv_dgr.le.-.766) then
         slr_rfr_ngl_dgr=0.0
      end if
      
c     Adjust elevation in degrees for refraction
      slr_elv_dgr=slr_elv_dgr+slr_rfr_ngl_dgr ! Note that 3.51823=1013.25 mb/288 C
      slr_zen_ngl_dgr=90.-slr_elv_dgr
      
c     Calculate distance to sun in A.U. & diameter in degrees
      slr_dst_au=1.00014-.01671*cos(mean_anm)-.00014*cos(2.*mean_anm)
      slr_dmt_dgr=.5332/slr_dst_au
      
c     Convert slr_azi and lat to degrees before returning
      slr_azi_dgr=slr_azi*rdn2dgr
      slr_hr_ngl_dgr=slr_hr_ngl*rdn2dgr
      slr_dcl_dgr=slr_dcl*rdn2dgr

c     Compute equation of time diagnostic
      eqn_time_dgr=mean_lon_dgr-slr_rgt_asc*rdn2dgr
      eqn_time_sec=86400.*eqn_time_dgr/360.

      return
      end                       ! end slr_crd_Mic88()
      
      real function airmass(slr_elv_dgr)
c     Purpose:
c     Calculates air mass using Kasten's approximation to Bemporad's table
      
      implicit none
c------------------------------Parameters-------------------------------
c------------------------------Commons----------------------------------
c---------------------------Input Arguments-----------------------------
      real
     $     slr_elv_dgr          ! elevation angle (degrees)
c---------------------------Input/Output Arguments----------------------
c---------------------------Output Arguments----------------------------
c---------------------------Local workspace-----------------------------
      real
     $     slr_zen_ngl          ! zenith angle
c------------------------------Externals--------------------------------
c------------------------------Main code--------------------------------
      slr_zen_ngl=(90.0-slr_elv_dgr)*3.141592654/180.0
      airmass=1.0/(cos(slr_zen_ngl)+0.50572*(6.07995+slr_elv_dgr)**(-1.6364))
      return
      end                       ! end airmass()
