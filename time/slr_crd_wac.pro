;csz++

; solar coordinate routines from bill collins (wac) 961020

;Charlie:

;Call this first.  If you like, you can also return the "day" (Julian
;day of year) and "hour" (GMT hour) variables for your calculation of
;the solar zenith angle in CCM2. The position of the CART site
;is [-97.4841 E, 36.6048 N].

;Bill

;---------------------------------------------------------------------------
; compute_solar -- compute solar ephemeris information given 
;                  time (seconds since 1970), latitude, and longitude

; arguments:
;		time			Time in sec. since 1/1/1970 (GMT)
;		lat			local latitude in degrees 
;						(north is positive)
;		lon			local longitude (east of Greenwich 
;						is positive. 
;                         			i.e. Honolulu is 15.3, -157.8)
;		coszen			returned by ephemeris
;		solar_flux		returned by ephemeris
;		sol_azim		returned by ephemeris

; author: Bill Collins
; date: 10/24/95

pro arm2ltst,base_time,time_offset, $
	gmt_doy=gmt_doy, $
	ltst_doy=ltst_doy, $
	cos_slr_zen_ang=cos_slr_zen_ang, $
	lat=lat, $
	lon=lon, $
	dbg=dbg

; Example usage:
; ARESE SGP CART site TST noon 951011_1200:
; arm2ltst,813435392,[0],lat=36.6048,lon=-97.4841,dbg=1
; ARESE SGP CART site TST 9AM 951011_0900:
; arm2ltst,813424576,[0],lat=36.6048,lon=-97.4841,dbg=1

if n_elements(lat) eq 0 then lat_dgr_north=36.6048 else lat_dgr_north=lat
if n_elements(lon_dgr_east) eq 0 then lon_dgr_east=-97.4841 else lon_dgr_east=lon
if n_elements(dbg) eq 0 then dbg=0

	sec_since_700101=double(base_time)+double(time_offset)
; pro compute_solar, time, lat, lon, coszen, solar_flux, sol_azim, day, hour
	compute_solar,sec_since_700101,lat_dgr_north,lon_dgr_east,cos_slr_zen_ang,flx_slr,azi_ang_slr,day_jul_int,hr_gmt,yr
	gmt_doy=day_jul_int+hr_gmt/24.
	gmm=2.*!pi*(gmt_doy-1.)/365.
	eqn_time_min=229.18*(.000075+.001868*cos(gmm)-.032077*sin(gmm)-.014615*cos(2.*gmm)-.040849*sin(2.*gmm))
	eqn_time_day=eqn_time_min/1440.
	gtst_doy=gmt_doy+eqn_time_day ; true solar day since 950101 in greenwich
	ltst_doy=gtst_doy+(lon_dgr_east/360.) ; local true solar time (gtst_doy corrected for longitude)

if dbg then begin
	print,'Quantities from compute_solar:'
	print,'cos_slr_zen_ang = ',cos_slr_zen_ang
	print,'slr_zen_ang_dgr = ',180*acos(cos_slr_zen_ang)/!pi
	print,'flx_slr = ',flx_slr
	print,'flx_TOA = ',flx_slr*cos_slr_zen_ang
	print,'year = ',yr
	print,'day_jul_int = ',day_jul_int	
	print,'hr_gmt = ',hr_gmt
	print,'gmt_doy = ',gmt_doy
	print,'gtst_doy = ',gtst_doy
	print,'ltst_doy = ',ltst_doy
endif; endif dbg
	
end; end arm2ltst()

;pro compute_solar, time, lat, lon, coszen, solar_flux, sol_azim
pro compute_solar, time, lat, lon, coszen, solar_flux, sol_azim, day, hour, year
;csz--

n_time = n_elements(time)
year = indgen(n_time)
day = indgen(n_time)
hour = (time / 3600.0) mod 24.0
;
; Calculate year and day.  Since these variables are slowly varying with time,
;   and since caldat and julday are slow routines, calculate year and day
;   only when the time changes by 24 hours (using the IDL uniq routine).
;
jan170_base = julday(1, 1, 1970)
juldat = jan170_base + fix(time / 86400.0)
uniq_day = uniq(juldat)
index0 = long(0)
for k = 0, n_elements(uniq_day) - 1 do begin
    caldat, juldat(uniq_day(k)), mm, dd, yy
    year(index0:uniq_day(k)) = yy
    dd = juldat(uniq_day(k)) - julday(1, 1, yy) + 1
    day(index0:uniq_day(k)) = dd
    index0 = long(uniq_day(k)) + long(1)
endfor
solar_flux = fltarr(n_time)
sunae1, year, day, hour, lat(0), lon(0), sol_azim, coszen, $
        flux = solar_flux

return
end
;---------------------------------------------------------------------------

;---------------------------------------------------------------------------

; Routine: sunae1.pro

; Purpose: calculates azimuth, cos(solar zenith)[elevation], 
;          and TOA flux of sun

; References:
; (a) Michalsky, J. J., 1988, The Astronomical Almanac's algorithm for
;     approximate solar position (1950-2050), Solar Energy, 227---235, 1988

; (b) Spencer, J. W., 1989, Comments on The Astronomical
;     Almanac's algorithm for approximate solar position (1950-2050)
;     Solar Energy, 42, 353

; Input:
; year - the year number (e.g. 1977)
; day  - the day number of the year starting with 1 for
;        January 1
; time - decimal time. E.g. 22.89 (8.30am eastern daylight time is
;        equal to 8.5+5(hours west of Greenwich) -1 (for daylight savings
;        time correction
; lat -  local latitude in degrees (north is positive)
; lon -  local longitude (east of Greenwich is positive. 
;                         i.e. Honolulu is 15.3, -157.8)

; Output:
; a - azimuth angle of the sun (measured east from north 0 to 360)
; e - elevation of the sun
; flux - TOA incident solar flux

; Spencer correction introduced and 3 lines of Michalsky code
; commented out (after calculation of az)

      pro sunae1,year,day,hour,lat,long,az,cos_zen,flux = flux
      twopi = double(2.*!PI) 
      rad = double(1/!RADEG) 
; get the current Julian date
      delta = double(year-1949.) 
      leap = double(fix(delta/4.)) 
      jd = double(32916.5+delta*365.+leap+day+hour/24.) 
; calculate ecliptic coordinates
      time = double(jd-51545.0) 
; force mean longitude between 0 and 360 degs
      mnlong = double(280.460+0.9856474*time) 
      mnlong = double(mnlong mod 360.) 
      mnlong = double(mnlong + 360. * (mnlong lt 0))
; mean anomaly in radians between 0, 2*!PI
      mnanom = double(357.528+0.9856003*time) 
      mnanom = double(mnanom mod 360.) 
      mnanom = double(mnanom + 360. * (mnanom lt 0.))
      mnanom = double(mnanom*rad) 
; compute ecliptic longitude and obliquity of ecliptic
      eclong = double(mnlong+1.915*sin(mnanom)+0.020*sin(2.*mnanom)) 
      eclong = double(eclong mod  360.) 
      eclong = double(eclong + 360. * (eclong lt 0))
      oblqec = double(23.429-0.0000004*time) 
      eclong = double(eclong*rad) 
      oblqec = double(oblqec*rad) 
; calculate right ascention and declination
      num = double(cos(oblqec)*sin(eclong)) 
      den = double(cos(eclong)) 
      ra = double(atan(num/den)) 
; force ra between 0 and 2*!PI
      ra = double(ra + !PI * (den lt 0.))
      ra = double(ra + twopi * (den ge 0. and num lt 0.))
; dec in radians
      dec = double(asin(sin(oblqec)*sin(eclong))) 
; calculate Greenwich mean sidereal time in hours
      gmst = double(6.697375+0.0657098242*time+hour) 
; hour not changed to sidereal sine "time" includes the fractional day
      gmst = double(gmst mod 24.) 
      gmst = double(gmst + 24. * (gmst lt 0.))
; calculate local mean sidereal time in radians
      lmst = double(gmst+long/15.) 
      lmst = double(lmst mod 24.) 
      lmst = double(lmst + 24. * (lmst lt 0.))
      lmst = double(lmst*15.*rad) 
; calculate hour angle in radians between -!PI, !PI
      ha  = double( lmst -ra) 
      ha = double(ha + twopi * (ha lt -!PI))
      ha = double(ha - twopi * (ha gt !PI)) 
      lat = double(lat*rad) 
; calculate azimuth and elevation
      el = double(asin(sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(ha))) 
      az = double(asin(-cos(dec)*sin(ha)/cos(el))) 
; this puts azimuth between 0 and 2*!PI radians
; add J. W. Spencer code (next 5 lines)
      azspenc = double(sin(dec) - sin(el)*sin(lat))
      azfix1 = where(azspenc  ge  0. and sin(az)  lt  0.)
      azfix2 = where(azspenc  lt  0.)
      if (azfix1(0) ne -1) then az(azfix1) = double(az(azfix1)+twopi) 
      if (azfix2(0) ne -1) then az(azfix2) = double(!PI-az(azfix2)) 
; end Spencer's corrections
;cc        elc = double(asin(sin(dec)/sin(lat))) 
;cc        if(el ge elc) az = double(!PI-az) 
;cc        if(el le elc  .and. ha gt 0.) az = double(twopi+az) 

; calculate refraction correction for US stand. atm.
      el = double(el/rad) 
      refrac = dblarr(n_elements(el))
      refrac(*) = 0.56
      i_r = where(el gt -0.56)
      if (i_r(0) ne -1) then begin
         refrac(i_r) = double(3.51561*(0.1594+0.0196*el(i_r)+$
                       0.00002*el(i_r)^2)/(1.+0.505*el(i_r)+$
                       0.0845*el(i_r)^2)) 
      endif
      cos_zen = double(sin((el+refrac)*rad))
      az = double(az/rad) 
      lat = double(lat/rad) 
      if keyword_set(flux) then begin
;
; Calculate flux using Michalsky's formula for earth-sun distance
;
; csz++
;         s0 = 1368.0
         s0 = 1367.0
; csz--
         R = 1.00014 - 0.01671 * cos(mnanom) - $
             0.00014 * cos(2.0 * mnanom)
         flux = double(s0 * (cos_zen gt 0.0) / R^2)
      endif
;
      return
      end
;---------------------------------------------------------------------------


