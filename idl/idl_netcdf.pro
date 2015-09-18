;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CVS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: idl_netcdf.pro,v $
; $Source: /home/zender/cvs/idl/idl_netcdf.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; Purpose: Serves as a template for full featured IDL programs.
;
; $Log: not supported by cvs2svn $
; Revision 1.4  2000-01-10 19:33:40  zender
; *** empty log message ***
;
; Revision 1.3  2000/01/01 01:55:52  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/03 16:52:04  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:50  zender
; Imported sources
;
; Revision 1.1  1993/07/06  01:04:15  zender
; Initial revision
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
pro nc, $
	filename=filename, $
	num_layers=num_layers, $
	pause=pause, $
	ps=ps
;
;Example usage:
;
; nc,filename='/data/zender/aca/nrel/950418_0000_dia_rsd.nc'
;
if n_elements(filename) eq 0 then filename = 'clouds.nc'
;
sys_time = systime(1) 
close,/all
erase
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read in binary data from a NetCDF file 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nc_id=ncdf_open(filename)
;
;now get the one-dimensional arrays
;
ncdf_varget,nc_id,'lev',lev
ncdf_varget,nc_id,'p',p
ncdf_varget,nc_id,'t',t
ncdf_varget,nc_id,'RH_liq',RH_liq
;
; say bye-bye
ncdf_close,nc_id
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of NetCDF commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot the Pressure vs. Altitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
plot,t-273.15,p/100.0, $
	tit='!5Temperature vs. Pressure', $
        xtit='!5Temperature !8T!5 (C)', $
        ytit='!5Pressure !8p!5 (mb)', $
        yrange=[1000.0,0.0], $
	/ynozero, $
        linestyle=0, $
	xstyle=0, $
	charsize=2.0, $
	ystyle=0, $
	thick=3.
;xyouts,.3,0.8,'DIA Sounding 950418 0000Z',size=1.75,/NORMAL

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot the Pressure vs. RH
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
plot,RH_liq*100.0,p/100.0, $
	tit='!5Relative Humidity vs. Pressure', $
        xtit='!5Relative Humidity !8RH!5 (%)', $
        ytit='!5Pressure !8p!5 (mb)', $
        yrange=[1000.0,0.0], $
	/ynozero, $
        linestyle=0, $
	xstyle=0, $
	charsize=2.0, $
	ystyle=0, $
	thick=3.
;xyouts,.3,0.8,'DIA Sounding 950418 0000Z',size=1.75,/NORMAL

;***************enforced exit********************************
goto,exit_gracefully

exit_gracefully: foo=1

end
