;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: cld_fgr.pro,v $
; $Source: /home/zender/cvs/idl/cld_fgr.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; NB: As a reminder when contouring, the number of contour levels and labels
; and annotations etc. should be one greater than the number of colors.
;
; Purpose: All the IDL figures for the shape and size sensitivity paper.
;
; $Log: not supported by cvs2svn $
; Revision 1.7  2000-01-15 02:07:49  zender
; *** empty log message ***
;
; Revision 1.4  2000/01/01 01:55:48  zender
; *** empty log message ***
;
; Revision 1.3  1999/12/31 20:12:45  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 02:09:36  zender
; *** empty log message ***
;
; Revision 1.1  1999/12/31 00:18:15  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/30 19:34:24  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:46  zender
; Imported sources
;
; Revision 2.0  1994/03/08  00:54:07  zender
;  the wonderful widgeting cloud data viewer
;
; Revision 1.8  1993/08/26  19:10:50  zender
; this is the version sent to the reviewers, fully double spaced,
; still has glitches on some of the color plots.
;
; Revision 1.7  1993/08/06  17:37:51  zender
; last version before california. text looks nice. minor problems
; with the color figures include lousy axis numbers and some
; color spills.
;
; Revision 1.6  1993/08/02  17:14:24  zender
; fixed the contour plotting error, i was dividing by crystal_length
; not the delta_length! changed from plotting the effective_radius
; (which was bogus), to finding and plotting effective_length
;
; Revision 1.5  1993/07/28  20:09:14  zender
; added text legends to line graphs
;
; Revision 1.4  1993/07/13  21:30:02  zender
; have mad many little adjustments
;
; Revision 1.3  1993/07/12  13:40:30  zender
; another all nighter, improved figure 6 slightly, looking
; at possible systematic error in contour plotting the
; log of the crystal distributions.
;
; Revision 1.2  1993/07/09  19:17:26  zender
; Have added many more figures. still no color plots.
;
; Revision 1.1  1993/07/05  23:31:07  zender
; Initial revision
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Dens Temp Pres
; Plot the environmental density, temperature, and pressure as a function
; of height
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro dens_temp_pres
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"dz",dz
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"env_density",env_density
ncdf_varget,cdfid,"env_pressure",env_pressure
ncdf_varget,cdfid,"env_temperature",env_temperature
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
axz_skip_nrm=0.15
num_y_axz=3
;
x_min=num_y_axz*axz_skip_nrm
y_min=0.11
x_max=0.95
y_max=0.90
plt_rgn_nrm=[x_min,y_min,x_max,y_max]
;
title_title='!7q!8(z)!5, !8p(z)!5, and !8T(z)!5'
xyouts,.5,0.95,title_title,size=2.5,alignment=0.5,/NORMAL
;
plot, $
	altitude(1:num_layer)/1000.0, $
	env_density(1:num_layer), $
	xtit='Altitude !8z!5 (km)', $
	ytit='Density !7q!8(z)!5 (kg-m!E-3!N)', $
	xstyle=1, $
	ystyle=0, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	position=plt_rgn_nrm, $
	linestyle=0, $
	/noerase
;
plot, $
	altitude(1:num_layer)/1000.0, $
	env_pressure(1:num_layer)/1.0e2, $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	linestyle=1, $
	position=plt_rgn_nrm, $
	/noerase
;
axis, $
	x_min-1.*axz_skip_nrm, $
	y_min, $
	/normal, $
	yaxis=0, $
	ystyle=0, $
	/ynozero, $
	ytitle='!5Pressure !8p(z)!5 (mb)', $
	charsize=2.
;
plot, $
	altitude(1:num_layer)/1000.0, $
	env_temperature(1:num_layer)-273.15, $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	linestyle=2, $
	position=plt_rgn_nrm, $
	/noerase
;
axis, $
	x_min-2.*axz_skip_nrm, $
	y_min, $
	/normal, $
	yaxis=0, $
	ystyle=0, $
	/ynozero, $
	ytit='Temperature !8T(z)!5 (!E!12_!N!5C)', $
	charsize=2.
;
ln_lgn_x1=.7
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.5
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!7q!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8p!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!8T!5',size=txt_lgn_sz,/NORMAL
;
end; end dens_temp_pres()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Dens Temp Pres commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Heating Rate Sensitivity.
; Plot the (primary - secondary)/abs(primary) heating rate profile 
; (i.e., the size or shape sensitivities)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro heat_rate_sns
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename_1=primary_fl
filename_2=secondary_fl
cdfid_1=ncdf_open(filename_1)
cdfid_2=ncdf_open(filename_2)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_layer",num_layer_1
ncdf_varget1,cdfid_2,"num_layer",num_layer_2
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"altitude",altitude_1
ncdf_varget,cdfid_1,"heating_rate_LW",LW_heating_1
ncdf_varget,cdfid_1,"heating_rate_SW",SW_heating_1
ncdf_varget,cdfid_1,"heating_rate_net",net_heating_1
ncdf_varget,cdfid_2,"altitude",altitude_2
ncdf_varget,cdfid_2,"heating_rate_LW",LW_heating_2
ncdf_varget,cdfid_2,"heating_rate_SW",SW_heating_2
ncdf_varget,cdfid_2,"heating_rate_net",net_heating_2
;
; say bye-bye
ncdf_close,cdfid_1
ncdf_close,cdfid_2
; End of NetCDF commands
;
; Compute the difference due to the shape effect
net_2_1_compare=altitude_2
LW_2_1_compare=altitude_2
SW_2_1_compare=altitude_2
epsilon=0.1/3600.
;
net_difference=net_heating_2-net_heating_1
net_2_1_compare=100.*net_difference*abs(net_heating_1)/ $
	(net_heating_1*net_heating_1+epsilon*epsilon)
LW_difference=LW_heating_2-LW_heating_1
LW_2_1_compare=100.*LW_difference*abs(LW_heating_1)/ $
	(LW_heating_1*LW_heating_1+epsilon*epsilon)
SW_difference=SW_heating_2-SW_heating_1
SW_2_1_compare=100.*SW_difference*abs(SW_heating_1)/ $
	(SW_heating_1*SW_heating_1+epsilon*epsilon)
;
;net_2_1_compare=net_difference*3600.
;LW_2_1_compare=LW_difference*3600.
;SW_2_1_compare=SW_difference*3600.
;
; Allow the y-axis to be within 10% of the curves
data_min=min([min(LW_2_1_compare),min(SW_2_1_compare),min(net_2_1_compare)])
data_max=max([max(LW_2_1_compare),max(SW_2_1_compare),max(net_2_1_compare)])
data_min=(data_min-10.) + abs(data_min mod 10.)
data_max=(data_max+10.) - abs(data_max mod 10.)
;
plot, $
	net_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	tit='!5Sensitivity of Heating Rates', $
;	xtit='Percent Change = 100*(Secondary - Primary)/!10"!5Primary!10"!5', $
	xtit='Percent Change (%)', $
;	xtit='Secondary - Primary (K hr!E-1!N)', $
	ytit='Altitude !8z!5 (km)', $
	xrange=[data_min,data_max], $
;	yrange=[0.1,0.6], $
	xmargin=[5.5,1.2], $	;[10,3] is [left,right] default
	ymargin=[3,1.5], $	;[4,2] is [bottom,top] default
	xstyle=1, $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	SW_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	thick=3.0, $
	linestyle=1
;
oplot,	$
	LW_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	thick=1.0, $
	linestyle=2
;
ln_lgn_x1=.3
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.5
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5SW',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5LW',size=txt_lgn_sz,/NORMAL
;
end; end heat_rate_sns()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Heating Rate Sensitivity.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Multi Radiative Forcing
; Plot the time evolution of radiative forcing for the two main cases.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro multi_rad_forcing
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename_1=primary_fl
cdfid_1=ncdf_open(filename_1)
filename_2=secondary_fl
cdfid_2=ncdf_open(filename_2)
filename_3=tertiary_fl
cdfid_3=ncdf_open(filename_3)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
ncdf_varget1,cdfid_2,"num_step",num_step_2
ncdf_varget1,cdfid_2,"rad_step",rad_step_2
ncdf_varget1,cdfid_3,"num_step",num_step_3
ncdf_varget1,cdfid_3,"rad_step",rad_step_3
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"LW_cloud_forcing",LW_forcing_1
ncdf_varget,cdfid_1,"SW_cloud_forcing",SW_forcing_1
ncdf_varget,cdfid_1,"time_array",time_array_1
ncdf_varget,cdfid_2,"LW_cloud_forcing",LW_forcing_2
ncdf_varget,cdfid_2,"SW_cloud_forcing",SW_forcing_2
ncdf_varget,cdfid_2,"time_array",time_array_2
ncdf_varget,cdfid_3,"LW_cloud_forcing",LW_forcing_3
ncdf_varget,cdfid_3,"SW_cloud_forcing",SW_forcing_3
ncdf_varget,cdfid_3,"time_array",time_array_3
;
; say bye-bye
ncdf_close,cdfid_1
ncdf_close,cdfid_2
ncdf_close,cdfid_3
; End of NetCDF commands
;
if n_elements(time_array_1) le 1 then begin
	print,"Requested data not in NetCDF file, returning."
	return
endif
;
; Convert to Diurnally Averaged forcings...
; The 1/pi factor is the diurnal average, and the 2/sqrt(3) factor converts
; the 30 degree zenith angle simulations into pseudo local-noon quantities.
SW_forcing_1=SW_forcing_1*2./(sqrt(3.)*!pi)
SW_forcing_2=SW_forcing_2*2./(sqrt(3.)*!pi)
SW_forcing_3=SW_forcing_3*2./(sqrt(3.)*!pi)
;
; ...and compute net forcings
net_forcing_1=SW_forcing_1+LW_forcing_1
net_forcing_2=SW_forcing_2+LW_forcing_2
net_forcing_3=SW_forcing_3+LW_forcing_3
;
plot, $
	time_array_1(0:num_step_1-1)/60.0, $
	net_forcing_1(0:num_step_1-1), $
	tit='!5Diurnally Averaged Cloud Forcing', $
	xtit='Elapsed Time (min.)', $
	ytit='Forcing (W-m!E-2!N)', $
	yrange=[-150.0,250.0], $
	xstyle=1, $
	ystyle=1, $
	thick=3.0, $
	xmargin=[7.8,.5], $	;[10,3] is [left,right] default
	ymargin=[3,1.5], $	;[4,2] is [bottom,top] default
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	time_array_2(0:num_step_2-1)/60.0, $
	net_forcing_2(0:num_step_2-1), $
	thick=3.0, $
	linestyle=0
;
;oplot,	$
;	time_array_3(0:num_step_3-1)/60.0, $
;	net_forcing_3(0:num_step_3-1), $
;	thick=3.0, $
;	linestyle=0
;
; Now overplot the SW forcings
oplot,	$
	time_array_1(0:num_step_1-1)/60.0, $
	SW_forcing_1(0:num_step_1-1), $
	thick=3.0, $
	linestyle=1
;
oplot,	$
	time_array_2(0:num_step_2-1)/60.0, $
	SW_forcing_2(0:num_step_2-1), $
	thick=3.0, $
	linestyle=1
;
;oplot,	$
;	time_array_3(0:num_step_3-1)/60.0, $
;	SW_forcing_3(0:num_step_3-1), $
;	thick=3.0, $
;	linestyle=1
;
; Now overplot the LW forcings
oplot,	$
	time_array_1(0:num_step_1-1)/60.0, $
	LW_forcing_1(0:num_step_1-1), $
	thick=3.0, $
	linestyle=2
;
oplot,	$
	time_array_2(0:num_step_2-1)/60.0, $
	LW_forcing_2(0:num_step_2-1), $
	thick=3.0, $
	linestyle=2
;
;oplot,	$
;	time_array_3(0:num_step_3-1)/60.0, $
;	LW_forcing_3(0:num_step_3-1), $
;	thick=3.0, $
;	linestyle=2
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.5
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
;xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(1),'!5SWCF',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(2),'!5LWCF',size=txt_lgn_sz,/NORMAL
;
txt_lgn_x=10
lgn_y=[125,35,75, $
	-55,-126,-120, $
	200,165]
;
xyouts,txt_lgn_x,lgn_y(0),'!5Truncated NetCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
;xyouts,txt_lgn_x,lgn_y(1),'!5Spherical NetCF',size=txt_lgn_sz,/DATA, $
;	ORIENTATION=-1
xyouts,txt_lgn_x,lgn_y(2),'!5Control NetCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1 
;
xyouts,txt_lgn_x,lgn_y(3),'!5Truncated SWCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
;xyouts,txt_lgn_x,lgn_y(4),'!5Spherical SWCF',size=txt_lgn_sz,/DATA, $
;	ORIENTATION=-1
xyouts,txt_lgn_x,lgn_y(5),'!5Control SWCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
;
xyouts,txt_lgn_x,lgn_y(6),'!5Control LWCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
xyouts,txt_lgn_x,lgn_y(7),'!5Truncated LWCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
;
end; end multi_rad_forcing()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Multi Radiative Forcing commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Multi Albedo Emissivity.
; Plot the comparison of albedo and emissivity for the three main cases.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro multi_alb_emis
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename_1=primary_fl
cdfid_1=ncdf_open(filename_1)
filename_2=secondary_fl
cdfid_2=ncdf_open(filename_2)
filename_3=tertiary_fl
cdfid_3=ncdf_open(filename_3)
filename_4=quaternary_fl
cdfid_4=ncdf_open(filename_4)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
ncdf_varget1,cdfid_2,"num_step",num_step_2
ncdf_varget1,cdfid_2,"rad_step",rad_step_2
ncdf_varget1,cdfid_3,"num_step",num_step_3
ncdf_varget1,cdfid_3,"rad_step",rad_step_3
ncdf_varget1,cdfid_4,"num_step",num_step_4
ncdf_varget1,cdfid_4,"rad_step",rad_step_4
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"albedo_of_time",albedo_1
ncdf_varget,cdfid_1,"emissivity_of_time",emissivity_1
ncdf_varget,cdfid_2,"albedo_of_time",albedo_2
ncdf_varget,cdfid_2,"emissivity_of_time",emissivity_2
ncdf_varget,cdfid_3,"albedo_of_time",albedo_3
ncdf_varget,cdfid_3,"emissivity_of_time",emissivity_3
ncdf_varget,cdfid_4,"albedo_of_time",albedo_4
ncdf_varget,cdfid_4,"emissivity_of_time",emissivity_4
ncdf_varget,cdfid_1,"time_array",time_array
;
; say bye-bye
ncdf_close,cdfid_1
ncdf_close,cdfid_2
ncdf_close,cdfid_3
ncdf_close,cdfid_4
; End of NetCDF commands
;
if n_elements(time_array) le 1 then begin
	print,"Requested data not in NetCDF file, returning."
	return
endif
;
plot, $
	emissivity_1(0:num_step_1-1), $
	albedo_1(0:num_step_1-1), $
	tit='!5Albedo vs. Emissivity', $
	xtit='Emissivity !7e!5', $
	ytit='Albedo !8A!5', $
	xrange=[.94,0.99], $
	xtickformat='(F4.2)', $
	yrange=[.20,.4], $
	xmargin=[6.5,1.5], $	;[10,3] is [left,right] default
	ymargin=[3,1.5], $	;[4,2] is [bottom,top] default
	xstyle=0, $
	ystyle=0, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0

oplot,	$
	emissivity_2(0:num_step_2-1), $
	albedo_2(0:num_step_2-1), $
	thick=3.0, $
	linestyle=1
;
oplot,	$
	emissivity_3(0:num_step_3-1), $
	albedo_3(0:num_step_3-1), $
	thick=3.0, $
	linestyle=2
;
oplot,	$
	emissivity_4(0:num_step_4-1), $
	albedo_4(0:num_step_4-1), $
	thick=3.0, $
	linestyle=3
;
ln_lgn_x1=.25
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.85
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Control Cloud',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5Truncated',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5Spherical Cloud',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(3),'!5Spherical Truncated',size=txt_lgn_sz,/NORMAL
;
end; end multi_alb_emis()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Multi Albedo Emissivity commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Humidity Level
; Plot the saturations and temperature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro humidity
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"saturation_ice",saturation_ice
ncdf_varget,cdfid,"saturation_liquid",saturation_liquid
ncdf_varget,cdfid,"env_temperature",env_temperature
ncdf_varget,cdfid,"env_pressure",env_pressure
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
plt_rgn_nrm=[ $
	.1, $ ; x_min
	.11, $ ; y_min
	.90, $ ; x_max
	.90] ; y_max
;
plot, $
	saturation_ice(1:num_layer)*100.0, $
	altitude(1:num_layer)/1000.0, $
;	tit='!5 ', $
	xtit='!5Saturation !8S!5 (%)', $
	ytit='!5Altitude !8z!5 (km)', $
	xrange=[20.0,120.0], $
;	yrange=[8.,18.], $
	xstyle=8, $
	ystyle=8, $
	/ynozero, $
	thick=3.0, $
	ytick_get=y_tick_coords, $
	charsize=2.0, $
;	position=plt_rgn_nrm, $
	xmargin=[5.4,5.6], $ ;[10,3] is [left,right] default
	ymargin=[3,2.3], $  ;[4,2] is [bottom,top] default
	linestyle=0
;
oplot,	$
	saturation_liquid(1:num_layer)*100.0, $
	altitude(1:num_layer)/1000.0, $
	thick=3.0, $
	linestyle=1
;
plot, $
	env_temperature(1:num_layer)-273.15, $
	altitude(1:num_layer)/1000.0, $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=2, $
;	position=plt_rgn_nrm, $
	xmargin=[5.4,5.6], $ ;[10,3] is [left,right] default
	ymargin=[3,2.3], $  ;[4,2] is [bottom,top] default
	/noerase
;
axis, $
	xaxis=1, $
	xtitle='!5Temperature !8T!5 (!E!12_!N!5C)', $
	xmargin=[5.4,5.6], $ ;[10,3] is [left,right] default
	charsize=2.
;
num_y_ticks=n_elements(y_tick_coords)
;print,num_y_ticks,y_tick_coords,' km'
ystring=strarr(num_y_ticks)
for tick=0,num_y_ticks-1 do begin
	foo=min(abs(altitude(1:num_layer)/1000.-y_tick_coords(tick)),altitude_idx)
	real_idx=altitude_idx+1
	if real_idx lt num_layer then begin	
;	do a forward interpolation
	y_interp=env_pressure(real_idx)+ $
		(1000.*y_tick_coords(tick)-altitude(real_idx))* $		
		((env_pressure(real_idx)-env_pressure(real_idx+1))/ $
		(altitude(real_idx)-altitude(real_idx+1)))
	endif else begin
;	do a backward interpolation
	y_interp=env_pressure(real_idx)+ $
		(1000.*y_tick_coords(tick)-altitude(real_idx))* $	
		((env_pressure(real_idx-1)-env_pressure(real_idx))/ $
		(altitude(real_idx-1)-altitude(real_idx)))
	endelse
	ystring(tick)=string(format='(I3)',y_interp/100.)
endfor
;print,num_y_ticks,ystring,' mb'
;
axis, $
	yaxis=1, $
	yticks=num_y_ticks-1, $
;	yminor=4, $
	ytickn=ystring, $
	ytitle='!5Pressure !8P!5 (mb)', $
	/ytype, $ ;this sets the axis to logarithmic
	ymargin=[3,2.3], $  ;[4,2] is [bottom,top] default
	charsize=2.
;
ln_lgn_x1=.45
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.85
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8S!I!5ice!N',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8S!I!5liquid!N',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!8T!5',size=txt_lgn_sz,/NORMAL
;
end; humidity()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Humidity Level commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Mass Mixing Ratio
; Plot the saturations and temperature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mass_mix_ratio
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"IWC",IWC
ncdf_varget,cdfid,"vapor_density",vapor_density
ncdf_varget,cdfid,"env_temperature",env_temperature
ncdf_varget,cdfid,"env_pressure",env_pressure
ncdf_varget,cdfid,"env_density",env_density
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
plt_rgn_nrm=[ $
	.1, $ ; x_min
	.11, $ ; y_min
	.90, $ ; x_max
	.90] ; y_max
;
qice=[0,IWC(1:num_layer)/env_density(1:num_layer),0]
qvapor=[0,vapor_density(1:num_layer)/env_density(1:num_layer),0]
;
plot, $
	qice(1:num_layer)*1000.0, $
	altitude(1:num_layer)/1000.0, $
;	tit='!5 ', $
	xtit='!5Specific Humidity !8q!5 (g/kg)', $
	ytit='!5Altitude !8z!5 (km)', $
;	xrange=[20.0,120.0], $
;	yrange=[8.,18.], $
	xstyle=8, $
	ystyle=8, $
	/ynozero, $
	thick=3.0, $
	ytick_get=y_tick_coords, $
	charsize=2.0, $
;	position=plt_rgn_nrm, $
	xmargin=[7,7], $ ;[10,3] is [left,right] default
	ymargin=[3,3], $  ;[4,2] is [bottom,top] default
	linestyle=0
;
oplot,	$
	qvapor(1:num_layer)*1000.0, $
	altitude(1:num_layer)/1000.0, $
	thick=3.0, $
	linestyle=1
;
plot, $
	env_temperature(1:num_layer)-273.15, $
	altitude(1:num_layer)/1000.0, $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=2, $
;	position=plt_rgn_nrm, $
	xmargin=[7,7], $ ;[10,3] is [left,right] default
	ymargin=[3,3], $  ;[4,2] is [bottom,top] default
	/noerase
;
axis, $
	xaxis=1, $
	xtitle='!5Temperature !8T!5 (!E!12_!N!5C)', $
	charsize=2.
;
num_y_ticks=n_elements(y_tick_coords)
;print,num_y_ticks,y_tick_coords,' km'
ystring=strarr(num_y_ticks)
for tick=0,num_y_ticks-1 do begin
	foo=min(abs(altitude(1:num_layer)/1000.-y_tick_coords(tick)),altitude_idx)
	real_idx=altitude_idx+1
	if real_idx lt num_layer then begin	
;	do a forward interpolation
	y_interp=env_pressure(real_idx)+ $
		(1000.*y_tick_coords(tick)-altitude(real_idx))* $		
		((env_pressure(real_idx)-env_pressure(real_idx+1))/ $
		(altitude(real_idx)-altitude(real_idx+1)))
	endif else begin
;	do a backward interpolation
	y_interp=env_pressure(real_idx)+ $
		(1000.*y_tick_coords(tick)-altitude(real_idx))* $	
		((env_pressure(real_idx-1)-env_pressure(real_idx))/ $
		(altitude(real_idx-1)-altitude(real_idx)))
	endelse
	ystring(tick)=string(format='(I3)',y_interp/100.)
endfor
;print,num_y_ticks,ystring,' mb'
;
axis, $
	yaxis=1, $
	yticks=num_y_ticks-1, $
;	yminor=4, $
	ytickn=ystring, $
	ytitle='!5Pressure !8P!5 (mb)', $
	/ytype, $ ;this sets the axis to logarithmic
	charsize=2.
;
ln_lgn_x1=0.50
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8q!I!5ice!N',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8q!I!5vapor!N',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!8T!5',size=txt_lgn_sz,/NORMAL
;
end; mass_mix_ratio()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Mass Mixing Ratio commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Local IWC
; Plot the local IWC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro local_IWC
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"dz",dz
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"IWC",IWC
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
IWP_above_here=IWC
IWP_below_here=IWC
for layer=1,num_layer do begin
	IWP_above_here(layer)=total(IWC(layer:num_layer))*dz
	IWP_below_here(layer)=total(IWC(1:layer))*dz
endfor
;
axz_skip_nrm=0.15
num_y_axz=3
;
x_min=num_y_axz*axz_skip_nrm
y_min=0.11
x_max=0.95
y_max=0.90
plt_rgn_nrm=[x_min,y_min,x_max,y_max]
;
title_title='!8IWC(z)!5, !8IWP > z!5, and !8IWP < z!5'
xyouts,.5,0.95,title_title,size=2.5,alignment=0.5,/NORMAL
;
plot, $
	altitude(1:num_layer)/1000.0, $
	IWC(1:num_layer)*1.0e6, $
	xtit='Altitude !8z!5 (km)', $
	ytit='Ice Water Content !8IWC!5 (mg-m!E-3!N)', $
	xstyle=1, $
	ystyle=0, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	position=plt_rgn_nrm, $
	linestyle=0, $
	/noerase
;
plot, $
	altitude(1:num_layer)/1000.0, $
	IWP_above_here(1:num_layer)*1.0e3, $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	linestyle=1, $
	position=plt_rgn_nrm, $
	/noerase
;
axis, $
	x_min-1.*axz_skip_nrm, $
	y_min, $
	/normal, $
	yaxis=0, $
	ystyle=0, $
	/ynozero, $
	ytit='Ice Water Path Above !8IWP > z!5 (g-m!E-2!N)', $
	charsize=2.
;
plot, $
	altitude(1:num_layer)/1000.0, $
	IWP_below_here(1:num_layer)*1.0e3, $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	linestyle=2, $
	position=plt_rgn_nrm, $
	/noerase
;
axis, $
	x_min-2.*axz_skip_nrm, $
	y_min, $
	/normal, $
	yaxis=0, $
	ystyle=0, $
	/ynozero, $
	ytit='Ice Water Path Below !8IWP < z!5 (g-m!E-2!N)', $
	charsize=2.
;
ln_lgn_x1=.7
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.5
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8IWC!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8IWP > z!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!8IWP < z!5',size=txt_lgn_sz,/NORMAL
;
end; end local_IWC()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Local IWC commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure IWP Reff Albedo Emissivity commands
; Plot the IWP and effective radius as a function of time.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro IWP_reff_alb_em
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename_1=primary_fl
cdfid_1=ncdf_open(filename_1)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"IWP_of_time",IWP_1
ncdf_varget,cdfid_1,"albedo_of_time",albedo_of_time_1
ncdf_varget,cdfid_1,"emissivity_of_time",emissivity_of_time_1
;ncdf_varget,cdfid_1,"saturation_ice_of_time",saturation_ice_of_time_1
ncdf_varget,cdfid_1,"cloud_effective_radius",cloud_effective_radius_1
ncdf_varget,cdfid_1,"time_array",time_array
;
; say bye-bye
ncdf_close,cdfid_1
; End of NetCDF commands
;
if n_elements(time_array) le 1 then begin
	print,"Requested data not in NetCDF file, returning."
	return
endif
;
axz_skip_nrm=0.15
num_y_axz=4
;
x_min=num_y_axz*axz_skip_nrm
y_min=0.11
x_max=0.95
y_max=0.90
plt_rgn_nrm=[x_min,y_min,x_max,y_max]
;
title_title='!8IWP(t)!5, !8r!5!Ie!N!8(t)!5, !8A(t)!5, and !7e!8(t)!5'
xyouts,.5,0.95,title_title,size=2.5,alignment=0.5,/NORMAL
;
plot, $
	time_array(0:num_step_1)/3600.0, $
	IWP_1(0:num_step_1)*1000.0, $
	xtit='Time !8t!5 (hrs.)', $
	ytit='Ice Water Path !8IWP!5 (g-m!E-2!N)', $
	xstyle=1, $
	ystyle=0, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	position=plt_rgn_nrm, $
	linestyle=0, $
	/noerase
;
plot, $
	time_array(0:num_step_1)/3600.0, $
	cloud_effective_radius_1(0:num_step_1)*1.0e6, $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	linestyle=1, $
	position=plt_rgn_nrm, $
	/noerase
;
axis, $
	x_min-1.*axz_skip_nrm, $
	y_min, $
	/normal, $
	yaxis=0, $
	ystyle=0, $
	/ynozero, $
	ytitle='!5Effective Radius !8r!5!Ie!N (!7l!5m)', $
	charsize=2.
;
plot, $
	time_array(0:num_step_1-1)/3600.0, $
	albedo_of_time_1(0:num_step_1-1), $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	linestyle=2, $
	position=plt_rgn_nrm, $
	/noerase
;
axis, $
	x_min-2.*axz_skip_nrm, $
	y_min, $
	/normal, $
	yaxis=0, $
	ystyle=0, $
	/ynozero, $
	ytitle='!5Albedo !8A!5', $
	charsize=2.
;
plot, $
	time_array(0:num_step_1-1)/3600.0, $
	emissivity_of_time_1(0:num_step_1-1), $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	linestyle=3, $
	position=plt_rgn_nrm, $
	/noerase
;
axis, $
	x_min-3.*axz_skip_nrm, $
	y_min, $
	/normal, $
	yaxis=0, $
	ystyle=0, $
	/ynozero, $
	ytitle='!5Emissivity !7e!5', $
	charsize=2.
;
ln_lgn_x1=.7
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.5
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8IWP!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8r!5!Ie!N',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!8A!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(3),'!7e!5',size=txt_lgn_sz,/NORMAL
;
end; end IWP_reff_alb_em()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure IWP Reff Albedo Emissivity commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Mass Distribution Color
; Contour plot the final Mass Distribution vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_size] (in C)
; is accessed as         foo(1..num_size,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mass_dst_color
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_size",num_size
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"IWC",IWC
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"crystal_mass",crystal_mass
ncdf_varget,cdfid,"effective_radius",effective_radius
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
data=concentration ; #/m^3
for layer=1,num_layer do begin
	data(1:num_size,layer)=concentration(1:num_size,layer)/ $
				delta_length(1:num_size)	; #/m^3/m
	data(1:num_size,layer)=data(1:num_size,layer)* $
				crystal_mass(1:num_size)	; kg/m^3/m
endfor
;
data=data(1:num_size,1:num_layer) ; kg/m^3/m
data=data*1000. ;g/m^3/m
data=data/1.0e6  ;g/m^3/um
data=data*1.0e3  ;mg/m^3/um
abc=crystal_length(1:num_size)*1.0e6
ord=altitude(1:num_layer)/1000.
title_title='!5Mass Distribution !8m(L,z)!5 (mg-m!E-3!N-!7l!5m!E-1!N)'
x_axz_title='Crystal length !8L!5 (!7l!5m)'
y_axz_title='Altitude !8z!5 (km)'
abc_min=0.
abc_max=650.
ord_min=14.5
ord_max=17.5
;
cntr_lvl_nbr=22
cntr_ntv=2.
cntr_lvl_min=1.0e-6
;
data_min=min(data)
data_max=max(data)
cntr_lvls=fltarr(cntr_lvl_nbr)
cntr_lvls(0)=0.
cntr_lvls(1)=cntr_lvl_min
for i=2,cntr_lvl_nbr-1 do begin
	cntr_lvls(i)=cntr_lvls(i-1)*cntr_ntv
endfor
cntr_which_lbl=indgen(cntr_lvl_nbr)*0	;which levels to label
cntr_ntt=string(format='(E7.1)',cntr_lvls) ;actual text labels
cbar_lgn=cntr_ntt
;
print,"maximum data = ",data_max
print,"maximum contour = ",cntr_lvls(cntr_lvl_nbr-1)
;
;cntr_lvls=dist_cntr_lvls
;cntr_lvl_nbr=n_elements(cntr_lvls)
;cntr_which_lbl=indgen(cntr_lvl_nbr)*0	;which levels to label
;cntr_ntt=dist_cntr_lbls			;actual text labels
;cbar_lgn=dist_cntr_lbls
;
num_cbar_colors=min([!d.table_size-1,cntr_lvl_nbr-1])
cntr_fll_idx=indgen(num_cbar_colors)+2
cbar_idx=indgen(num_cbar_colors)+2
;
cbar_fmt='(A5)'
cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=1.3
cbar_unit=""
cbar_lbl_sz=1.5
;
cbar_psn=[ $
	.88, $ ; x_min
	.10, $ ; y_min
	.92, $ ; x_max
	.90] ; y_max
;
plt_rgn_nrm=[ $
	.1, $ ; x_min
	.1, $ ; y_min
	.87, $ ; x_max
	.90] ; y_max
;
plt_rgn_dvc=convert_coord( $
	[plt_rgn_nrm(0), $
	plt_rgn_nrm(2)], $
	[plt_rgn_nrm(1), $
	plt_rgn_nrm(3)], $
	/normal, $
	/to_device)
;
plt_rgn_dvc=[ $
	plt_rgn_dvc(0), $
	plt_rgn_dvc(1), $
	plt_rgn_dvc(3), $
	plt_rgn_dvc(4)]
;
clr_mk,num_cbar_colors,color_order,pll_idx,1
;
tvlct,r_curr,g_curr,b_curr
;
clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=num_cbar_colors, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz
;
title_sz=1.6
;
contour, $
	data, $
	abc, $
	ord, $
	tit=title_title, $
	xtit=x_axz_title, $
	ytit=y_axz_title, $
	level=cntr_lvls, $                 
	c_color=cntr_fll_idx, $		
	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
;	font=0, $				
	xrange=[abc_min,abc_max], $
;	yrange=[ord_min,ord_max], $
	xstyle=1, $
	ystyle=1, $
	charsize=title_sz, $
	position=plt_rgn_nrm, $
	ycharsize=1.1, $
	xcharsize=1.1, $
	/closed, $				
	/fill, $				
	/noerase
;
contour, $
	data, $
	abc, $
	ord, $
	level=cntr_lvls, $                 
	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
	/closed, $				
	/overplot
;
effective_length=effective_radius
for layer=1,num_layer do begin
	effective_length(layer)=equiv_rad_to_hex_clm_length(effective_radius(layer))
endfor
;
cloud_eff_length=equiv_rad_to_hex_clm_length(cloud_eff_radius(num_step))
;
in_cloud=where(IWC gt .5e-6)
;
oplot,	$
	effective_length(in_cloud)*1.0e6, $
	altitude(in_cloud)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
	[cloud_eff_length,cloud_eff_length]*1.0e6, $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
end; end mass_dst_color()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Mass Distribution Color commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Mass Distribution B&W
; Contour plot the final Mass Distribution vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_size] (in C)
; is accessed as         foo(1..num_size,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mass_dst_bw
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_size",num_size
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"IWC",IWC
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"crystal_mass",crystal_mass
ncdf_varget,cdfid,"effective_radius",effective_radius
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
data=concentration ; #/m^3
for layer=1,num_layer do begin
	data(1:num_size,layer)=concentration(1:num_size,layer)/ $
				delta_length(1:num_size)	; #/m^3/m
	data(1:num_size,layer)=data(1:num_size,layer)* $
				crystal_mass(1:num_size)	; kg/m^3/m
endfor
;
data=data(1:num_size,1:num_layer) ; kg/m^3/m
data=data*1000. ;g/m^3/m
data=data/1.0e6  ;g/m^3/um
data=data*1.0e3  ;mg/m^3/um
abc=crystal_length(1:num_size)*1.0e6
ord=altitude(1:num_layer)/1000.
title_title='!5Mass Distribution !8m(L,z)!5 (mg-m!E-3!N-!7l!5m!E-1!N)'
x_axz_title='Crystal length !8L!5 (!7l!5m)'
y_axz_title='Altitude !8z!5 (km)'
abc_min=0.
abc_max=650.
ord_min=14.5
ord_max=17.5
;
cntr_lvl_nbr=22
cntr_ntv=2.
cntr_lvl_min=1.0e-6
;
data_min=min(data)
data_max=max(data)
cntr_lvls=fltarr(cntr_lvl_nbr)
cntr_lvls(0)=0.
cntr_lvls(1)=cntr_lvl_min
for i=2,cntr_lvl_nbr-1 do begin
	cntr_lvls(i)=cntr_lvls(i-1)*cntr_ntv
endfor
;
contour, $
	data, $
	abc, $
	ord, $
	tit=title_title, $
	xtit=x_axz_title, $
	ytit=y_axz_title, $
	level=cntr_lvls, $                 
	xrange=[abc_min,abc_max], $
;	yrange=[ord_min,ord_max], $
	xstyle=1, $
	ystyle=1, $
	/closed, $
	/noerase
;
effective_length=effective_radius
for layer=1,num_layer do begin
	effective_length(layer)=equiv_rad_to_hex_clm_length(effective_radius(layer))
endfor
;
cloud_eff_length=equiv_rad_to_hex_clm_length(cloud_eff_radius(num_step))
;
in_cloud=where(IWC gt .5e-6)
;
oplot,	$
	effective_length(in_cloud)*1.0e6, $
	altitude(in_cloud)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
	[cloud_eff_length,cloud_eff_length]*1.0e6, $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
end; mass_dst_bw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Mass Distribution B&W commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Fall Speed B&W.
; Contour plot the Fall Speed vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_size] (in C)
; is accessed as         foo(1..num_size,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fall_speed_bw
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_size",num_size
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"crystal_length",crystal_length
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"daltitude_dtime",daltitude_dtime
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
data=fltarr(num_size,num_layer)
for layer=1,num_layer do begin
	data(*,layer-1)=daltitude_dtime(layer,1:num_size)
endfor
data=-data*1.0e2
abc=crystal_length(1:num_size)*1.0e6
ord=altitude(1:num_layer)/1000.
title_title='!5Net Fall speed !8dz/dt!5 (cm-s!E-1!N)'
x_axz_title='Crystal length !8L!5 (!7l!5m)'
y_axz_title='Altitude !8z!5 (km)'
abc_min=0.
abc_max=650.
ord_min=14.5
ord_max=17.5
;
cntr_lvl_nbr=30
grd_ntv=5.
cntr_lvls=fltarr(cntr_lvl_nbr)
cntr_ln_sty=fltarr(cntr_lvl_nbr)
cntr_thk=fltarr(cntr_lvl_nbr)
;
data_min=min(data)
data_max=max(data)
cntr_lvl_min=(data_min-grd_ntv) + abs(data_max mod grd_ntv)
;max_level=(data_max+grd_ntv) - abs(data_max mod grd_ntv)
max_level=100.
cntr_ntv=(max_level-cntr_lvl_min)/(cntr_lvl_nbr-1)
cntr_lvls=cntr_lvl_min+findgen(cntr_lvl_nbr)*cntr_ntv
for i=0,cntr_lvl_nbr-1 do begin
	if cntr_lvls(i) lt 0. then begin
		cntr_ln_sty(i)=1 
		cntr_thk(i)=2 
	endif else begin
		cntr_ln_sty(i)=0 
		cntr_thk(i)=1 
	endelse
endfor
;
contour, $
	data, $
	abc, $
	ord, $
	tit=title_title, $
	xtit=x_axz_title, $
	ytit=y_axz_title, $
	level=cntr_lvls, $                 
	c_linestyle=cntr_ln_sty, $	
	c_thick=cntr_thk, $		
	xrange=[abc_min,abc_max], $
;	yrange=[ord_min,ord_max], $
	xstyle=1, $
	ystyle=1, $
	/closed					
;
end; end fall_speed_bw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Fall Speed B&W commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Size Distribution color.
; Contour plot the initial size distribution vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_size] (in C)
; is accessed as         foo(1..num_size,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sz_dst_color
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_size",num_size
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"IWC",IWC
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"effective_radius",effective_radius
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
distribution=concentration
for layer=1,num_layer do begin
	distribution(1:num_size,layer)=concentration(1:num_size,layer)/ $
		delta_length(1:num_size) ; #/m^3/m
endfor
;
data=distribution(1:num_size,1:num_layer)/1.0e9
abc=crystal_length(1:num_size)*1.0e6
ord=altitude(1:num_layer)/1000.
title_title='!5Size Distribution !8n(L,z)!5 (!8l!5!E-1!N-!7l!5m!E-1!N)'
x_axz_title='Crystal length !8L!5 (!7l!5m)'
y_axz_title='Altitude !8z!5 (km)'
abc_min=0.
abc_max=650.
ord_min=14.5
ord_max=17.5
x_map_cbar_mrg=[4,4]	;[10,3] is [left,right] default
y_map_cbar_mrg=[3,2]	;[4,2] is [bottom,top] default
;
cntr_lvl_nbr=25
cntr_ntv=2.
cntr_lvl_min=0.5e-3
;
data_min=min(data)
data_max=max(data)
cntr_lvls=fltarr(cntr_lvl_nbr)
cntr_lvls(0)=0.
cntr_lvls(1)=cntr_lvl_min
for i=2,cntr_lvl_nbr-1 do begin
	cntr_lvls(i)=cntr_lvls(i-1)*cntr_ntv
endfor
cntr_which_lbl=indgen(cntr_lvl_nbr)*0	;which levels to label
cntr_ntt=string(format='(E7.1)',cntr_lvls) ;actual text labels
cbar_lgn=cntr_ntt
;
print,"maximum data = ",data_max
print,"maximum contour = ",cntr_lvls(cntr_lvl_nbr-1)
;
;cntr_lvls=dist_cntr_lvls
;cntr_lvl_nbr=n_elements(cntr_lvls)
;cntr_which_lbl=indgen(cntr_lvl_nbr)*0	;which levels to label
;cntr_ntt=dist_cntr_lbls			;actual text labels
;cbar_lgn=dist_cntr_lbls
;
num_cbar_colors=min([!d.table_size-1,cntr_lvl_nbr-1])
cntr_fll_idx=indgen(num_cbar_colors)+2
cbar_idx=indgen(num_cbar_colors)+2
;
cbar_fmt='(A5)'
cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=1.3
cbar_unit=""
cbar_lbl_sz=1.5
;
cbar_psn=[ $
	.85, $ ; x_min
	.10, $ ; y_min
	.89, $ ; x_max
	.95] ; y_max
;
plt_rgn_nrm=[ $
	.11, $ ; x_min
	.1, $ ; y_min
	.84, $ ; x_max
	.95] ; y_max
;
plt_rgn_dvc=convert_coord( $
	[plt_rgn_nrm(0), $
	plt_rgn_nrm(2)], $
	[plt_rgn_nrm(1), $
	plt_rgn_nrm(3)], $
	/normal, $
	/to_device)
;
plt_rgn_dvc=[ $
	plt_rgn_dvc(0), $
	plt_rgn_dvc(1), $
	plt_rgn_dvc(3), $
	plt_rgn_dvc(4)]
;	
;
clr_mk,num_cbar_colors,color_order,pll_idx,1
;
tvlct,r_curr,g_curr,b_curr
;
clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=num_cbar_colors, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz
;
title_sz=1.6
;
contour, $
	data, $
	abc, $
	ord, $
	tit=title_title, $
	xtit=x_axz_title, $
	ytit=y_axz_title, $
	level=cntr_lvls, $                 
	c_color=cntr_fll_idx, $		
	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
;	font=0, $				
	xrange=[abc_min,abc_max], $
;	yrange=[ord_min,ord_max], $
	xstyle=1, $
	ystyle=1, $
	charsize=title_sz, $
	position=plt_rgn_nrm, $
	ycharsize=1.1, $
	xcharsize=1.1, $
	xmargin=x_map_cbar_mrg, $
	ymargin=y_map_cbar_mrg, $
	/closed, $				
	/fill, $				
	/noerase
;
contour, $
	data, $
	abc, $
	ord, $
	level=cntr_lvls, $                 
	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
	xmargin=x_map_cbar_mrg, $
	ymargin=y_map_cbar_mrg, $
	/closed, $				
	/overplot
;
effective_length=effective_radius
for layer=1,num_layer do begin
	effective_length(layer)=equiv_rad_to_hex_clm_length(effective_radius(layer))
endfor
;
cloud_eff_length=equiv_rad_to_hex_clm_length(cloud_eff_radius(num_step))
;
in_cloud=where(IWC gt .5e-6)
;
oplot,	$
	effective_length(in_cloud)*1.0e6, $
	altitude(in_cloud)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
	[cloud_eff_length,cloud_eff_length]*1.0e6, $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
end; end sz_dst_color()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Size Distribution color commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Size Distribution B&W
; Contour plot the final size distribution vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_size] (in C)
; is accessed as         foo(1..num_size,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro sz_dst_bw
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_size",num_size
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"IWC",IWC
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"effective_radius",effective_radius
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
distribution=concentration
for layer=1,num_layer do begin
	distribution(1:num_size,layer)=concentration(1:num_size,layer)/ $
		delta_length(1:num_size) ; #/m^3/m
endfor
;
data=distribution(1:num_size,1:num_layer)/1.0e9
abc=crystal_length(1:num_size)*1.0e6
ord=altitude(1:num_layer)/1000.
title_title='!5Size Distribution !8n(L,z)!5 (!8l!5!E-1!N-!7l!5m!E-1!N)'
x_axz_title='Crystal length !8L!5 (!7l!5m)'
y_axz_title='Altitude !8z!5 (km)'
abc_min=0.
abc_max=650.
ord_min=14.5
ord_max=17.5
;
cntr_lvl_nbr=25
cntr_ntv=2.
cntr_lvl_min=0.5e-3
;
cntr_lvls=fltarr(cntr_lvl_nbr)
cntr_lvls(0)=0.
cntr_lvls(1)=cntr_lvl_min
for i=2,cntr_lvl_nbr-1 do begin
	cntr_lvls(i)=cntr_lvls(i-1)*cntr_ntv
endfor
;
contour, $
	data, $
	abc, $
	ord, $
	tit=title_title, $
	xtit=x_axz_title, $
	ytit=y_axz_title, $
	level=cntr_lvls, $                 
	xrange=[abc_min,abc_max], $
;	yrange=[ord_min,ord_max], $
	xstyle=1, $
	ystyle=1, $
	/closed, $				
	/noerase
;
effective_length=effective_radius
for layer=1,num_layer do begin
	effective_length(layer)=equiv_rad_to_hex_clm_length(effective_radius(layer))
endfor
cloud_eff_length=equiv_rad_to_hex_clm_length(cloud_eff_radius(num_step))
;
in_cloud=where(IWC gt .5e-6)
;
oplot,	$
	effective_length(in_cloud)*1.0e6, $
	altitude(in_cloud)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
	[cloud_eff_length,cloud_eff_length]*1.0e6, $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
end; end sz_dst_bw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Size Distribution B&W commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure IWC color.
; Contour plot the IWC vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_size] (in C)
; is accessed as         foo(1..num_size,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro IWC_color
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_size",num_size
ncdf_varget1,cdfid,"num_step",num_step
ncdf_varget1,cdfid,"num_frame",num_frame
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"time_snapshot",time_snapshot
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"IWC_snapshot",IWC_snapshot
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
if n_elements(time_snapshot) le 1 then begin
	print,"Requested data not in NetCDF file, returning."
	return
endif
;
data=fltarr(num_frame+1,num_layer)
for layer=1,num_layer do begin
data(*,layer-1)=IWC_snapshot(layer,*)
endfor
data=data*1.0e6
abc=time_snapshot(*)/3600.
ord=altitude(1:num_layer)/1000.
title_title='!5Ice Water Content !8IWC(t,z)!5 (mg-m!E-3!N)'
x_axz_title='Time !8t!5 (hrs.)'
y_axz_title='Altitude !8z!5 (km)'
abc_min=0.
abc_max=650.
ord_min=14.5
ord_max=17.5
;
cntr_lvl_nbr=22
grd_ntv=5.
cntr_lvl_min=0.
first_level=0.5
;
data_min=min(data)
data_max=max(data)
max_level=(data_max+grd_ntv) - abs(data_max mod grd_ntv)
cntr_ntv=(max_level-cntr_lvl_min)/(cntr_lvl_nbr-2)
cntr_lvls=fltarr(cntr_lvl_nbr)
cntr_lvls(0)=cntr_lvl_min
if first_level lt data_min then first_level=cntr_ntv
cntr_lvls(1)=first_level
for i=2,cntr_lvl_nbr-1 do begin
	cntr_lvls(i)=cntr_lvls(i-1)+cntr_ntv
endfor
cntr_which_lbl=indgen(cntr_lvl_nbr)*0         ;which levels to label
cntr_ntt=string(format='(F4.1)',cntr_lvls) ;actual text labels
cbar_lgn=cntr_ntt
;
print,"maximum data = ",data_max
print,"maximum contour = ",cntr_lvls(cntr_lvl_nbr-1)
;
num_cbar_colors=min([!d.table_size-1,cntr_lvl_nbr-1])
cntr_fll_idx=indgen(num_cbar_colors)+2
cbar_idx=indgen(num_cbar_colors)+2
;
cbar_fmt='(A5)'
cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=1.3
cbar_unit=""
cbar_lbl_sz=1.5
;
cbar_psn=[ $
	.90, $ ; x_min
	.10, $ ; y_min
	.94, $ ; x_max
	.95] ; y_max
;
plt_rgn_nrm=[ $
	.11, $ ; x_min
	.1, $ ; y_min
	.89, $ ; x_max
	.95] ; y_max
;
plt_rgn_dvc=convert_coord( $
	[plt_rgn_nrm(0), $
	plt_rgn_nrm(2)], $
	[plt_rgn_nrm(1), $
	plt_rgn_nrm(3)], $
	/normal, $
	/to_device)
;
plt_rgn_dvc=[ $
	plt_rgn_dvc(0), $
	plt_rgn_dvc(1), $
	plt_rgn_dvc(3), $
	plt_rgn_dvc(4)]
;	
;
clr_mk,num_cbar_colors,color_order,pll_idx,1
;
tvlct,r_curr,g_curr,b_curr
;
clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=num_cbar_colors, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz
;
title_sz=1.6
;
contour, $
	data, $
	abc, $
	ord, $
	tit=title_title, $
	xtit=x_axz_title, $
	ytit=y_axz_title, $
	level=cntr_lvls, $                 
	c_color=cntr_fll_idx, $		
	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
;	font=0, $				
;	xrange=[abc_min,abc_max], $
;	yrange=[ord_min,ord_max], $
	xstyle=1, $
	ystyle=1, $
	charsize=title_sz, $
	position=plt_rgn_nrm, $
	ycharsize=1.1, $
	xcharsize=1.1, $
	/closed, $				
	/fill, $				
	/noerase
;
contour, $
	data, $
	abc, $
	ord, $
	level=cntr_lvls, $                 
	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
	/closed, $				
	/overplot
;
end; end IWC_color()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure IWC color commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure IWC B&W.
; Contour plot the IWC vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_size] (in C)
; is accessed as         foo(1..num_size,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro IWC_bw
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_size",num_size
ncdf_varget1,cdfid,"num_step",num_step
ncdf_varget1,cdfid,"num_frame",num_frame
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"time_snapshot",time_snapshot
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"IWC_snapshot",IWC_snapshot
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
if n_elements(time_snapshot) le 1 then begin
	print,"Requested data not in NetCDF file, returning."
	return
endif
;
data=fltarr(num_frame+1,num_layer)
for layer=1,num_layer do begin
data(*,layer-1)=IWC_snapshot(layer,*)
endfor
data=data*1.0e6
abc=time_snapshot(*)/3600.
ord=altitude(1:num_layer)/1000.
title_title='!5Ice Water Content !8IWC(t,z)!5 (mg-m!E-3!N)'
x_axz_title='Time !8t!5 (hrs.)'
y_axz_title='Altitude !8z!5 (km)'
abc_min=0.
abc_max=650.
ord_min=14.5
ord_max=17.5
;
cntr_lvl_nbr=22
grd_ntv=5.
cntr_lvl_min=0.
first_level=0.5
;
data_min=min(data)
data_max=max(data)
max_level=(data_max+grd_ntv) - abs(data_max mod grd_ntv)
cntr_ntv=(max_level-cntr_lvl_min)/(cntr_lvl_nbr-2)
cntr_lvls=fltarr(cntr_lvl_nbr)
cntr_lvls(0)=cntr_lvl_min
if first_level lt data_min then first_level=cntr_ntv
cntr_lvls(1)=first_level
for i=2,cntr_lvl_nbr-1 do begin
	cntr_lvls(i)=cntr_lvls(i-1)+cntr_ntv
endfor
;
contour, $
	data, $
	abc, $
	ord, $
	tit=title_title, $
	xtit=x_axz_title, $
	ytit=y_axz_title, $
	level=cntr_lvls, $                 
;	xrange=[abc_min,abc_max], $
;	yrange=[ord_min,ord_max], $
	xstyle=1, $
	ystyle=1, $
	/closed					
;
end; end IWC_bw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure IWC B&W commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Heating Rates
; Plot the vertical heating rate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro heating_rates
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename_1=primary_fl
cdfid_1=ncdf_open(filename_1)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_layer",num_layer_1
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"altitude",altitude_1
ncdf_varget,cdfid_1,"heating_rate_LW",LW_heating_1
ncdf_varget,cdfid_1,"heating_rate_SW",SW_heating_1
ncdf_varget,cdfid_1,"heating_rate_net",net_heating_1
;
; say bye-bye
ncdf_close,cdfid_1
; End of NetCDF commands
;
; Allow the y-axis to be within .5 degrees of the curves
data_min=3600.*min([min(LW_heating_1),min(SW_heating_1),min(net_heating_1)])
data_max=3600.*max([max(LW_heating_1),max(SW_heating_1),max(net_heating_1)])
data_min=(data_min-.5) + abs(data_min mod .5)
data_max=(data_max+.5) - abs(data_max mod .5)
;
plot, $
	net_heating_1(1:num_layer_1)*3600.0, $
	altitude_1(1:num_layer_1)/1000.0, $
	tit='!5 Heating Rates', $
	xtit='Heating Rate (K-hr!E-1!N)', $
	ytit='Altitude !8z!5 (km)', $
	xrange=[data_min,data_max], $
;	yrange=[8.,18.], $
	xmargin=[5.5,1.2], $	;[10,3] is [left,right] default
	ymargin=[3,1.5], $	;[4,2] is [bottom,top] default
	xstyle=1, $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	SW_heating_1(1:num_layer_1)*3600.0, $
	altitude_1(1:num_layer_1)/1000.0, $
	thick=3.0, $
	linestyle=1
;
oplot,	$
	LW_heating_1(1:num_layer_1)*3600.0, $
	altitude_1(1:num_layer_1)/1000.0, $
	thick=3.0, $
	linestyle=2
;
ln_lgn_x1=0.60
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.5
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5SW',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5LW',size=txt_lgn_sz,/NORMAL
;
end; end heating_rates()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Heating Rates commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Radiative Forcing
; Plot the time evolution of radiative forcing for the three main cases.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rad_forcing
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename_1=primary_fl
cdfid_1=ncdf_open(filename_1)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"LW_cloud_forcing",LW_forcing_1
ncdf_varget,cdfid_1,"SW_cloud_forcing",SW_forcing_1
ncdf_varget,cdfid_1,"time_array",time_array_1
;
; say bye-bye
ncdf_close,cdfid_1
; End of NetCDF commands
;
if n_elements(time_array_1) le 1 then begin
	print,"Requested data not in NetCDF file, returning."
	return
endif
;
; Convert to Diurnally Averaged forcings...
; The 1/pi factor is the diurnal average, and the 2/sqrt(3) factor converts
; the 30 degree zenith angle simulations into pseudo local-noon quantities.
SW_forcing_1=SW_forcing_1*2./(sqrt(3.)*!pi)
;
; ...and compute net forcings
net_forcing_1=SW_forcing_1+LW_forcing_1
;
; Allow the top and bottom boundaries to be within 50 W/m^2 of the curves
data_min=min([min(SW_forcing_1),min(LW_forcing_1),min(net_forcing_1)])
data_max=max([max(LW_forcing_1),max(SW_forcing_1),max(net_forcing_1)])
data_min=(data_min-50.) + abs(data_min mod 50.)
data_max=(data_max+50.) - abs(data_max mod 50.)
;
plot, $
	time_array_1(0:num_step_1-1)/3600.0, $
	net_forcing_1(0:num_step_1-1), $
	tit='!5Diurnally Averaged Cloud Forcing', $
	xtit='!5Elapsed Time !8t!5 (min.)', $
	ytit='!5Forcing (W-m!E-2!N)', $
	xstyle=1, $
	yrange=[data_min,data_max], $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
; Now overplot the SW forcings
oplot,	$
	time_array_1(0:num_step_1-1)/3600.0, $
	SW_forcing_1(0:num_step_1-1), $
	thick=3.0, $
	linestyle=1
;
; Now overplot the LW forcings
oplot,	$
	time_array_1(0:num_step_1-1)/3600.0, $
	LW_forcing_1(0:num_step_1-1), $
	thick=3.0, $
	linestyle=2
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.5
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5SWCF',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5LWCF',size=txt_lgn_sz,/NORMAL
;
end; end rad_forcing()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Radiative Forcing commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Heating Distribution 
; Contour plot the final integrated heating rate distribution vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_size] (in C)
; is accessed as         foo(1..num_size,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro heat_dst_bw
@cld_cmn.com
erase
;
neg_dst_cntr_lvls=[ $
  -128,-64,-32,-16,-8, $
  -4,-2,-1,-.5,-.25, $
  0.0, $
  .25,.5,1,2,4, $
  8,16,32,64,128, $
  256,512,1024,2048,4096, $
  8192]
;
heat_dst_cntr_lvls=[ $
  -128,-64,-32,-16,-1, $
  -.5,-.4,-.3,-.2,-.1, $
  0.0, $
  .1,0.2,.3,.4,.5, $
  1,2,5,10,15, $
  20,25,30,35,40, $
  8192]
;
heat_which_lbl=[ $
  1,1,1,1,1, $	
  1,1,1,1,1, $	
  1, $
  1,1,1,1,1, $	
  1,1,1,1,0, $	
  1,0,1,0,1, $	
  1]
;
heat_ntt=[ $
  '-128','-64','-32','-16','-1', $
  '-.5','-.4','-.3','-.2','-.1', $
  '0', $
  '.1','.2','.3','.4','.5', $
  '1','2','5','10','15', $
  '20','25','30','35','40', $
  '8192']
;
neg_dst_ln_sty=[ $
  1,1,1,1,1, $
  1,1,1,1,1, $
  0, $
  0,0,0,0,0, $
  0,0,0,0,0, $
  0,0,0,0,0, $
  0]
;
neg_dst_thick=[ $
  2,2,2,2,2, $
  2,2,2,2,2, $
  2, $
  1,1,1,1,1, $
  1,1,1,1,1, $
  1,1,1,1,1, $
  1]
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_size",num_size
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"effective_radius",effective_radius
ncdf_varget,cdfid,"env_density",env_density
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
ncdf_varget,cdfid,"crystal_heat_SW",SW_int_heat
ncdf_varget,cdfid,"crystal_heat_LW",LW_int_heat
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
spec_heat_const_pres=1005.	;J/kg/K
seconds_per_hour=3600.		
distribution=concentration
net_int_heat=concentration
for layer=1,num_layer do begin
	distribution(1:num_size,layer)=concentration(1:num_size,layer)/ $
		delta_length(1:num_size) ; #/m^3/m
	net_int_heat(1:num_size,layer)=distribution(1:num_size,layer)* $
		(SW_int_heat(1:num_size,layer)+LW_int_heat(1:num_size,layer))
	net_int_heat(1:num_size,layer)=net_int_heat(1:num_size,layer)/ $
		(env_density(layer)*spec_heat_const_pres)	; K/s/m
endfor
net_int_heat=net_int_heat*seconds_per_hour/1.0e6	; K/hour/micron
net_int_heat=net_int_heat*1.0e3			; mK/hour/micron
;print,max(net_int_heat)
;print,min(net_int_heat)
;net_int_heat=net_int_heat*1.0e6/1.0e6	; micro-Watts/m^3/micron
;
contour,net_int_heat(1:num_size,1:num_layer), $
	alog10(crystal_length(1:num_size)*1.0e6), $
;	crystal_length(1:num_size)*1.0e6, $
	altitude(1:num_layer)/1000.0, $	
	tit='!5 Heating Rate Distribution (!5mK-hr!E-1!N-!7l!5m!E-1!N) in High Thin Cirrus', $
;	xtit='Crystal Length !8L!5 (!7l!5m)', $
	xtit='Log!I10!N(Crystal Length !8L!5 (!7l!5m))', $
	ytit='Altitude !8z!5 (km)', $
	level=heat_dst_cntr_lvls, $       
	c_labels=heat_which_lbl, $		
	c_linestyle=neg_dst_ln_sty, $	
	c_thick=neg_dst_thick, $		
	c_annotation=heat_ntt, $	
;	font=0, $				
	/follow, $				
;	xrange=[0.0,1000.0], $
	xrange=[.5,3.], $
	xstyle=1, $
	ystyle=1, $
;	/min_curve_surf, $
	/noerase
;
effective_length=effective_radius
for layer=1,num_layer do begin
;	effective_length(layer)=equiv_rad_to_bullet_length(effective_radius(layer))
	effective_length(layer)=equiv_rad_to_hex_clm_length(effective_radius(layer))
	if altitude(layer)/1000. gt 17.2 then effective_length(layer)=1.0e-6
endfor
;cloud_eff_length=equiv_rad_to_bullet_length(cloud_eff_radius(num_step))
cloud_eff_length=equiv_rad_to_hex_clm_length(cloud_eff_radius(num_step))
;
oplot,	$
;	effective_length(1:num_layer)*1.0e6, $
	alog10(effective_length(1:num_layer)*1.0e6), $
	altitude(1:num_layer)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
;	[cloud_eff_length,cloud_eff_length]*1.0e6, $
	alog10([cloud_eff_length,cloud_eff_length]*1.0e6), $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
end; end heat_dst_bw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Heating Distribution commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Temperature
; Plot the temperature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro temperature
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"env_temperature",env_temperature
ncdf_varget,cdfid,"env_pressure",env_pressure
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
plot, $
	env_temperature(1:num_layer)-273.15, $
	altitude(1:num_layer)/1000.0, $
	tit='!5 Env_Temperature', $
	xtit='!5Env_Temperature (C)', $
	ytit='!5Altitude !8z!5 (km)', $
;	xrange=[0.0,1.0], $
;	yrange=[8.,18.], $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	env_pressure(1:num_layer)/100.0, $
	altitude(1:num_layer)/1000.0, $
	thick=3.0, $
	linestyle=1
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5env_temperature',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5env_pressure',size=txt_lgn_sz,/NORMAL
;
end; end temperature()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure temperature commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CCN_Activated
; Plot the CCN_activated
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro CCN_activated
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"CCN_activated",CCN_activated
ncdf_varget,cdfid,"saturation_liquid",saturation_liquid
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
plot, $
	CCN_activated(1:num_layer), $
	altitude(1:num_layer)/1000.0, $
	tit='!5 CCN_Activated', $
	xtit='!5CCN_Activated ()', $
	ytit='!5Altitude !8z!5 (km)', $
;	xrange=[0.0,1.0], $
;	yrange=[8.,18.], $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	saturation_liquid(1:num_layer), $
	altitude(1:num_layer)/1000.0, $
	thick=3.0, $
	linestyle=1
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5CCN_activated',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5saturation_liquid',size=txt_lgn_sz,/NORMAL
;
end; end CCN_activated()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CCN_Activated commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Albedo Emissivity.
; Plot albedo and emissivity.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro alb_emis
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename_1=primary_fl
cdfid_1=ncdf_open(filename_1)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"albedo_of_time",albedo_1
ncdf_varget,cdfid_1,"emissivity_of_time",emissivity_1
ncdf_varget,cdfid_1,"time_array",time_array
;
; say bye-bye
ncdf_close,cdfid_1
; End of NetCDF commands
;
if n_elements(time_array) le 1 then begin
	print,"Requested data not in NetCDF file, returning."
	return
endif
;
plot, $
	emissivity_1(0:num_step_1-1), $
	albedo_1(0:num_step_1-1), $
	tit='!5 Emissivity vs. Albedo', $
	xtit='Emissivity !7e!5', $
	ytit='Albedo !8A!5', $
;	font=0, $				
;	xrange=[.85,1.01], $
;	yrange=[0.1,0.6], $
	xstyle=1, $
	ystyle=0, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
ln_lgn_x1=.3
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Cloud',size=txt_lgn_sz,/NORMAL
;
end; end alb_emis()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Albedo Emissivity commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Albedo Emissivity of time.
; Plot the albedo and emissivity as a function of time.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro alb_emis_of_time
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename_1=primary_fl
cdfid_1=ncdf_open(filename_1)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"albedo_of_time",albedo_1
ncdf_varget,cdfid_1,"emissivity_of_time",emissivity_1
ncdf_varget,cdfid_1,"time_array",time_array
;
; say bye-bye
ncdf_close,cdfid_1
; End of NetCDF commands
;
if n_elements(time_array) le 1 then begin
	print,"Requested data not in NetCDF file, returning."
	return
endif
;
plt_rgn_nrm=[ $
	.13, $ ; x_min
	.11, $ ; y_min
	.87, $ ; x_max
	.90] ; y_max
;
plot, $
	time_array(0:num_step_1-1)/3600.0, $
	albedo_1(0:num_step_1-1), $
	tit='!5 Albedo and Emissivity vs. Time', $
	xtit='Time !8t!5 (hrs.)', $
	ytit='Albedo !8A!5', $
	xstyle=1, $
	ystyle=8, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	position=plt_rgn_nrm, $
	linestyle=0
;
plot, $
	time_array(0:num_step_1-1)/3600.0, $
	emissivity_1(0:num_step_1-1), $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=1, $
	position=plt_rgn_nrm, $
	/noerase
;
axis, $
	yaxis=1, $
	ystyle=0, $
	/ynozero, $
	ytitle='!5Emissivity !7e!5', $
	charsize=2.
;
ln_lgn_x1=.3
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Albedo',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5Emissivity',size=txt_lgn_sz,/NORMAL
;
end; end alb_emis_of_time()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Albedo Emissivity of time commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure Foo
; Plot the foo
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro foo
@cld_cmn.com
erase
;
; Read in binary data from the NetCDF file 
filename=primary_fl
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"foo",foo
ncdf_varget,cdfid,"foo",foo
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
plot, $
	foo(1:num_layer), $
	altitude(1:num_layer)/1000.0, $
	tit='!5 Foo', $
	xtit='!5Foo ()', $
	ytit='!5Altitude !8z!5 (km)', $
;	xrange=[0.0,1.0], $
;	yrange=[8.,18.], $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	foo(1:num_layer), $
	altitude(1:num_layer)/1000.0, $
	thick=3.0, $
	linestyle=1
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5foo',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5foo',size=txt_lgn_sz,/NORMAL
;
end; end foo()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Foo commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




