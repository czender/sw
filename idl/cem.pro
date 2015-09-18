; $Id$

; Purpose: offline analysis routines for looking at CEM data. 
; These routines are heavily reliant on netcdf.

; $Log: not supported by cvs2svn $
; Revision 1.11  2000-10-10 22:52:14  zender
; *** empty log message ***
;
; Revision 1.10  2000/01/15 02:07:48  zender
; *** empty log message ***
;
; Revision 1.7  2000/01/01 01:55:47  zender
; *** empty log message ***
;
; Revision 1.6  1999/12/31 20:12:44  zender
; *** empty log message ***
;
; Revision 1.5  1999/12/31 02:09:35  zender
; *** empty log message ***
;
; Revision 1.4  1999/12/31 00:18:14  zender
; *** empty log message ***
;
; Revision 1.3  1999/12/30 19:34:23  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/04 23:37:11  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:45  zender
; Imported sources
;
; Revision 1.4  1995/12/28  00:12:21  zender
; added _bch() routines to all procedures with a graphic output
; for the spcp_01 paper.
;
; Revision 1.3  1995/11/10  23:56:10  zender
; starting to write thesis
;
; Revision 1.2  1994/10/21  01:50:33  zender
; added lag correlations and lots of functionality to cem_time.
; and there's an animation routine as well.
;
; Revision 1.1  1994/09/21  02:13:34  zender
; Initial revision

@ibp_clr.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin CEM Batch
; This batch routine calls cem_cmfmc with all the interesting
; pressure levels.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_bch
;
approx_p=[850.0,700.0,500.0,400.0,300.0,250.0,200.0]
lev=[6,10,16,19,23,26,29]
lev_nbr=n_elements(lev)
;
open_ps,fl_nm='/data/zender/ps/cem_mc_hst_500mb.eps',x_sz=5.1,y_sz=3.0,/eps
cem_cmfmc,start_time_min=4320,end_time_min=34500,incr_time_min=60,level=16,ttl_sng='CEM',rng_y=[0,30]
close_ps,fl_nm='/data/zender/ps/cem_mc_hst_500mb.eps'
;
open_ps,fl_nm='/data/zender/ps/cem_mc_hst_200mb.eps',x_sz=5.1,y_sz=3.0,/eps
cem_cmfmc,start_time_min=4320,end_time_min=34500,incr_time_min=60,level=29,ttl_sng='CEM',rng_y=[0,45]
close_ps,fl_nm='/data/zender/ps/cem_mc_hst_200mb.eps'
;
;for level=0,lev_nbr-1 do begin
;
;	cem_cmfmc,start_time_min=4320,end_time_min=34500,incr_time_min=60,level=lev(level)
;;	cem_cmfmc,start_time_min=0,end_time_min=1430,incr_time_min=10,level=lev(level)
;;
;endfor
;
end; end cem_bch()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End CEM Batch
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM Time
; Analyze the time series of (semi) arbitrary CEM fields
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_time_bch
;
open_ps,fl_nm='/data/zender/ps/cem_mc_diwpdt_a.eps',/one_clm,/eps
cem_time,start_time_min=0,end_time_min=1430,incr_time_min=10,coplot=1,scatplot=0,lagplot=0,slice_lbl=0
close_ps,fl_nm='/data/zender/ps/cem_mc_diwpdt_a.eps'
;
open_ps,fl_nm='/data/zender/ps/cem_mc_diwpdt_b.eps',/one_clm,/eps
cem_time,start_time_min=0,end_time_min=150,incr_time_min=10,coplot=1,scatplot=0,lagplot=0,slice_lbl=0
close_ps,fl_nm='/data/zender/ps/cem_mc_diwpdt_b.eps'
;
open_ps,fl_nm='/data/zender/ps/cem_mc_diwpdt_c.eps',/one_clm,/eps
cem_time,start_time_min=0,end_time_min=150,incr_time_min=10,coplot=0,scatplot=1,lagplot=0,slice_lbl=0
close_ps,fl_nm='/data/zender/ps/cem_mc_diwpdt_c.eps'
;
end; end cem_time_bch()
pro cem_time, $
	IWC=IWC, $
	alpha=alpha, $
	beta=beta, $
	chr_sz=chr_sz, $
	clean=clean, $
	coplot=coplot, $
	delta=delta, $
	end_time_min=end_time_min, $
	gamma=gamma, $
	iciwp=iciwp, $
	incr_time_min=incr_time_min, $
	info=info, $
	lag_time_min=lag_time_min, $
	lagplot=lagplot, $
	level=level, $
	name_0=name_0, $
	name_1=name_1, $
	fld_nbr=fld_nbr, $
	same=same, $
	scatplot=scatplot, $
	slice_lbl=slice_lbl, $
	solo=solo, $
	start_time_min=start_time_min
;
@cem_sct.pro
;
; cem_time,start_time_min=0,end_time_min=1430,incr_time_min=10
; cem_time,start_time_min=0,end_time_min=150,incr_time_min=10
; cem_time,start_time_min=4320,end_time_min=34500,incr_time_min=60
; cem_time,name_0='conv_avg_mass_flux',name_1='dom_avg_CWP_partial'
; cem_time,name_0='dom_avg_IWP_partial',start_time_min=240,end_time_min=520,lag_time_min=-60
; cem_time,name_0='dom_avg_RH_ice',IWC=1,start_time_min=0,end_time_min=1430,level=[16,16,16,16,16,0]
; cem_time,name_0='dom_avg_q_vapor',IWC=0,start_time_min=600,end_time_min=1000,level=[23,16,16,16,16,0]
; cem_time,fld_nbr=5,level=[16,16,16,16,16,0]
; cem_time,fld_nbr=2,level=[16,16],name_0='cloud_frac',name_1='conv_frac'
;
if n_elements(start_time_min) eq 0 then start_time_min=0
;if n_elements(end_time_min) eq 0 then end_time_min=1430
if n_elements(end_time_min) eq 0 then if start_time_min lt 3000 then end_time_min=1430 else end_time_min=34500
;if n_elements(incr_time_min) eq 0 then incr_time_min=10
if n_elements(incr_time_min) eq 0 then begin
	if start_time_min lt 3000 then incr_time_min=10
	if start_time_min gt 3000 then incr_time_min=60
endif
if n_elements(lag_time_min) eq 0 then lag_time_min=0
if n_elements(level) eq 0 then level=[16,16,16,16,16,0] ; layer 16 is 500 mb
if n_elements(fld_nbr) eq 0 then fld_nbr=5
if n_elements(clean) eq 0 then clean=1
if n_elements(IWC) eq 0 then IWC=0
if n_elements(solo) eq 0 then solo=0
if n_elements(iciwp) eq 0 then iciwp=0
if n_elements(slice_lbl) eq 0 then slice_lbl=1
if n_elements(same) eq 0 then same=0
if n_elements(chr_sz) eq 0 then chr_sz=1.8
if n_elements(coplot) eq 0 then coplot=1
if n_elements(lagplot) eq 0 then lagplot=1
if n_elements(scatplot) eq 0 then scatplot=1
if n_elements(info) eq 0 then info=1
if n_elements(alpha) eq 0 then alpha=1. ; [g/kg]
if n_elements(beta) eq 0 then beta=1. ; [something weird]
if n_elements(gamma) eq 0 then gamma=1. ; [s^-1]
if n_elements(delta) eq 0 then delta=1. ; [something weird]
;
if n_elements(name_0) eq 0 then name_0='conv_avg_mass_flux'
if n_elements(name_1) eq 0 then begin
if IWC then begin
	name_1='dom_avg_IWC'
endif else begin
	name_1='dom_avg_IWP_partial'
endelse
endif
if n_elements(name_2) eq 0 then name_2='cloud_frac'
if n_elements(name_3) eq 0 then name_3='dom_avg_q_vapor'
if n_elements(name_4) eq 0 then name_4='dom_avg_q_sat_ice'
if n_elements(name_5) eq 0 then name_5='dom_avg_precip'
;
DRW_FIELDS_SOLO=solo
SAME_AXZ=same
ICIWP=iciwp
lev_nbr=52
time_nbr=(end_time_min-start_time_min)/incr_time_min+1
;
name=strarr(fld_nbr)
fld=intarr(fld_nbr)
data=fltarr(fld_nbr,time_nbr)
time=fltarr(time_nbr)
;
name(0)=name_0
name(1)=name_1
if fld_nbr gt 2 then name(2)=name_2
if fld_nbr gt 3 then name(3)=name_3
if fld_nbr gt 4 then name(4)=name_4
if fld_nbr gt 5 then name(5)=name_5
print,'name(0) = ',name(0)
print,'name(1) = ',name(1)
;
abb_sng=strarr(fld_nbr)
axz_ttl_sng=strarr(fld_nbr)
bin_sz=fltarr(fld_nbr)
data_max=fltarr(fld_nbr)
data_min=fltarr(fld_nbr)
data_nbr=fltarr(fld_nbr)
scale=fltarr(fld_nbr)
slice_sng=strarr(fld_nbr)
sym_sng=strarr(fld_nbr)
unit_sng=strarr(fld_nbr)
;
; Loop over the requested time period
;
end_time_min=long(end_time_min)
start_time_min=long(start_time_min)
incr_time_min=long(incr_time_min)
time_min=end_time_min
;
time_idx=0
for time_min=start_time_min,end_time_min,incr_time_min do begin
day=((time_min - (time_min mod 1440))/1440)+1;
minute_of_day=time_min-1440*(day-1);
hour_of_day=(minute_of_day-(minute_of_day mod 60))/60;
;
; Construct the file name
;
day_sng=string(format='(I2.2)',day)
minute_of_day_sng=string(format='(I4.4)',minute_of_day)
fl_nm='cem.'+day_sng+'.'+minute_of_day_sng+'.nc'
;
; Read in binary data from the NetCDF file 
;
fl_path='/data/zender/_aux0_/cem'
fl_nm=fl_path+'/'+fl_nm
nc_id=ncdf_open(fl_nm)
;
for fld_idx=0,fld_nbr-1 do begin
;
; get the one-dimensional arrays
;
ncdf_varget,nc_id,name(fld_idx),tmp_data
;
if n_elements(tmp_data) gt 1 then data(fld_idx,time_idx)=tmp_data(level(fld_idx)) else data(fld_idx,time_idx)=tmp_data
;
endfor; end loop over fld_idx
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
time(time_idx)=time_min
print,'time = ',auto_sng(time_min,0),', Day = ',auto_sng(day,0),' minute = ',auto_sng(minute_of_day,0),', data(0) = ', auto_sng(data(0,time_idx),2),', data(1) = ', auto_sng(data(1,time_idx),2) ; ,', data(2) = ', auto_sng(data(2,time_idx),2)
;
time_idx=time_idx+1
;
endfor ; end loop over times
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Lag loop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; If a lag time is specified then we need to shift the data for
; the second field to the right by the specified amount. So if the
; specified lag is, say, J elements, then right shift the second field 
; by J elements so its last J elements wrap around to become its first
; five elements. Overwrite the J wrapped-around elements of the second
; field with J new values from the J times preceding the start time of
; the first field. This is known as a positive lag of J. For positive
; shifts the first field is lagging (behind in space/time) the second
; field. For negative shifts the first field is leading (ahead in 
; space/time) the second field. 
;
if lag_time_min ne 0 then begin
;
; First, shift the data. Data that is wrapped around will be 
; overwritten by the new data during the next read segment.
;
lag_fld_idx=1 	; i.e., shifting the IWP shifts the dIWP/dt
lag_shift=fix(lag_time_min/incr_time_min) ; lag_shift is J
data(lag_fld_idx,*)=shift(data(lag_fld_idx,*),lag_shift)
;
; If lag_shift is negative, i.e., negative lag, then start writing
; J elements from the end of the array. Start writing at element 0
; for positive lags.
;
if lag_shift lt 0 then lag_idx=time_nbr+lag_shift else lag_idx=0
print,'lag_shift = ',lag_shift
print,'lag_idx = ',lag_idx
;
time_idx=0
for time_min=incr_time_min,abs(lag_time_min),incr_time_min do begin
if lag_shift lt 0 then real_time_min=end_time_min+time_min else real_time_min=start_time_min-lag_time_min-incr_time_min+time_min
day=((real_time_min - (real_time_min mod 1440))/1440)+1;
minute_of_day=real_time_min-1440*(day-1);
hour_of_day=(minute_of_day-(minute_of_day mod 60))/60;
;
; Construct the file name
;
day_sng=string(format='(I2.2)',day)
minute_of_day_sng=string(format='(I4.4)',minute_of_day)
fl_nm='cem.'+day_sng+'.'+minute_of_day_sng+'.nc'
;
; Read in binary data from the NetCDF file 
;
fl_path='/data/zender/_aux0_/cem'
fl_nm=fl_path+'/'+fl_nm
nc_id=ncdf_open(fl_nm)
;
ncdf_varget,nc_id,name(lag_fld_idx),tmp_data
;
if n_elements(tmp_data) gt 1 then data(lag_fld_idx,time_idx+lag_idx)=tmp_data(level(lag_fld_idx)) else data(lag_fld_idx,time_idx+lag_idx)=tmp_data
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
print,'real_time = ',auto_sng(real_time_min,0),', data(0) = ', auto_sng(data(0,time_idx+lag_idx),2),', data(1) = ', auto_sng(data(1,time_idx+lag_idx),2)
;
time_idx=time_idx+1
;
endfor ; end loop over times
;
endif ; end if lag is specified
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Lag loop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; After the time loop, gather some useful information that
; does not change (too much!) with each passing file
;
nc_id=ncdf_open(fl_nm)
ncdf_varget,nc_id,'dom_avg_p',dom_avg_p
ncdf_close,nc_id
;
; If in-cloud condensate is to be used as the prognostic variable,
; then post-process the domain averaged condensate to obtain it here.
;
if ICIWP then begin
	no_cloud_idxs=where(data(2,*) lt 1.0e-2)
	cloud_idxs=where(data(2,*) ge 1.0e-2)
	data(1,no_cloud_idxs)=0.
	data(1,cloud_idxs)=data(1,cloud_idxs)/data(2,cloud_idxs)
endif; end if ICIWP
;
; The fields from disk are ready to be overplotted and correlated. 
; Here is our chance to compute some derived fields while
; all input data is stil in SI unit.
; Find the rate of production of IWP, and various possible
; source and sink terms to it, such as
; cloud fraction, convective mass flux, and saturation.
;
dt=incr_time_min*60.
two_dt=2.*dt
dIWdt=(shift(data(1,*),-1)-shift(data(1,*),1))/two_dt
dIWdt(0)=(data(1,1)-data(1,0))/dt
dIWdt(time_nbr-1)=(data(1,time_nbr-1)-data(1,time_nbr-2))/dt
;
;dqdt=(shift(data(0,*),-1)-shift(data(0,*),1))/two_dt
;dqdt(0)=(data(0,1)-data(0,0))/dt
;dqdt(time_nbr-1)=(data(0,time_nbr-1)-data(0,time_nbr-2))/dt
;
;conv=alpha*data(0,*)
;conv=2.37e-3*data(0,*)
;cond=max([beta*(data(3,*)/data(4,*)-.65),0.0])
;cond=beta*data(3,*)/data(4,*)
;cond=1.310e-3*data(3,*)/data(4,*)
;precip=gamma*data(1,*)
;precip=2.18e-4*data(1,*)
;evap=delta*(data(4,*)-data(3,*))*sqrt(data(5,*))
;conv_precip=conv-precip
;conv_cond_precip=conv+cond-precip
;
;name(0)='dqdt'
;data(0,*)=dqdt
;level(0)=level(0)	; set the slice string to the q layer
;
;name(0)='conv'
;data(0,*)=conv
;level(0)=level(0)	; set the slice string to the M_c layer
;
;name(0)='cond'
;data(0,*)=cond
;level(0)=level(3)	; set the slice string to the q, q* layer
;
;name(0)='precip'
;data(0,*)=precip
;level(0)=level(2)	; set the slice string to the cloud_frac layer
;
;name(0)='evap'
;data(0,*)=evap
;level(0)=level(3)	; set the slice string to the q, q* layer
;
;name(0)='conv_precip'
;data(0,*)=conv_precip
;level(0)=level(0)	; set the slice string to the M_c layer
;
;name(0)='conv_cond_precip'
;data(0,*)=conv_cond_precip
;level(0)=level(0)	; set the slice string to the M_c layer
;
;name(0)=name(1)
;data(0,*)=data(1,*)
;level(0)=level(1)	; set the slice string to the IW layer
;
;name(0)=name(2)
;data(0,*)=data(2,*)
;level(0)=level(2)	; set the slice string to the cloud_frac layer
;
if IWC then begin
	name(1)='dIWCdt'
endif else begin
	name(1)='dIWPdt'
endelse
data(1,*)=dIWdt
level(1)=level(1)	; set the slice string to the IW layer
;
; Assign the field-specific labeling info. First find the index
; of each field into the field structure.
;
for fld_idx=0,fld_nbr-1 do begin
located_fld=False
idx=0
while ((not located_fld) and (idx lt fld_nbr_in_database)) do begin
	if (name(fld_idx) eq fld_lst(idx).name) then begin
		located_fld=True
		fld(fld_idx)=idx
	endif	
	idx=idx+1
endwhile; end loop over idx
if located_fld eq False then print,'Error in cem_time: fld '+name(fld_idx)+' not found in database.'
;; Draw the profiles of each fld loaded before comparing the
; two fields.
;
sym_sng(fld_idx)=fld_lst(fld(fld_idx)).symbol
slice_sng(fld_idx)=auto_sng(round(dom_avg_p(level(fld_idx))/100),0)+' mb'
if not slice_lbl then slice_sng(fld_idx)=''
abb_sng(fld_idx)=fld_lst(fld(fld_idx)).abbrev
axz_ttl_sng(fld_idx)=slice_sng(fld_idx)+' '+fld_lst(fld(fld_idx)).symbol+' ('+fld_lst(fld(fld_idx)).unit+')'
bin_sz(fld_idx)=fld_lst(fld(fld_idx)).bin_sz
scale(fld_idx)=fld_lst(fld(fld_idx)).scale
;
; this is a debugging measure to make pzns. easier
;
;scale(fld_idx)=1.
;scale(fld_idx)=1200.e3
data(fld_idx,*)=data(fld_idx,*)*scale(fld_idx)
data_max(fld_idx)=max(data(fld_idx,*))
data_min(fld_idx)=min(data(fld_idx,*))
data_nbr(fld_idx)=n_elements(data(fld_idx,*))
unit_sng(fld_idx)=fld_lst(fld(fld_idx)).unit
;
if end_time_min gt 3000 then begin
	abc=time/1440
	x_ttl='Model Time !8t!5 (days)'
endif else begin
	abc=time
	x_ttl='Model Time !8t!5 (min.)'
endelse
;
if DRW_FLD_LST_SOLO then begin
;
generic_hst, $
	abb_sng=abb_sng(fld_idx), $
	bin_sz=bin_sz(fld_idx), $
	chr_sz=chr_sz, $
	clean=clean, $
	data=data(fld_idx,*), $
	data_max=data_max(fld_idx), $
	data_min=data_min(fld_idx), $
	data_nbr=data_nbr(fld_idx), $
	rgn_sng='Marshall Islands', $
	slice_sng=slice_sng(fld_idx), $
	sub_ttl_sng='RAMS Instantaneous', $
	ttl_sng=abb_sng(fld_idx), $
	sym_sng=sym_sng(fld_idx), $
	unit_sng=unit_sng(fld_idx)
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
;
plot, $
	abc, $
	data(fld_idx,*), $
	tit=fld_lst(fld(fld_idx)).abbrev+' !8vs.!5 Time', $
	xtit=x_ttl, $
	ytit=axz_ttl_sng(fld_idx), $
	xstyle=1, $
	ystyle=0, $
	/ynozero, $
	thick=3.0, $
	charsize=chr_sz, $
	position=[0.13,.11,0.87,0.90], $ ; x_min,y_min,x_max,y_max
	linestyle=0
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
;
endif; endif DRW_FLD_LST_SOLO
endfor; end loop over fld_idx
if coplot then begin
;
; Compare the first two fields by overplotting their time series
; and doing a regression on their scatterplot.
;
print,abc
print,data(0,*)
plot, $
	abc, $
	data(0,*), $
	tit='Evolution of '+abb_sng(0)+', '+abb_sng(1), $
	xtit=x_ttl, $
	ytit=axz_ttl_sng(0), $
	xstyle=0, $
	ystyle=8, $
	/ynozero, $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[6.5,7.1], $ ;[10,3] is [left,right] default
	ymargin=[3.5,2], $  ;[4,2] is [bottom,top] default
	linestyle=0
;
if SAME_AXZ then begin
;
oplot, $
	abc, $
	data(1,*), $
	thick=3.0, $
	linestyle=2
;
endif else begin; end if plotting two fields on same axes
;
plot, $
	abc, $
	data(1,*), $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[6.5,7.1], $ ;[10,3] is [left,right] default
	ymargin=[3.5,2], $  ;[4,2] is [bottom,top] default
	linestyle=2, $
	/noerase
;
endelse; end if plotting two fields on different axes
;
axis, $
	yaxis=1, $
	ytitle=axz_ttl_sng(1), $
	ymargin=[3.5,2], $  ;[4,2] is [bottom,top] default
	charsize=chr_sz
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif; endif coplot
;
if scatplot then begin
generic_scat, $
	abb_sng=abb_sng(0:1), $
	bin_sz=bin_sz(0:1), $
	chr_sz=chr_sz, $
	clean=clean, $
	data=data(0:1,*), $
	info=info, $
	data_max=data_max(0:1), $
	data_min=data_min(0:1), $
	data_nbr=data_nbr(0:1), $
	rgn_sng='Marshall Islands', $
	slice_sng=slice_sng(0:1), $
	sub_ttl_sng=['RAMS Instantaneous','RAMS Instantaneous'], $
	sym_sng=sym_sng(0:1), $
	unit_sng=unit_sng(0:1)
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
endif; endif scatplot
;
if lagplot then begin
generic_corr, $
	abb_sng=abb_sng(0:1), $
	data=data(0:1,*), $
	delta_time=incr_time_min, $
	data_max=data_max(0:1), $
;	max_mag_lag=max_mag_lag, $
	data_min=data_min(0:1), $
	data_nbr=data_nbr(0:1), $
	rgn_sng='Marshall Islands', $
	slice_sng=slice_sng(0:1), $
	sub_ttl_sng=['RAMS Instantaneous','RAMS Instantaneous'], $
	sym_sng=sym_sng(0:1), $
	unit_sng=unit_sng(0:1)
;
endif; endif lagplot
;
; Summarize the relationship between the two data series by printing
; out their respective time averages.
;
print,'Average of data(0,*)= ',total(data(0,*))/time_nbr
print,'Average of data(1,*)= ',total(data(1,*))/time_nbr
;
end; end cem_time()
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM Time commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM PRECT
; Plot the time series of domain-averaged precipitation.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_prect, $
	start_time_min=start_time_min, $
	end_time_min=end_time_min, $
	incr_time_min=incr_time_min, $
	level=level
;
; cem_prect,start_time_min=0,end_time_min=1430,incr_time_min=10
; cem_prect,start_time_min=4320,end_time_min=34500,incr_time_min=60
;
if n_elements(chr_sz) eq 0 then chr_sz=1.2
if n_elements(start_time_min) eq 0 then start_time_min=0
if n_elements(end_time_min) eq 0 then end_time_min=50
if n_elements(incr_time_min) eq 0 then incr_time_min=10
;if n_elements(foo) eq 0 then foo=1
;
data=0.
time=0.
time_nbr=(end_time_min-start_time_min)/incr_time_min+1
;
; Loop over the requested time period
;
end_time_min=long(end_time_min)
start_time_min=long(start_time_min)
incr_time_min=long(incr_time_min)
time_min=end_time_min
;
for time_min=start_time_min,end_time_min,incr_time_min do begin
day=((time_min - (time_min mod 1440))/1440)+1;
minute_of_day=time_min-1440*(day-1);
hour_of_day=(minute_of_day-(minute_of_day mod 60))/60;
;
; Construct the file name
;
day_sng=string(format='(I2.2)',day)
minute_of_day_sng=string(format='(I4.4)',minute_of_day)
fl_nm='cem.'+day_sng+'.'+minute_of_day_sng+'.nc'
;
; Read in binary data from the NetCDF file 
;
fl_path='/data/zender/_aux0_/cem'
fl_nm=fl_path+'/'+fl_nm
nc_id=ncdf_open(fl_nm)
;
; get the one-dimensional arrays
;
ncdf_varget,nc_id,'dom_avg_precip',dom_avg_precip
;
; say bye-bye
ncdf_close,nc_id
; End of NetCDF commands
;
data=[data,dom_avg_precip]
time=[time,1.+time_min/1440.0]
print,'time = ',auto_sng(time_min,0),', Day = ',auto_sng(day,0),' minute = ',auto_sng(minute_of_day,0),', P_t = ', auto_sng(dom_avg_precip*8.64e7,2),' mm/day.'
;
endfor ; end loop over times
;
data=data*8.64e7 ; convert to mm/day
data=data(1:time_nbr) ; truncate excess space--recall data(0) is garbage
time=time(1:time_nbr) ; truncate excess space--recall time(0) is garbage
;
generic_hst, $
	abb_sng='Domain Avg. Precip', $
	bin_sz=1.0, $
	chr_sz=chr_sz, $
	data=data, $
	data_max=max(data), $
	data_min=min(data), $
	data_nbr=n_elements(data), $
	rgn_sng='Marshall Islands', $
	slice_sng='', $
	sub_ttl_sng='RAMS Instantaneous', $
	sym_sng='!8P!5!IT!N', $
	ttl_sng='Domain Avg. Precip', $
	unit_sng='!5mm day!E-1!N'
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
;
plt_rgn_nrm=[ $
	.13, $ ; x_min
	.11, $ ; y_min
	.87, $ ; x_max
	.90] ; y_max
;
plot, $
	time(0:time_nbr-1), $
	data(0:time_nbr-1), $
	tit='!5Domain-Averaged Precipitation vs. Time', $
	xtit='Model Time !8t!5 (min.)', $
	ytit=''+'!8P!5!IT!N'+' ('+'!5mm day!E-1!N'+')', $
	xstyle=1, $
	ystyle=0, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	position=plt_rgn_nrm, $
	linestyle=0
;
print,' Hit any key to continue...'
junk = get_kbrd(1)
;
generic_pow_spec, $
	data=data, $
	time_nbr=time_nbr, $
	delta_time=incr_time_min*60.0, $
	x_ttl_sng=spec_x_ttl_sng, $
	y_ttl_sng=spec_y_ttl_sng
;
end; end cem_prect()
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM PRECT commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM CMFMC
; Histograph the convective mass flux.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_cmfmc, $
	start_time_min=start_time_min, $
	end_time_min=end_time_min, $
	incr_time_min=incr_time_min, $
	ttl_sng=ttl_sng, $
	rng_y=rng_y, $
	level=level
;
; cem_cmfmc,start_time_min=300,end_time_min=350,incr_time_min=10
; cem_cmfmc,start_time_min=0,end_time_min=1430,incr_time_min=10
; cem_cmfmc,start_time_min=4320,end_time_min=34500,incr_time_min=60
;
if n_elements(chr_sz) eq 0 then chr_sz=1.2
if n_elements(start_time_min) eq 0 then start_time_min=0
if n_elements(end_time_min) eq 0 then end_time_min=50
if n_elements(incr_time_min) eq 0 then incr_time_min=10
if n_elements(ttl_sng) eq 0 then ttl_sng=''
if n_elements(rng_y) le 1 then rng_y=0
if n_elements(level) eq 0 then level=23
;if n_elements(foo) eq 0 then foo=1
;
data=0.
time_nbr=(end_time_min-start_time_min)/incr_time_min+1
;
; Loop over the requested time period
;
end_time_min=long(end_time_min)
start_time_min=long(start_time_min)
incr_time_min=long(incr_time_min)
time_min=end_time_min
;
for time_min=start_time_min,end_time_min,incr_time_min do begin
day=((time_min - (time_min mod 1440))/1440)+1;
minute_of_day=time_min-1440*(day-1);
hour_of_day=(minute_of_day-(minute_of_day mod 60))/60;
;
; Construct the file name
;
day_sng=string(format='(I2.2)',day)
minute_of_day_sng=string(format='(I4.4)',minute_of_day)
fl_nm='cem.'+day_sng+'.'+minute_of_day_sng+'.nc'
;
; Read in binary data from the NetCDF file 
fl_path='/data/zender/_aux0_/cem'
fl_nm=fl_path+'/'+fl_nm
nc_id=ncdf_open(fl_nm)
;
; Get the two-dimensional arrays
;
;offset=[0,level]
;count=[1,1]
ncdf_varget,nc_id,'conv_avg_mass_flux',conv_avg_mass_flux
;ncdf_varget1,nc_id,'conv_avg_mass_flux',conv_avg_mass_flux,offset=offset
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
data=[data,conv_avg_mass_flux(level)] ; level 23 is at 10 km and ~ 300 mb
;data=[data,conv_avg_mass_flux]
print,'time = ',auto_sng(time_min,0),', Day = ',auto_sng(day,0),' minute = ',auto_sng(minute_of_day,0),', M_c = ', auto_sng(conv_avg_mass_flux(level)*1000.0,2),' g/m^2/s.'
;
endfor ; end loop over times
;
; After the time loop, gather some useful information that
; does not change (too much!) with each passing file
;
nc_id=ncdf_open(fl_nm)
ncdf_varget,nc_id,'dom_avg_p',dom_avg_p
ncdf_close,nc_id
;
data=data*1000. ; convert to g/m^2/s
data=data(1:time_nbr) ; truncate excess space--recall data(0) is garbage
;
generic_hst, $
	abb_sng='', $
	bin_sz=2.5, $
	chr_sz=chr_sz, $
	clean=1, $
	data=data, $
	data_max=max(data), $
	data_min=min(data), $
	data_nbr=n_elements(data), $
	rgn_sng='Marshall Islands', $
	slice_sng=auto_sng(round(dom_avg_p(level)/100.),0)+' mb', $
	sub_ttl_sng='RAMS Instantaneous', $
	sym_sng='!8M!5!Ic!N', $
	ttl_sng=ttl_sng, $
	unit_sng='!5g m!e-2!N s!E-1!N', $
	rng_y=rng_y
;
end; end cem_cmfmc()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM CMFMC commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM omega
; Line plot the time-average convective mass flux, the large
; scale mass flux, and the residual mass flux in pressure coordinates.
; This procedure produces a plot to compare to YEC73 fig. 13.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_omega_bch
open_ps,fl_nm='/data/zender/ps/mc_profile.eps',/half,/eps
cem_omega,start_time_min=4320,end_time_min=34500,incr_time_min=60
close_ps,fl_nm='/data/zender/ps/mc_profile.eps'
end; end cem_omega_bch()
pro cem_omega, $
	start_time_min=start_time_min, $
	end_time_min=end_time_min, $
	incr_time_min=incr_time_min, $
	level=level
;
; cem_omega,start_time_min=4320,end_time_min=34500,incr_time_min=60

; cem_omega,start_time_min=4320,end_time_min=5000,incr_time_min=60
;
if n_elements(start_time_min) eq 0 then start_time_min=4320
if n_elements(end_time_min) eq 0 then end_time_min=34500
if n_elements(incr_time_min) eq 0 then incr_time_min=60
;if n_elements(foo) eq 0 then foo=1
;
lev_nbr=52
time_nbr=(end_time_min-start_time_min)/incr_time_min+1
print,'time_nbr = ',time_nbr
print,'lev_nbr = ',lev_nbr
time_dom_avg=fltarr(lev_nbr)*0.
time_conv_avg=fltarr(lev_nbr)*0.
time=0.
;
end_time_min=long(end_time_min)
start_time_min=long(start_time_min)
incr_time_min=long(incr_time_min)
time_min=end_time_min
;
; Loop over the requested time period
;
for time_min=start_time_min,end_time_min,incr_time_min do begin
day=((time_min - (time_min mod 1440))/1440)+1;
minute_of_day=time_min-1440*(day-1);
hour_of_day=(minute_of_day-(minute_of_day mod 60))/60;
;
; Construct the file name
;
day_sng=string(format='(I2.2)',day)
minute_of_day_sng=string(format='(I4.4)',minute_of_day)
fl_in='cem.'+day_sng+'.'+minute_of_day_sng+'.nc'
;
; Read in binary data from the NetCDF file 
;
fl_path='/data/zender/_aux0_/cem'
fl_in=fl_path+'/'+fl_in
nc_id=ncdf_open(fl_in)
;
; get the two-dimensional arrays
;
ncdf_varget,nc_id,'conv_frac',conv_frac
ncdf_varget,nc_id,'conv_avg_omega',conv_avg_omega
ncdf_varget,nc_id,'dom_avg_omega',dom_avg_omega
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
; Accumulate the data for later averaging, print file info. to stdout
;
time_dom_avg=time_dom_avg+dom_avg_omega
time_conv_avg=time_conv_avg+conv_avg_omega
time=[time,1.+time_min/1440.0] ; time is in days
;
print,'time = ',auto_sng(time_min,0),', Day = ',auto_sng(day,0),' minute = ',auto_sng(minute_of_day,0),', fraction 300 mb conv. active grdpts. = ',auto_sng(conv_frac(23)*100.0,0)
;
endfor ; end loop over times
;
; After the time loop, gather some useful information that
; does not change (too much!) with each passing file
;
nc_id=ncdf_open(fl_in)
ncdf_varget,nc_id,'w_thr',w_thr
ncdf_varget,nc_id,'omega_forcing',omega_forcing
ncdf_varget,nc_id,'dom_avg_p',dom_avg_p
ncdf_close,nc_id
if n_elements(dom_avg_p) ne lev_nbr then print,'Something wrong with levels'
;
time=time(1:time_nbr) ; truncate excess space--recall time(0) is garbage
time_dom_avg=time_dom_avg/time_nbr
time_conv_avg=time_conv_avg/time_nbr
;
time_dom_avg=time_dom_avg*3600.; convert to Pa/hour
time_conv_avg=time_conv_avg*3600.; convert to Pa/hour
time_dom_avg=time_dom_avg/100.; convert to mb/hour
time_conv_avg=time_conv_avg/100.; convert to mb/hour
;
abc_1=-time_conv_avg(0:lev_nbr-1)
abc_2=-time_dom_avg(0:lev_nbr-1)
abc_3=abc_2-abc_1
abc_4=-omega_forcing*3600./100.
;
ord=dom_avg_p(0:lev_nbr-1)/100.
;
if 1 then begin
nc_id=ncdf_open('/data/zender/_aux0_/omega0_2/omega_txyavg_0114.nc')
;nc_id=ncdf_open('/data/zender/_aux0_/omega0_2/omega_0114_CMFMC_avg.nc')
ncdf_varget,nc_id,'CMFMC',omega_CMFMC
ncdf_close,nc_id
omega_CMFMC=reform(omega_CMFMC)
;
; Convert the omega mass fluxes from kg/m2/s to mb/hr
;
omega_CMFMC(0)=0
omega_CMFMC=omega_CMFMC*3600.*9.8/100.
omega_CMFMC(0)=1.0e36
print,omega_CMFMC
abc_5=omega_CMFMC
ord_5=[1000.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0]
;
; _6 is eyeballed from YEC73 Fig. 13
;
abc_6=[8.,6.5,5.4,5.4,6.,5.8,4.8,2.7,0.85,0.0]
ord_6=[950.0,850.0,700.0,500.0,400.0,300.0,250.0,200.0,150.0,100.0]
;
endif; end if 1
;
plot, $
	abc_1, $
	ord, $
	tit='!5Mean Convective Mass Flux', $
	xtit='!8M!5(!8p!5) (mb hr!E-1!N)', $
	ytit='!5Pressure !8p!5 (mb)', $
	max_value=1.0e20, $
;	xrange=[0.0,0.06], $
;	xrange=[-10.0,10.0], $
	xrange=[0.0,10.0], $
	yrange=[1000.0,100.0], $
	xstyle=1, $
	ystyle=1, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	xmargin=[7.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3.2,2], $  ;[4,2] is [bottom,top] default
;	ticklen=1.0, $
	linestyle=2
;
;oplot,	$
;	abc_2, $
;	ord, $
;	thick=3.0, $
;	linestyle=1
;;
;oplot,	$
;	abc_3, $
;	ord, $
;	thick=3.0, $
;	linestyle=2
;;
;oplot,	$
;	abc_4, $
;	ord, $
;	thick=3.0, $
;	linestyle=3
;;
if 1 then begin
oplot,	$
	abc_5(1:9), $
	ord_5(1:9), $
	thick=3.0, $
	linestyle=1
oplot,	$
	abc_6, $
	ord_6, $
	thick=3.0, $
	linestyle=0
endif; end if 1
;
ln_lgn_x1=0.65
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=.7
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5YEC73',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5GCM',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5CEM',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(0),'!8M!I!5c!N !7X!I.2!5',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(1),'!8M!I!5t!N (Total)',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(2),'!8M!I!5r!N (Residual)',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(3),'!8M!I!5f!N (Forcing)',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(4),'!5Days '+auto_sng(round(min(time)),0)+' - '+auto_sng(fix(max(time))),size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(5),'!8w!I!5min!N = '+auto_sng(w_thr,1)+' m/s',size=txt_lgn_sz,/NORMAL
;
end; end cem_omega()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM omega commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM vz
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_vz_bch
open_ps,fl_nm='/data/zender/ps/cem_forcing.eps',/one_clm,/eps
cem_vz,fl_nm='/data/zender/_aux0_/cem/cem_forcing.nc',scale=1.0,rng_y=20.0,ttl=''
close_ps,fl_nm='/data/zender/ps/cem_forcing.eps'
end; end cem_vz_bch()
pro cem_vz, $
	fl_nm=fl_nm, $
	fld_nm=fld_nm, $
	rng_y=rng_y, $
	ttl=ttl, $
	scale=scale, $
	level=level
;
if n_elements(scale) eq 0 then scale=1.
if n_elements(ttl) eq 0 then ttl=''
if n_elements(rng_y) eq 0 then rng_y=20.
if n_elements(fl_nm) eq 0 then fl_nm='/data/zender/_aux0_/cem/cem_forcing.nc'
if n_elements(fld_nm) eq 0 then fld_nm='q_condensed'
;
erase
;
; Get the requested field
;
nc_id=ncdf_open(fl_nm)
ncdf_varget,nc_id,'u',u
ncdf_varget,nc_id,'w_forcing',w
ncdf_varget,nc_id,'lev',lev
ncdf_close,nc_id
;
;nc_id=ncdf_open('/data/zender/_aux0_/cem/cem_0424_avg.nc')
;ncdf_varget,nc_id,'u',u
;ncdf_varget,nc_id,'w',w
;ncdf_varget,nc_id,'lev',lev
;ncdf_close,nc_id
;
lev_nbr=n_elements(lev)
;
u=reform(u)
w=reform(w)
;
abc_1=u
abc_2=w*100.
ord=lev/1000.
;
plot, $
	abc_1, $
	ord, $
;	tit='!5Mean Quantity', $
	xtit='!5Imposed Zonal Wind !8u!5 (m s!E-1!N)', $
	ytit='!5Altitude !8z!5 (km)', $
	xstyle=8, $
	ystyle=1, $
	yrange=[0.0,rng_y], $
	thick=3.0, $
	charsize=2.0, $
	xmargin=[5.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3.5,2.5], $  ;[4,2] is [bottom,top] default
	linestyle=0
;
plot, $
	abc_2, $
	ord, $
	xstyle=12, $
	ystyle=12, $
	xrange=[-3.0,3.], $
	yrange=[0.0,rng_y], $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=1, $
;	position=plt_rgn_nrm, $
	xmargin=[5.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3.5,2.5], $  ;[4,2] is [bottom,top] default
	/noerase
;
axis, $
	xaxis=1, $
	xtit='!5Imposed Vertical Updraft !8w!5 (cm s!E-1!N)', $
	charsize=2.
;
ln_lgn_x1=.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.6
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8u!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8w!5',size=txt_lgn_sz,/NORMAL
;
end; end cem_vz()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM vz commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM tz
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_tz, $
	fl_nm=fl_nm, $
	fld_nm=fld_nm, $
	rng_y=rng_y, $
	ttl=ttl, $
	scale=scale, $
	level=level
;
; cem_tz,fl_nm='/data/zender/_aux0_/cem/cem.10.0300.nc',scale=1000,ttl='!5Condensate Mixing Ratio (g kg!E-1!N)',fld_nm='q_condensed'
;
if n_elements(scale) eq 0 then scale=1.
if n_elements(ttl) eq 0 then ttl=''
if n_elements(rng_y) eq 0 then rng_y=0
if n_elements(fl_nm) eq 0 then fl_nm='/data/zender/_aux0_/cem/cem.01.1380.nc'
if n_elements(fld_nm) eq 0 then fld_nm='q_condensed'
;
erase
;
; Get the requested field
;
nc_id=ncdf_open(fl_nm)
ncdf_varget,nc_id,fld_nm,data_1
ncdf_varget,nc_id,'u',u
ncdf_varget,nc_id,'w',w
ncdf_varget,nc_id,'lev',lev
ncdf_close,nc_id
;
lev_nbr=n_elements(lev)
;
; Accumulate the data for later averaging, print file info. to stdout
;
abc_1=data_1*scale
abc_2=data_2*scale
ord=lev/1000.
;
plot, $
	abc_1, $
	ord, $
	tit='!5Mean Quantity', $
	xtit='!5Altitude !8z!5 (km)', $
	ytit='!5Speed (m s!E-1!N)', $
	xstyle=1, $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	xmargin=[7.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3,2], $  ;[4,2] is [bottom,top] default
	linestyle=0
;
oplot,	$
	abc_2, $
	dom_avg_p(0:lev_nbr-1)/100.0, $
	thick=3.0, $
	linestyle=1
;
ln_lgn_x1=0.60
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
xyouts,txt_lgn_x,lgn_y(0),'!8M!I!5c!N (Convective)',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8M!I!5t!N (Total)',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(4),'!5Days '+auto_sng(round(min(time)),0)+' - '+auto_sng(fix(max(time))),size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(5),'!8w!I!5min!N = '+auto_sng(w_thr,1)+' m/s',size=txt_lgn_sz,/NORMAL
;
end; end cem_tz()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM tz commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM IWC
; Line plot the time-Mean IWC, and the time-Mean IWC as a % of
; total IWC in pressure coordinates. The stencil for this procedure
; was CEM omega. This routine will be used to construct a pzn. to 
; redistribute prognosed IWC in the CCM into layers like an anvil.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_IWC_bch
open_ps,fl_nm='/data/zender/ps/cem_fice.eps',/one_clm,/eps
cem_IWC,start_time_min=4320,end_time_min=34500,incr_time_min=60
close_ps,fl_nm='/data/zender/ps/cem_fice.eps'
end; end cem_IWC_bch()
pro cem_IWC, $
	start_time_min=start_time_min, $
	end_time_min=end_time_min, $
	incr_time_min=incr_time_min
;
; cem_IWC,start_time_min=4320,end_time_min=34500,incr_time_min=60
;
if n_elements(start_time_min) eq 0 then start_time_min=4320
if n_elements(end_time_min) eq 0 then end_time_min=34500
if n_elements(incr_time_min) eq 0 then incr_time_min=60
;if n_elements(foo) eq 0 then foo=1
;
lev_nbr=52
time_nbr=(end_time_min-start_time_min)/incr_time_min+1
print,'time_nbr = ',time_nbr
print,'lev_nbr = ',lev_nbr
time_dom_avg_CWC=fltarr(lev_nbr)*0.
time_dom_avg_IWC=fltarr(lev_nbr)*0.
time_dom_avg_LWC=fltarr(lev_nbr)*0.
time=0.
;
end_time_min=long(end_time_min)
start_time_min=long(start_time_min)
incr_time_min=long(incr_time_min)
time_min=end_time_min
;
; Loop over the requested time period
;
for time_min=start_time_min,end_time_min,incr_time_min do begin
day=((time_min - (time_min mod 1440))/1440)+1;
minute_of_day=time_min-1440*(day-1);
hour_of_day=(minute_of_day-(minute_of_day mod 60))/60;
;
; Construct the file name
;
day_sng=string(format='(I2.2)',day)
minute_of_day_sng=string(format='(I4.4)',minute_of_day)
fl_in='cem.'+day_sng+'.'+minute_of_day_sng+'.nc'
;
; Read in binary data from the NetCDF file 
;
fl_path='/data/zender/_aux0_/cem'
fl_in=fl_path+'/'+fl_in
nc_id=ncdf_open(fl_in)
;
; get the two-dimensional arrays
;
ncdf_varget,nc_id,'conv_frac',conv_frac
ncdf_varget,nc_id,'dom_avg_CWC',dom_avg_CWC
ncdf_varget,nc_id,'dom_avg_IWC',dom_avg_IWC
ncdf_varget,nc_id,'dom_avg_LWC',dom_avg_LWC
;
; End of NetCDF commands
;
ncdf_close,nc_id
;
; Accumulate the data for later averaging, print file info. to stdout
;
time_dom_avg_CWC=time_dom_avg_CWC+dom_avg_CWC
time_dom_avg_IWC=time_dom_avg_IWC+dom_avg_IWC
time_dom_avg_LWC=time_dom_avg_LWC+dom_avg_LWC
time=[time,1.+time_min/1440.0] ; time is in days
;
print,'time = ',auto_sng(time_min,0),', Day = ',auto_sng(day,0),' minute = ',auto_sng(minute_of_day,0),', fraction 500 mb conv. active grdpts. = ',auto_sng(conv_frac(16)*100.0,1)
;
endfor ; end loop over times
;
; After the time loop, gather some useful information that
; does not change (too much!) with each passing file
;
nc_id=ncdf_open(fl_in)
ncdf_varget,nc_id,'w_thr',w_thr
ncdf_varget,nc_id,'dom_avg_p',dom_avg_p
ncdf_varget,nc_id,'dom_avg_t',dom_avg_t
ncdf_varget,nc_id,'lev',lev
ncdf_close,nc_id
if n_elements(dom_avg_p) ne lev_nbr then print,'Something wrong with levels'
;
time=time(1:time_nbr) ; truncate excess space--recall time(0) is garbage
;
time_dom_avg_CWC=time_dom_avg_CWC/time_nbr
time_dom_avg_IWC=time_dom_avg_IWC/time_nbr
time_dom_avg_LWC=time_dom_avg_LWC/time_nbr
;
time_dom_avg_CWC=time_dom_avg_CWC*1000.; convert to g/m^3
time_dom_avg_IWC=time_dom_avg_IWC*1000.; convert to g/m^3
time_dom_avg_LWC=time_dom_avg_LWC*1000.; convert to g/m^3
;
abc_1=time_dom_avg_CWC(0:lev_nbr-1)
abc_2=time_dom_avg_IWC(0:lev_nbr-1)
abc_3=time_dom_avg_LWC(0:lev_nbr-1)
;
;ord=dom_avg_p(0:lev_nbr-1)/100.
ord=dom_avg_t(0:lev_nbr-1)-273.16
;ord=lev/1000.
;
;ice_idx=where(abc_1 gt 0,count)
;abc_2=100.*abc_1/total(abc_1)
;abc_2=100.*time_dom_avg_CWC(0:lev_nbr-1)/total(time_dom_avg_CWC(0:lev_nbr-1))
;abc_2=100.*time_dom_avg_IWC(0:lev_nbr-1)/total(time_dom_avg_IWC(0:lev_nbr-1))
;abc_2=100.*time_dom_avg_LWC(0:lev_nbr-1)/total(time_dom_avg_LWC(0:lev_nbr-1))
;
erase
;
plot, $
	abc_1, $
	ord, $
	tit='!5Mean Condensate by Phase', $
;	xtit='!5!8CWC!5(!8p!5), !8IWC!5(!8p!5), !8LWC!5(!8p!5) (g m!E-3!N)', $
	xtit='!5!8CWC!5, !8IWC!5, !8LWC!5 (g m!E-3!N)', $
;	xtit='!5!8CWC!5(!8z!5), !8IWC!5(!8z!5), !8LWC!5(!8z!5) (g m!E-3!N)', $
;	xtit='!5Condensed Water Content !8CWC!5(!8z!5) (g m!E-3!N)', $
;	xtit='!5Ice Water Content !8IWC!5(!8z!5) (g m!E-3!N)', $
;	xtit='!5Liquid Water Content !8LWC!5(!8z!5) (g m!E-3!N)', $
;	ytit='!5Pressure !8p!5 (mb)', $
	ytit='!5Temperature !8T!5 (!E!12_!5!NC)', $
;	ytit='!5Height !8z!5 (km)', $
;	yrange=[1000.0,100.0], $
	yrange=[28.,-75.], $
;	yrange=[0.0,15.], $
	xstyle=0, $
	ystyle=1, $
	/ynozero, $
	thick=3.0, $
	charsize=2.0, $
	xmargin=[7.5,2], $ ;[10,3] is [left,right] default
	ymargin=[3.1,2], $  ;[4,2] is [bottom,top] default
;	ticklen=1.0, $
	linestyle=0
;
oplot, $
	abc_2, $
	ord, $
	thick=3.0, $
	linestyle=1
;
oplot, $
	abc_3, $
	ord, $
	thick=3.0, $
	linestyle=2
;
;plot, $
;	abc_2, $
;	dom_avg_p(0:lev_nbr-1)/100.0, $
;	yrange=[1000.0,100.0], $
;	xstyle=13, $
;	ystyle=12, $
;	/ynozero, $
;	thick=3.0, $
;	charsize=2.0, $
;	xmargin=[7.5,2], $ ;[10,3] is [left,right] default
;	ymargin=[3.1,2], $  ;[4,2] is [bottom,top] default
;	linestyle=1, $
;	/noerase, $
;	/nodata
;;
;axis,	$
;	xaxis=1, $
;	charsize=2.0, $
;	xtit='!5% of Total !8CWC!5', $
;;	xtit='!5% of Total !8IWC!5', $
;;	xtit='!5% of Total !8LWC!5', $
;	xmargin=[7.5,2], $ ;[10,3] is [left,right] default
;	ymargin=[3.1,2] ;[4,2] is [bottom,top] default
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.80
lgn_dy=0.075
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!8CWC!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!8IWC!5',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!8LWC!5',size=txt_lgn_sz,/NORMAL
;
t_500mb=-6.
;xyouts,!x.crange(1),t_500mb,'!8P!5 = 500 mb',size=txt_lgn_sz,alignment=1.0,/DATA
;plots,[!x.crange(0),!x.crange(1)],[t_500mb,t_500mb],linestyle=0,thick=3.0,/DATA
arrow,!x.crange(0)+1.05*(!x.crange(1)-!x.crange(0)),t_500mb,!x.crange(1),t_500mb,thick=2.0,/DATA
arrow,!x.crange(0)-0.05*(!x.crange(1)-!x.crange(0)),t_500mb,!x.crange(0),t_500mb,thick=2.0,/DATA
;
;xyouts,txt_lgn_x,lgn_y(3),'!5Days '+auto_sng(round(min(time)),0)+' - '+auto_sng(fix(max(time))),size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(4),'!8w!I!5min!N = '+auto_sng(w_thr,1)+' m/s',size=txt_lgn_sz,/NORMAL
;
end; end cem_iwc()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM IWC commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM Contour
; Contour plot the mass mixing ratio vs. altitude and longitude.
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..lev_nbr][1..nbr_sz] (in C)
; is accessed as         foo(1..nbr_sz,1..lev_nbr)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_contour_bch
open_ps,fl_nm='/data/zender/ps/cem_snapshot.eps',x_sz=6.5,y_sz=2.5,/eps,/color
cem_contour,fl_nm='/data/zender/_aux0_/cem/cem.01.1380.nc',scale=1000,palette=1,order=1,ttl='!5Condensate Mixing Ratio (g kg!E-1!N)',fld_nm='q_condensed',rng_y=[0,15],cntr_lvl_nbr=13,chr_sz=0.8,prn=1
close_ps,fl_nm='/data/zender/ps/cem_snapshot.eps'
end; end cem_contour_bch()
pro cem_contour, $
	fl_nm=fl_nm, $
	palette=palette, $
	order=order, $
	cntr_lvl_nbr=cntr_lvl_nbr, $
	rng_y=rng_y, $
	prn=prn, $
	ttl=ttl, $
	chr_sz=chr_sz, $
	scale=scale, $
	fld_nm=fld_nm
;
; cem_contour,ttl='!5Condensate Mixing Ratio (g kg!E-1!N)',scale=1000,palette=4,order=1,rng_y=[0,15]
;
@ibp_clr.com
;
if n_elements(scale) eq 0 then scale=1.
if n_elements(prn) eq 0 then prn=0
if n_elements(ttl) eq 0 then ttl=''
if n_elements(chr_sz) eq 0 then chr_sz=1.
if n_elements(cntr_lvl_nbr) eq 0 then cntr_lvl_nbr=22
if n_elements(rng_y) eq 0 then rng_y=0
if n_elements(order) ne 0 then color_order=order else color_order=1
if n_elements(palette) eq 0 then palette=3
if n_elements(fl_nm) eq 0 then fl_nm='/data/zender/_aux0_/cem/cem.01.1380.nc'
if n_elements(fld_nm) eq 0 then fld_nm='q_condensed'
;
if n_elements(clr_tbl) eq 0 then clr_tbl=34
clr_rfr,clr_tbl,"",0
;
erase
;
; Read in binary data from the NetCDF file 
;
nc_id=ncdf_open(fl_nm)
;
; get the scalars 
;
lev_id=ncdf_dimid(nc_id,'lev')
ncdf_diminq,nc_id,lev_id,foo,lev_nbr
lon_id=ncdf_dimid(nc_id,'lon')
ncdf_diminq,nc_id,lon_id,foo,lon_nbr
;
; get the one-dimensional arrays
;
ncdf_varget,nc_id,'lev',lev
ncdf_varget,nc_id,'lon',lon
;
; get the two-dimensional arrays
;
ncdf_varget,nc_id,fld_nm,tmp_data
;
; say bye-bye
ncdf_close,nc_id
; End of NetCDF commands
;
if n_elements(tmp_data) le 1 then begin
	print,'Requested data not in NetCDF file, returning.'
	return
endif
;
;ncdf_varget,nc_id,'dom_avg_p',dom_avg_p
;ncdf_varget,nc_id,'dom_avg_t',dom_avg_t
;print,format='(4(a8,1x))','level','z','p','T'
;for level=0,lev_nbr-1 do begin
;	print,format='(4(a8,1x))',auto_sng(level,0),auto_sng(lev(level)/1000,2),auto_sng(dom_avg_p(level)/100,2),auto_sng(dom_avg_t(level)-273.15,2)
;endfor
;
tmp_data=reform(tmp_data)
data=fltarr(lon_nbr,lev_nbr)
for level=0,lev_nbr-1 do begin
	data(*,level)=tmp_data(level,*)
endfor
data=data*scale
abc=lon/1000.
ord=lev/1000.
if n_elements(ttl) le 0 then ttl=fld_nm
x_axz_ttl='Longitude !8x!5 (km)'
y_axz_ttl='Altitude !8z!5 (km)'
abc_min=0.
abc_max=906.
if n_elements(rng_y) le 1 then begin
	ord_min=0.
	ord_max=20.
endif else begin
	ord_min=rng_y(0)
	ord_max=rng_y(1)
endelse
;
grd_ntv=0.
;cntr_lvl_min=0.
cntr_lvl_min=min(data)
;first_level=-1.
;
data_min=min(data)
data_max=max(data)
;max_level=(data_max+grd_ntv) - abs(data_max mod grd_ntv)
max_level=data_max
cntr_ntv=(max_level-cntr_lvl_min)/(cntr_lvl_nbr-2)
cntr_lvl=fltarr(cntr_lvl_nbr)
cntr_lvl(0)=cntr_lvl_min
;if first_level lt data_min then first_level=cntr_ntv
;cntr_lvl(1)=first_level
for i=1,cntr_lvl_nbr-1 do begin
	cntr_lvl(i)=cntr_lvl(i-1)+cntr_ntv
endfor
cntr_which_lbl=indgen(cntr_lvl_nbr)*0	;which levels to label
cntr_ntt=auto_sng(cntr_lvl,2)	;actual text labels
cbar_lgn=cntr_ntt
;
print,'maximum data = ',data_max
print,'maximum contour = ',cntr_lvl(cntr_lvl_nbr-1)
;
cbar_clr_nbr=min([!d.table_size-1,cntr_lvl_nbr-1])
cntr_fll_idx=indgen(cbar_clr_nbr)+2
cbar_idx=indgen(cbar_clr_nbr)+2
;
cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=0.8*chr_sz
cbar_unit=''
cbar_lbl_sz=0.8*chr_sz
;
if prn then begin
plt_rgn_nrm=[.08,.15,0.9,0.9]; [x_min,y_min,x_max,y_max]
cbar_psn=[.91,.15,0.95,0.9]; [x_min,y_min,x_max,y_max]
endif else begin; not prn
plt_rgn_nrm=[0.1,.1,0.87,0.90]; [x_min,y_min,x_max,y_max]
cbar_psn=[.88,0.10,0.92,0.90]; [x_min,y_min,x_max,y_max]
endelse; not prn
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
clr_mk,cbar_clr_nbr,color_order,palette,1
;
tvlct,r_curr,g_curr,b_curr
;
clr_bar_drw, $
	bar_psn=cbar_psn, $
	bar_clr_nbr=cbar_clr_nbr, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz
;
contour, $
	data, $
	abc, $
	ord, $
	tit=ttl, $
	xtit=x_axz_ttl, $
	ytit=y_axz_ttl, $
	level=cntr_lvl, $                 
	c_color=cntr_fll_idx, $		
	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
;	font=0, $				
	xrange=[abc_min,abc_max], $
	yrange=[ord_min,ord_max], $
	xstyle=1, $
	ystyle=1, $
	charsize=chr_sz, $
	position=plt_rgn_nrm, $
	ycharsize=1.1, $
	xcharsize=1.1, $
	/closed, $				
	/fill, $				
	/noerase
;
;contour, $
;	data, $
;	abc, $
;	ord, $
;	level=cntr_lvl, $                 
;	c_labels=cntr_which_lbl, $	
;;	c_annotation=cntr_ntt, $	
;	/closed, $				
;	/overplot
;
end; end cem_contour()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM Contour commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure CEM Animate
; Analyze the time series of (semi) arbitrary CEM fields.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cem_anim, $
	start_time_min=start_time_min, $
	end_time_min=end_time_min, $
	incr_time_min=incr_time_min, $
	fld_nm=fld_nm, $
	palette=palette, $
	order=order, $
	rng_y=rng_y, $
	ttl=ttl, $
	scale=scale, $
	rate=rate
;
; Remember not to try to contour the condensate fields at time 0 
; because they're all 0 and will cause an error.
;
; cem_anim,start_time_min=4320,end_time_min=34500,incr_time_min=60
; cem_anim,start_time_min=60,end_time_min=1430,incr_time_min=60
; cem_anim,fld_nm='q_condensed'
;
if n_elements(start_time_min) eq 0 then start_time_min=60
if n_elements(end_time_min) eq 0 then end_time_min=1430
if n_elements(incr_time_min) eq 0 then incr_time_min=60
if n_elements(fld_nm) eq 0 then fld_nm='q_condensed'
if n_elements(order) ne 0 then color_order=order else color_order=1
if n_elements(palette) eq 0 then palette=3
if n_elements(rate) eq 0 then rate=25
if n_elements(scale) eq 0 then scale=1.
if n_elements(ttl) eq 0 then ttl=''
if n_elements(rng_y) eq 0 then rng_y=0
;if n_elements(foo) eq 0 then foo=''
;
; Loop over the requested time period
;
cntr_wdw_x_sz=1200
cntr_wdw_y_sz=500
;
end_time_min=long(end_time_min)
start_time_min=long(start_time_min)
incr_time_min=long(incr_time_min)
time_min=end_time_min
;
time_nbr=(end_time_min-start_time_min)/incr_time_min+1
window,1,xsize=cntr_wdw_x_sz,ysize=cntr_wdw_y_sz
wset,1
xinteranimate,set=[cntr_wdw_x_sz,cntr_wdw_y_sz,time_nbr],/track
;
time_idx=0
time=fltarr(time_nbr)
for time_min=start_time_min,end_time_min,incr_time_min do begin
day=((time_min - (time_min mod 1440))/1440)+1;
minute_of_day=time_min-1440*(day-1);
hour_of_day=(minute_of_day-(minute_of_day mod 60))/60;
;
; Construct the file name
;
day_sng=string(format='(I2.2)',day)
minute_of_day_sng=string(format='(I4.4)',minute_of_day)
fl_nm='cem.'+day_sng+'.'+minute_of_day_sng+'.nc'
;
; Call the contouring routine
;
fl_path='/data/zender/_aux0_/cem'
fl_nm=fl_path+'/'+fl_nm
;
print,'time = '+auto_sng(time_min,0)+' min., '+fl_nm+': '+fld_nm
;
cem_contour,fl_nm=fl_nm,palette=palette,order=order,rng_y=rng_y,scale=scale,ttl=ttl,fld_nm=fld_nm
;
xinteranimate,frame=time_idx,window=1
;
time(time_idx)=time_min
;
time_idx=time_idx+1
;
endfor ; end loop over times
;
xinteranimate,rate
xinteranimate,/close
;
end; end cem_anim()
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure CEM Animate commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

