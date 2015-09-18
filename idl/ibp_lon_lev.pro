; $Id$

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin lon_lev() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lon_lev_event,event

@ibp_cmn.com

event_type=tag_names(event,/structure)
;print,event_type,event

; In case event is not handled by the below
command_sng=''

; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value

if (usr_value eq 'Wlon_lev_gph') then begin
	data_crd=convert_coord( $
		event.x, $
		event.y, $
		/device, $
		/to_data)
	if ((event.type eq 0) and (event.press eq 2)) then begin
; middle button was pressed
		ptr_psn(0)=data_crd(0)
		ptr_psn(1)=data_crd(1)
	endif
	if ((event.type eq 1) and (event.release eq 2)) then begin
; middle button was released
		ptr_psn(2)=data_crd(0)
		ptr_psn(3)=data_crd(1)
		inverse_lon_lev
	endif
endif; endif Wlon_lev_gph

if command_sng ne '' then print,command_sng

end; end lon_lev() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End lon_lev() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin make_lon_lev() 
; If the current device is X then this routine will create the 
; widgets and install the event handler for the Longitude Level graphs if
; they don't already exist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_lon_lev

@ibp_cmn.com

if !d.name eq 'X' then begin

; If the window is already up then OK else start it up.
; be careful because old widget ID's stay in common memory 
; after crashes!

widget_xst=1
if n_elements(Wlon_lev_gph) ne 0 then begin
	widget_xst=widget_info(Wlon_lev_gph,/valid_id)
endif

if ((widget_xst) and (n_elements(Wlon_lev_gph) ne 0)) then begin
curr_sz=indgen(2)
widget_control,Wlon_lev_base,tlb_get_size=curr_sz

; If the window has been resized, then destroy the widget hierarchy
; so that the widget can be recreated at the new size.
wset,lon_lev_wdw

if ((lon_lev_wdw_x_sz ne curr_sz(0)) or (lon_lev_wdw_y_sz ne curr_sz(1))) then begin
	lon_lev_wdw_x_sz=curr_sz(0)
	lon_lev_wdw_y_sz=curr_sz(1)
	widget_control,Wlon_lev_base,/destroy
	widget_xst=0
endif

endif

if ((widget_xst ne 1) or (n_elements(Wlon_lev_gph) eq 0)) then begin

; Create the Longitude Level widget hierarchy
Wlon_lev_base=widget_base( $
	group_leader=Wtop_lvl, $
	title='Longitude Level')
Wlon_lev_gph=widget_draw( $
	Wlon_lev_base, $
		xsize=lon_lev_wdw_x_sz, $
		ysize=lon_lev_wdw_y_sz, $
	/button_events, $
	uvalue='Wlon_lev_gph')

widget_control,Wlon_lev_base,/realize
widget_control,get_value=lon_lev_wdw,Wlon_lev_gph

; Call the event manager for the new base widgets
xmanager, $
	'make_lon_lev', $
	Wlon_lev_base, $
	group_leader=Wtop_lvl, $
	event_handler='lon_lev_event'		

endif; end creating widgets

widget_control,Wlon_lev_base,iconify=0 ; de-iconify when drawing
wset,lon_lev_wdw

endif; endif device is X

end; end make_lon_lev() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End make_lon_lev() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure Longitude Level
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lon_lev
@ibp_cmn.com

erase

; Read in binary data from the NetCDF file 
nc_id=ncdf_open(fl_in)

; get scalars 
lev_id=ncdf_dimid(nc_id,lev_dim_nm)
ncdf_diminq,nc_id,lev_id,foo,lev_nbr
lat_id=ncdf_dimid(nc_id,'lat')
ncdf_diminq,nc_id,lat_id,foo,lat_nbr
lon_id=ncdf_dimid(nc_id,'lon')
ncdf_diminq,nc_id,lon_id,foo,lon_nbr

; get one-dimensional arrays
ncdf_varget,nc_id,lev_dim_nm,lev
ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,'lon',lon

; due to ECMWF idiosyncracies, ECMWF variables are defined on 
; interfaces, not levels, so swap the aint data into the lev 
; array for ECMWF
result=strpos(fl_in,'ecmwf')
if result ne -1 then begin
result=strpos(fl_in,'omega')
if result eq -1 then begin
	ncdf_attget,nc_id,'alev',alev,/GLOBAL
	ncdf_attget,nc_id,'amid',amid,/GLOBAL
	ncdf_attget,nc_id,'bmid',bmid,/GLOBAL
	ncdf_attget,nc_id,'aint',aint,/GLOBAL
	ncdf_attget,nc_id,'bint',bint,/GLOBAL
	level_pressure=amid*1000+bmid*1013
	interface_pressure=aint*1000+bint*1013
	print,'lev = ',lev
	print,'level_pressure = ',level_pressure
	print,'interface_pressure = ',interface_pressure
	lev=interface_pressure(1:lev_nbr)
endif; endif we're looking at ECMWF data
endif; endif we're looking at ECMWF data

; get three-dimensional arrays
ncdf_varget,nc_id,fld_nm,tmp_data

; say bye-bye
ncdf_close,nc_id
; End of NetCDF commands

tmp_data=reform(tmp_data)

if n_elements(tmp_data) le 1 then begin
	print,'Requested data not in netCDF file, returning.'
	goto,end_of_procedure
endif

dim_nbr=(size(tmp_data))(0)

if dim_nbr eq 4 then begin
; Find level and time indices
	lev_idx=-1
	for dim_idx=0,dim_nbr-1 do begin
		if dim_lst(fld_lst(fld_idx).dim(dim_idx)).name eq 'time' then time_idx=dim_idx
		if dim_lst(fld_lst(fld_idx).dim(dim_idx)).name eq 'lev' then lev_idx=dim_idx
	endfor
	if time_idx eq -1 then print,'ERROR: dimension time not found'
	if lev_idx eq -1 then print,'ERROR: dimension lev not found'

	print,'Time dimension is index = ',time_idx
	print,'Extracting time ',time_slc,' which is at ',time(time_slc)
	if time_idx eq 0 then begin
		tmp_data=tmp_data(time_slc,*,*,*)
	endif else if time_idx eq 1 then begin
		tmp_data=tmp_data(*,time_slc,*,*)
	endif else if time_idx eq 2 then begin
		tmp_data=tmp_data(*,*,time_slc,*)
	endif else if time_idx eq 3 then begin
		tmp_data=tmp_data(*,*,*,time_slc)
	endif else begin
		print,'Unprogrammed time slice'	
	endelse ; endif time_idx
	tmp_data=reform(tmp_data)

endif; endif 4 dimensions

dim_nbr=(size(tmp_data))(0)
if dim_nbr ne 3 then begin
	print,'Need 3 and only 3 dimensions, returning.'
	goto,end_of_procedure
endif

; Find level index
lev_idx=-1
for dim_idx=0,dim_nbr-1 do begin
	if dim_lst(fld_lst(fld_idx).dim(dim_idx)).name eq 'lev' then lev_idx=dim_idx
endfor
if lev_idx eq -1 then print,'ERROR: dimension lev not found'

; Put dimensions in canonical foo(lon,lev,lat) order in IDL notation.
if lev_idx ne 1 then print,'Changing dimension order to lon,lev,lat'

if lev_idx eq 2 then begin
; Change from lon,lat,lev to lon,lev,lat order
	tmp_data=transpose(tmp_data,[0,2,1])
endif; endif lev_idx eq 2

print,'Beginning data reformatting...'

if lon(lon_nbr-1) gt 90.0 and lon(0) ge 0.0 then begin
	if lon(0) eq 0.0 then nbr_shift=lon_nbr/2-1 else nbr_shift=lon_nbr/2
	lon=shift(lon,nbr_shift)
	west_lon_idx=where(lon gt 180.0,count)	
	if count gt 0 then lon(west_lon_idx)=-360.0+lon(west_lon_idx)
	tmp_data=shift(tmp_data,nbr_shift,0,0)
endif; endif lon seems to be in the right range

print,'Beginning Region Truncation...'

foo=min(abs(lon-lon_min),lon_min_idx)
foo=min(abs(lon-lon_max),lon_max_idx)
foo=min(abs(lat-lat_min),lat_min_idx)
foo=min(abs(lat-lat_max),lat_max_idx)
foo=min(abs(lev-lev_min),lev_min_idx)
foo=min(abs(lev-lev_max),lev_max_idx)

; If the data isn't stored in CCM altitude order (surface pressure
; stored at last index...), then re-order that dimension because
; the truncation routine needs lev_min_idx < lev_max_idx.

if lev_min_idx gt lev_max_idx then begin
	tmp_data=reverse(tmp_data,2)
	lev=reverse(lev,1)
	swap=lev_min_idx
	lev_min_idx=lev_max_idx
	lev_max_idx=swap
endif; endif reversing level order

print,'lon_nbr = ',lon_nbr,' min(lon) = ',min(lon),' max(lon) = ',max(lon)
print,'lon_min = ',lon_min,' lon_min_idx = ',lon_min_idx,' = ',lon(lon_min_idx)
print,'lon_max = ',lon_max,' lon_max_idx = ',lon_max_idx,' = ',lon(lon_max_idx)
print,''
print,'lat_nbr = ',lat_nbr,' min(lat) = ',min(lat),' max(lat) = ',max(lat)
print,'lat_min = ',lat_min,' lat_min_idx = ',lat_min_idx,' = ',lat(lat_min_idx)
print,'lat_max = ',lat_max,' lat_max_idx = ',lat_max_idx,' = ',lat(lat_max_idx)
print,''
print,'lev_nbr = ',lev_nbr,' min(lev) = ',min(lev),' max(lev) = ',max(lev)
print,'lev_min = ',lev_min,' lev_min_idx = ',lev_min_idx,' = ',lev(lev_min_idx)
print,'lev_max = ',lev_max,' lev_max_idx = ',lev_max_idx,' = ',lev(lev_max_idx)

; Truncate data to specified region
; Dimensions are assumed to be foo(lon,lev,lat) in IDL notation.

; Is there to be any level regionalization?
if ((lev_min_idx ne 0) or (lev_max_idx ne lev_nbr-1)) then begin
; Truncate the data outside the desired level region
	tmp_data=tmp_data(*,lev_min_idx:lev_max_idx,*)
	lev=lev(lev_min_idx:lev_max_idx)
endif; endif level needed truncating

; Is there to be any meridional regionalization?
if ((lat_min_idx ne 0) or (lat_max_idx ne lat_nbr-1)) then begin
; Truncate the data outside the desired longitude region
	tmp_data=tmp_data(*,*,lat_min_idx:lat_max_idx)
	lat=lat(lat_min_idx:lat_max_idx)
endif; endif longitude needed truncating

; Is there to be any zonal regionalization?
if ((lon_min_idx ne 0) or (lon_max_idx ne lon_nbr-1)) then begin

; If the zonal range crosses the date line.....
if lon_min gt lon_max then begin
; Find the current index of the desired left-most longitude
	shift_amt=lon_nbr-lon_min_idx
	tmp_data=shift(tmp_data,shift_amt,0,0)	
	lon=shift(lon,shift_amt)	

; Truncate the data to the right of the desired right-most longitude
	tmp_data=tmp_data(0:lon_max_idx+shift_amt,*,*)
	lon=lon(0:lon_max_idx+shift_amt)

endif; endif data region crossed the date line

; If the zonal range does not cross the date line.....
if lon_min lt lon_max then begin
; Truncate the data outside the desired longitude region
	tmp_data=tmp_data(lon_min_idx:lon_max_idx,*,*)
	lon=lon(lon_min_idx:lon_max_idx)

endif; endif data region did not cross the date line

endif; endif longitude needed regionalization

lat_nbr=n_elements(lat)
lev_nbr=n_elements(lev)
lon_nbr=n_elements(lon)
data_nbr=n_elements(tmp_data)

print,'Post rgn_trunc....'
print,'lon_nbr = ',lon_nbr
print,'lev_nbr = ',lev_nbr
print,'lat_nbr = ',lat_nbr
print,'data_nbr = ',data_nbr
print,'size(tmp_data) = ',size(tmp_data)

abc=lon
ord=lev

if (strlen(unit_jd(fgr_idx)) le 0) then unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
if (strlen(ttl_jd(fgr_idx)) le 0) then ttl_jd(fgr_idx)=src_sng+' '+fld_sng+', '+month_sng
if (strlen(xttl_jd(fgr_idx)) le 0) then xttl_jd(fgr_idx)='Longitude !7k!5 (!E!12_!5!N)'
if (strlen(yttl_jd(fgr_idx)) le 0) then yttl_jd(fgr_idx)='Hybrid Pressure !7g!5 (mb)'

abc_min=lon_min
abc_max=lon_max
ord_min=lev_min
ord_max=lev_max

; This nifty matrix trick fills each of the lat_nbr rows of the result
; with the corresponding element of lat_wgt. There are lon_nbr
; elements per row in the result.
lat_wgt_rgn=lat_wgt(lat_min_idx:lat_max_idx)

; Do meridional averaging 
data=fltarr(lon_nbr,lev_nbr)
good_data_nbr=0
for lon_idx=0,lon_nbr-1 do begin
for lev_idx=0,lev_nbr-1 do begin
	good_idx=where(tmp_data(lon_idx,lev_idx,*) lt very_big_nbr,count)
	vld_cnt_pct=100.0*count/lat_nbr
	if vld_cnt_pct lt vld_avg_thr_pct then begin
		print,'WARNING: Only '+auto_sng(vld_cnt_pct,1)+'% of lats valid in meridional average of '+fld_nm+' at lon('+auto_sng(lon_idx,0)+') = '+auto_sng(lon(lon_idx),1)+' degrees.'
	endif
	if count ne 0 then begin
		good_data=tmp_data(lon_idx,lev_idx,good_idx)
		good_lat_wgt_rgn=lat_wgt_rgn(good_idx)
		data(lon_idx,lev_idx)=total(good_data*good_lat_wgt_rgn)/total(good_lat_wgt_rgn)
		good_data_nbr=good_data_nbr+count
	endif else begin
		data(lon_idx,lev_idx)=1.0e36
	endelse
endfor; end loop over lev
endfor; end loop over lon

data=rotate(data,2)
ord=rotate(ord,2)
abc=rotate(abc,2)

scale_data
data_max=max(data)
data_min=min(data)

cbar_fmt=replicate('UNSET',1)
if not usr_cntr_jd(fgr_idx) then cntr_lvl=0.
form_cntr_lvl
comp_cntr_lvl
;bld_sngs
cbar_setup

cbar_psn=[ $
	0.88, $ ; x_min
	0.10, $ ; y_min
	0.92, $ ; x_max
	0.90] ; 7_max

plt_rgn_nrm=[ $
	0.10, $ ; x_min
	0.10, $ ; y_min
	0.87, $ ; x_max
	0.90] ; y_max

clr_mk,cbar_clr_nbr,color_order,pll_idx,background_jd(fgr_idx)

tvlct,r_curr,g_curr,b_curr

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

if !d.name ne 'X' then ttl_sz=1.6 else ttl_sz=1.0
if !d.name ne 'X' then chr_sz_y=1.1 else chr_sz_y=0.8
if !d.name ne 'X' then chr_sz_x=1.1 else chr_sz_x=0.8
 
contour, $
	data, $
	abc, $
	ord, $
	tit=ttl_jd(fgr_idx), $
	xtit=xttl_jd(fgr_idx), $
	ytit=yttl_jd(fgr_idx), $
	level=cntr_lvl, $                 
	c_color=cntr_fll_idx, $		
	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
;	font=0, $				
	xrange=[abc_min,abc_max], $
	yrange=[ord_max,ord_min], $
	xstyle=13, $
	ystyle=1, $
	charsize=ttl_sz, $
	position=plt_rgn_nrm, $
	ycharsize=chr_sz_y, $
	xcharsize=chr_sz_x, $
;	/closed, $				
	/fill, $				
	/noerase

; NB: Until IDL 4.0, this command was not needed because map_set() set x.crange and y.crange, which the axes routines need to do correct labeling
!x.range=[abc_min,abc_max]

lon_axz_drw, $
	lon_top=abc_max, $
	lon_btm=abc_min, $
	axz_vrt=0

if cntr_ovl_jd(fgr_idx) then $
contour, $
	data, $
	abc, $
	ord, $
	level=cntr_lvl, $                 
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $	
	c_annotation=cntr_ntt, $	
;	/closed, $				
	/overplot

lev_min_sng=auto_sng(round(lev_min),0)
lev_max_sng=auto_sng(round(lev_max),0)
lev_sng='!5'+lev_min_sng+' mb < !7g!5 < '+ $
		lev_max_sng+' mb'

if (info_jd(fgr_idx) gt 0) then begin
	xyouts,0.99,.025,sub_ttl_sng+' '+sym_sng+' in '+rgn_lst(rgn_idx).string+': '+lon_sng+', '+lat_sng+', '+lev_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
;	xyouts,1.0,.025,'Used '+auto_sng(good_data_nbr,0)+'/'+auto_sng(data_nbr,0)+' points.',alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
endif; endif info

end_of_procedure: foo=1
end; end lon_lev()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure Longitude Level
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
