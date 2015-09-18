;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CVS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: ibp_time_lat.pro,v $
; $Source: /home/zender/cvs/idl/ibp_time_lat.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; $Log: not supported by cvs2svn $
; Revision 1.14  2000-03-11 00:21:17  zender
; *** empty log message ***
;
; Revision 1.13  2000/03/01 02:28:35  zender
; *** empty log message ***
;
; Revision 1.12  2000/01/15 02:07:53  zender
; *** empty log message ***
;
; Revision 1.8  2000/01/01 01:55:52  zender
; *** empty log message ***
;
; Revision 1.6  1999/12/31 02:09:42  zender
; *** empty log message ***
;
; Revision 1.5  1999/12/31 00:18:27  zender
; *** empty log message ***
;
; Revision 1.4  1999/12/30 19:34:26  zender
; *** empty log message ***
;
; Revision 1.3  1999/10/12 07:37:43  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/03 16:52:04  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:50  zender
; Imported sources
;
; Revision 3.8  1998/10/19 19:01:36  zender
; maintenance check-in
;
; Revision 3.7  1995/09/22  17:43:49  zender
; removed lon_lat map from main screen, implemented adjustable
; window sizes. still have bug in color table 4.
;
; Revision 3.6  1995/06/15  17:31:35  zender
; added azimuthal polar plots, changed distribution slightly.
;
; Revision 3.5  1995/05/09  22:55:21  zender
; finished merging doetzl's mods. added a customize screen to
; each figure type. removed buttons from main window. added
; user-specified contour levels. merged contour level database
; into pre_fld structure.
;
; Revision 3.4  1995/05/08  17:14:46  zender
; finished merge with doetzl's customization mods. about to
; bring customization to all cases.
;
; Revision 3.3  1995/05/07  22:59:09  zender
; about to merge in doetzl's mods.
;
; Revision 3.2  1995/04/07  21:38:54  zender
; first version for general release. renamed everything ibp
; instead of ccm.
;
; Revision 3.1  1995/04/07  21:17:20  zender
; final version before ccm -> ibp name switch. diddled with
; color routines to add a black/white palette. switched off
; user input of cntr_lvls in ibp_lat_lev but now allow exact
; contour specs. in ibp_fgr for pre-defined fields.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin time_lat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro time_lat_event,event
;
@ibp_cmn.com
;
event_type=tag_names(event,/structure)
;print,event_type,event
;
; In case event is not handled by the below
command_sng=''
;
; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value
;
if (usr_value eq 'Wtime_lat_gph') then begin
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
		inverse_time_lat
	endif
endif; endif Wtime_lat_gph
;
if command_sng ne '' then print,command_sng
;
end; end time_lat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End time_lat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin make_time_lat() 
; If the current device is X then this routine will create the 
; widgets and install the event handler for the Time Latitude graphs if
; they don't already exist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_time_lat
;
@ibp_cmn.com
;
if !d.name eq 'X' then begin
;
; If the window is already up then OK else start it up.
; be careful because old widget ID's stay in common memory 
; after crashes!
widget_xst=1
if n_elements(Wtime_lat_gph) ne 0 then begin
	widget_xst=widget_info(Wtime_lat_gph,/valid_id)
endif
;
if ((widget_xst) and (n_elements(Wtime_lat_gph) ne 0)) then begin
curr_sz=indgen(2)
widget_control,Wtime_lat_base,tlb_get_size=curr_sz
;
; If the window has been resized, then destroy the widget hierarchy
; so that the widget can be recreated at the new size.
;
wset,time_lat_wdw
;
if ((time_lat_wdw_x_sz ne curr_sz(0)) or (time_lat_wdw_y_sz ne curr_sz(1))) then begin
	time_lat_wdw_x_sz=curr_sz(0)
	time_lat_wdw_y_sz=curr_sz(1)
	widget_control,Wtime_lat_base,/destroy
	widget_xst=0
endif
;
endif
;
if ((widget_xst ne 1) or (n_elements(Wtime_lat_gph) eq 0)) then begin
; Create the Time Latitude widget hierarchy
Wtime_lat_base=widget_base( $
	group_leader=Wtop_lvl, $
	title='Time Latitude')
Wtime_lat_gph=widget_draw( $
	Wtime_lat_base, $
		xsize=time_lat_wdw_x_sz, $
		ysize=time_lat_wdw_y_sz, $
	/button_events, $
	uvalue='Wtime_lat_gph')
;
widget_control,Wtime_lat_base,/realize
widget_control,get_value=time_lat_wdw,Wtime_lat_gph
;
; Call the event manager for the new base widgets
xmanager, $
	'make_time_lat', $
	Wtime_lat_base, $
	group_leader=Wtop_lvl, $
	event_handler='time_lat_event'		
;
endif; end creating widgets
;
widget_control,Wtime_lat_base,iconify=0 ; de-iconify when drawing
wset,time_lat_wdw
;
endif; endif device is X
;
end; end make_time_lat() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End make_time_lat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure Time Latitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro time_lat
@ibp_cmn.com
;
erase
;
; Read in binary data from the NetCDF file 
;
nc_id=ncdf_open(fl_in)
;
; get the scalars 
;
if lev_dim_nm ne '' then begin
	lev_id=ncdf_dimid(nc_id,lev_dim_nm)
	ncdf_diminq,nc_id,lev_id,foo,lev_nbr
endif else lev_nbr=0
lat_id=ncdf_dimid(nc_id,'lat')
ncdf_diminq,nc_id,lat_id,foo,lat_nbr
lon_id=ncdf_dimid(nc_id,'lon')
ncdf_diminq,nc_id,lon_id,foo,lon_nbr
time_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,time_id,foo,time_nbr
;
; get the one-dimensional arrays
;
if lev_nbr ne 0 then ncdf_varget,nc_id,lev_dim_nm,lev
ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,'lon',lon
ncdf_varget,nc_id,'time',time
;
; due to ECMWF idiosyncracies, ECMWF variables are defined on 
; interfaces, not levels, so swap the aint data into the lev 
; array for ECMWF
;
result=strpos(fl_in,'ecmwf')
if result ne -1 then begin
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
;
; get the three/four dimensional arrays
;
ncdf_varget,nc_id,fld_nm,tmp_data
;
; say bye-bye
ncdf_close,nc_id
; End of NetCDF commands
;
tmp_data=reform(tmp_data)
;
if n_elements(tmp_data) le 1 then begin
	print,'Requested data not in netCDF file, returning.'
	goto,end_of_procedure
endif
;
dim_nbr=(size(tmp_data))(0)
if dim_nbr eq 4 then begin
;
; Extract the correct vertical slice
;
	tmp_data=tmp_data(*,vrt_slc,*,*)
	tmp_data=reform(tmp_data)
endif
;
dim_nbr=(size(tmp_data))(0)
if dim_nbr ne 3 then begin
	print,'Need 3 and only 3 dimensions, returning.'
	goto,end_of_procedure
endif
;
print,'Beginning data reformatting...'
;
if lon(lon_nbr-1) gt 90. and lon(0) ge 0. then begin
	if lon(0) eq 0. then nbr_shift=lon_nbr/2-1 else nbr_shift=lon_nbr/2
	lon=shift(lon,nbr_shift)
	west_lon_idx=where(lon gt 180.0,count)	
	if count gt 0 then lon(west_lon_idx)=-360.+lon(west_lon_idx)
	tmp_data=shift(tmp_data,nbr_shift,0,0)
endif; endif lon seems to be in the right range
;
print,'Beginning Region Truncation...'
;
; DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 
;
time_min=time(0)
time_max=time(time_nbr-1)
;
foo=min(abs(lon-lon_min),lon_min_idx)
foo=min(abs(lon-lon_max),lon_max_idx)
foo=min(abs(lat-lat_min),lat_min_idx)
foo=min(abs(lat-lat_max),lat_max_idx)
foo=min(abs(time-time_min),time_min_idx)
foo=min(abs(time-time_max),time_max_idx)
;
; DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 
;
time_min_idx=0
time_max_idx=time_nbr-1
;
if time_min_idx gt time_max_idx then begin
	tmp_data=reverse(tmp_data,2)
	time=reverse(time,1)
	swap=time_min_idx
	time_min_idx=time_max_idx
	time_max_idx=swap
endif; endif reversing time order
;
print,'lon_nbr = ',lon_nbr,' min(lon) = ',min(lon),' max(lon) = ',max(lon)
print,'lon_min = ',lon_min,' lon_min_idx = ',lon_min_idx,' = ',lon(lon_min_idx)
print,'lon_max = ',lon_max,' lon_max_idx = ',lon_max_idx,' = ',lon(lon_max_idx)
print,''
print,'lat_nbr = ',lat_nbr,' min(lat) = ',min(lat),' max(lat) = ',max(lat)
print,'lat_min = ',lat_min,' lat_min_idx = ',lat_min_idx,' = ',lat(lat_min_idx)
print,'lat_max = ',lat_max,' lat_max_idx = ',lat_max_idx,' = ',lat(lat_max_idx)
print,''
print,'time_nbr = ',time_nbr,' min(time) = ',min(time),' max(time) = ',max(time)
print,'time_min = ',time_min,' time_min_idx = ',time_min_idx,' = ',time(time_min_idx)
print,'time_max = ',time_max,' time_max_idx = ',time_max_idx,' = ',time(time_max_idx)
;
; Truncate the data to the specified region
; Dimensions are assumed to be foo(lon,lat,time) in IDL notation.
;
; Is there to be any time regionalization?
if ((time_min_idx ne 0) or (time_max_idx ne time_nbr-1)) then begin
; Truncate the data outside the desired time region
	tmp_data=tmp_data(*,*,time_min_idx:time_max_idx)
	time=time(time_min_idx:time_max_idx)
endif; endif time needed truncating
;
; Is there to be any meridional regionalization?
if ((lat_min_idx ne 0) or (lat_max_idx ne lat_nbr-1)) then begin
; Truncate the data outside the desired latitude region
	tmp_data=tmp_data(*,lat_min_idx:lat_max_idx,*)
	lat=lat(lat_min_idx:lat_max_idx)
endif; endif latitude needed truncating
;
; Is there to be any zonal regionalization?
if ((lon_min_idx ne 0) or (lon_max_idx ne lon_nbr-1)) then begin
;
; If the zonal range crosses the date line.....
if lon_min gt lon_max then begin
; Find the current index of the desired left-most longitude
	shift_amt=lon_nbr-lon_min_idx
	tmp_data=shift(tmp_data,shift_amt,0,0)	
	lon=shift(lon,shift_amt)	
;
; Truncate the data to the right of the desired right-most longitude
	tmp_data=tmp_data(0:lon_max_idx+shift_amt,*,*)
	lon=lon(0:lon_max_idx+shift_amt)
;
endif; endif data region crossed the date line
;
; If the zonal range does not cross the date line.....
if lon_min lt lon_max then begin
; Truncate the data outside the desired longitude region
	tmp_data=tmp_data(lon_min_idx:lon_max_idx,*,*)
	lon=lon(lon_min_idx:lon_max_idx)
;
endif; endif data region did not cross the date line
;
endif; endif longitude needed regionalization
;
lat_nbr=n_elements(lat)
time_nbr=n_elements(time)
lon_nbr=n_elements(lon)
data_nbr=n_elements(tmp_data)
;
print,'Post rgn_trunc....'
print,'lon_nbr = ',lon_nbr
print,'lat_nbr = ',lat_nbr
print,'time_nbr = ',time_nbr
print,'data_nbr = ',data_nbr
print,'size(tmp_data) = ',size(tmp_data)
;
; The abc should be labelled by the month name.
;
;abc=time
abc=indgen(time_nbr)+.5
ord=lat
;
if (strlen(unit_jd(fgr_idx)) le 0) then unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
if (strlen(ttl_jd(fgr_idx)) le 0) then ttl_jd(fgr_idx)=src_sng+' '+fld_sng+', '+month_sng
if (strlen(xttl_jd(fgr_idx)) le 0) then xttl_jd(fgr_idx)='Time'
if (strlen(yttl_jd(fgr_idx)) le 0) then yttl_jd(fgr_idx)='Latitude !7u!5 (!E!12_!5!N)'
;
abc_min=time_min
abc_max=time_max
ord_min=lat_min
ord_max=lat_max
;
data=fltarr(time_nbr,lat_nbr)
good_data_nbr=0
for time_idx=0,time_nbr-1 do begin
for lat_idx=0,lat_nbr-1 do begin
	good_idx=where(tmp_data(*,lat_idx,time_idx) lt very_big_nbr,count)
	vld_cnt_pct=100.*count/lon_nbr
	if vld_cnt_pct lt vld_avg_thr_pct then begin
		print,'WARNING: Only '+auto_sng(vld_cnt_pct,1)+'% of lons valid in zonal average of '+fld_nm+' in lat('+auto_sng(lat_idx,0)+') = '+auto_sng(lat(lat_idx),1)+' degrees.'
	endif
	if count ne 0 then begin
		good_data=tmp_data(good_idx,lat_idx,time_idx) 
		data(time_idx,lat_idx)=total(good_data)/count
		good_data_nbr=good_data_nbr+count
	endif else begin
		data(time_idx,lat_idx)=1.0e36
	endelse
endfor; end loop over latitude
endfor; end loop over time
;
; We want to plot the Northernmost latitude at the top
;
if lat_min_idx gt lat_max_idx then begin
	data=reverse(data,2)
	lat=reverse(lat,1)
	swap=lat_min_idx
	lat_min_idx=lat_max_idx
	lat_max_idx=swap
endif; endif reversing lat order
;
;data=rotate(data,2)
;ord=rotate(ord,2)
;abc=rotate(abc,2)
;
scale_data
good_data=where(data lt very_big_nbr,count)
if count ne 0 then good_data=data(good_data) else good_data=data
data_max=max(good_data)
bin_max=data_max
data_min=min(data)
bin_min=data_min
;
cbar_fmt=replicate('UNSET',1)
if not usr_cntr_jd(fgr_idx) then cntr_lvl=0.
form_cntr_lvl
comp_cntr_lvl
cntr_lvl_min=min(cntr_lvl)
cntr_lvl_max=max(cntr_lvl)
;
;bld_sngs
cbar_setup
;
cbar_psn=[ $
	.88, $ ; x_min
	.10, $ ; y_min
	.92, $ ; x_max
	.90] ; 7_max
;
plt_rgn_nrm=[ $
	.10, $ ; x_min
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

plt_rgn_dvc=[ $
	plt_rgn_dvc(0), $
	plt_rgn_dvc(1), $
	plt_rgn_dvc(3), $
	plt_rgn_dvc(4)]

; Following must probably be used for non-cylindrical maps.
; using map_image has the effect of interpolating between data points 
; which results in a smoother image (and larger file) at printout
; data_img=map_image(data,startx,starty,scale=0.02)
; but this works best for cylindrical maps.
; using the raw data produces a grainier image and smaller file.
data_img=data

arr_sz=size(data_img)
x_arr_sz=arr_sz(1)
y_arr_sz=arr_sz(2)

; Convert x-y normal coordinate vector to device coords
img_psn_dvc=fix(plt_rgn_nrm*[ $
	!d.x_size, $
	!d.y_size, $
	!d.x_size, $
	!d.y_size])

; Get width and height in device units
dvc_wdt=img_psn_dvc(2)-img_psn_dvc(0)-1
dvc_hgt=img_psn_dvc(3)-img_psn_dvc(1)-1

;print,'lat_ctr = ',lat_ctr
;print,'lon_ctr = ',lon_ctr
;print,'map_lmt_dgr = ',map_lmt_dgr
;print,'img_psn_nrm = ',img_psn_nrm
;print,'plt_rgn_nrm = ',plt_rgn_nrm
;print,'lon_lat_wdw_x_sz = ',lon_lat_wdw_x_sz
;print,'lon_lat_wdw_y_sz = ',lon_lat_wdw_y_sz
;print,'img_psn_dvc = ',img_psn_dvc
;print,'dvc_wdt = ',dvc_wdt
;print,'dvc_hgt = ',dvc_hgt
;
; When scaling images it is important to place the top argument at
; no greater than !d.table_size-1 since there are at most 256 colors
; with PostScript and !d.table_size=256. Since the minimum value is
; always 0, then the colors that will show on the screen/paper are
; color indices 0..255.
if(!d.name eq 'PS') then begin
; Postscript will automatically scale the image up to xsize X ysize.
	tv, $  
	        bytscl( $
	                data_img, $
	                min=cntr_lvl_min, $
	                max=cntr_lvl_max, $
	                top=(!d.table_size-3))+2, $
	        img_psn_dvc(0)+1, $
	        img_psn_dvc(1)+1, $
		xsize=dvc_wdt, $
		ysize=dvc_hgt, $
		/device
endif else begin
; Other devices (pixel based systems) need to have the scaling done 
; for them in order to fit the device. 
	tv, $
		bytscl( $
			congrid( $
				data_img, $
				dvc_wdt, $
				dvc_hgt), $
	                min=cntr_lvl_min, $
	                max=cntr_lvl_max, $
	                top=(!d.table_size-3) $
			)+2, $
	        img_psn_dvc(0)+1, $
	        img_psn_dvc(1)+1, $
		/device
endelse

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

drw_ttl,ttl_jd(fgr_idx),''

print,'ord_max = ',ord_max
print,'ord_min = ',ord_min
print,'ord = ',ord
;
print,'Before'
print,'!x.crange(0) = ',!x.crange(0)
print,'!x.crange(1) = ',!x.crange(1)
print,'!y.crange(0) = ',!y.crange(0)
print,'!y.crange(1) = ',!y.crange(1)
;
plot, $
	abc, $	
	ord, $
	position=plt_rgn_nrm, $
	yrange=[ord_min,ord_max], $
	xstyle=13, $
	ystyle=13, $
	/noerase, $
	/nodata
;
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
	position=plt_rgn_nrm, $
	xstyle=13, $
	ystyle=13, $
	yrange=[ord_min,ord_max], $
	max_value=very_big_nbr, $
;	/closed, $				
	/noerase, $
	/overplot

print,'After'
print,'!x.crange(0) = ',!x.crange(0)
print,'!x.crange(1) = ',!x.crange(1)
print,'!y.crange(0) = ',!y.crange(0)
print,'!y.crange(1) = ',!y.crange(1)

lat_axz_drw, $
	lat_top=ord_max, $
	lat_btm=ord_min, $
	axz_vrt=1

time_axz_drw, $
	time_min=1, $
	time_max=n_elements(abc), $
	lbl_sty='none', $
	axz_vrt=0

time_min_sng=auto_sng(round(time_min),0)
time_max_sng=auto_sng(round(time_max),0)
time_sng='!5'+time_min_sng+' < !8t!5 < '+ $
		time_max_sng+' '

if (info_jd(fgr_idx) gt 0) then begin
	xyouts,0.99,.025,sub_ttl_sng+' '+sym_sng+' in '+rgn_lst(rgn_idx).string+': '+lon_sng+', '+lat_sng+', '+time_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
;	xyouts,1.0,.025,'Used '+auto_sng(good_data_nbr,0)+'/'+auto_sng(data_nbr,0)+' points.',alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
endif; endif info

print,'DEBUG: assigning data to tmp_data on exit from ibp_time_lat()'
data=tmp_data

end_of_procedure: foo=1
end; end time_lat()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure Time Latitude
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
