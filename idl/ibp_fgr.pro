; $Id$

; Purpose: Routines to manipulate global datasets

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin rfr_dim_lst()
; Strategy: Obtain and display the valid dimensions of the selected field.
; The routine also looks for default dimensions to select, i.e.,
; it selects lon for the x dimension and lat for the y dimension.
; No data is altered in this procedure.
; 
; Note that the longitude and latitude should always be re-read with the
; data so that they are both padded and regionalized in synchrony.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rfr_dim_lst
@ibp_cmn.com

print,'fld = '+auto_sng(fld,0)
print,'fld_lst(fld_idx).id = ',fld_lst(fld_idx).id
print,'fld_lst(fld_idx).name = ',fld_lst(fld_idx).name
print,'fld_lst(fld_idx).data_type = ',fld_lst(fld_idx).data_type
print,'fld_lst(fld_idx).dim_nbr = ',fld_lst(fld_idx).dim_nbr
for i=0,fld_lst(fld_idx).dim_nbr-1 do begin
	print,'fld_lst('+auto_sng(fld,0)+').dim('+auto_sng(i,0)+') = '+ $
	auto_sng(fld_lst(fld_idx).dim(i),0)+' = size of '+ $
	auto_sng(dim_lst(fld_lst(fld_idx).dim(i)).size,0)+' = '+ $
	dim_lst(fld_lst(fld_idx).dim(i)).name
endfor

dim_nbr=fld_lst(fld_idx).dim_nbr
dim_id_list=fld_lst(fld_idx).dim(0:dim_nbr-1)
dim_nm_list=dim_lst(dim_id_list).name
lev_dim_nm=''
for dim=0,dim_nbr-1 do begin
	dim_id=fld_lst(fld_idx).dim(dim)
	dim_nm=dim_lst(dim_id).name
	result=strpos(dim_nm,'lon')
	if result ne -1 then lon_dim_id=dim_id
	result=strpos(dim_nm,'lat')
	if result ne -1 then lat_dim_id=dim_id
	result=strpos(dim_nm,'ilev')
	if result ne -1 then begin
		lev_dim_id=dim_id
		lev_dim_nm='ilev'
	endif else begin
		result=strpos(dim_nm,'lev')
		if result ne -1 then begin
			lev_dim_id=dim_id
			lev_dim_nm='lev'
		endif
	endelse
	result=strpos(dim_nm,'time')
	if result ne -1 then begin 
		time_dim_id=dim_id
;		print,'DBG: found time dim, defining time_dim_id = ',time_dim_id
	endif
endfor

if n_elements(lon_dim_id) eq 0 then begin
	x_dim=0 
	x_dim_nm_list_loc=0 
	lon_nbr=1
endif else begin
	x_dim=lon_dim_id
	x_dim_nm_list_loc=(where(dim_id_list eq x_dim))(0)
	lon_nbr=dim_lst(lon_dim_id).size
endelse 
if n_elements(lat_dim_id) eq 0 then begin 
	y_dim=1 
	y_dim_nm_list_loc=0 
	lat_nbr=1
endif else begin 
	y_dim=lat_dim_id 
	lat_nbr=dim_lst(lat_dim_id).size
	y_dim_nm_list_loc=(where(dim_id_list eq y_dim))(0)
endelse
if n_elements(lev_dim_id) eq 0 then begin 
	z_dim=1 
	lev_nbr=1
endif else begin 
	z_dim=lev_dim_id 
	lev_nbr=dim_lst(lev_dim_id).size
endelse
;print,"DBG: checking time_dim_id to define time_nbr: time_dim_id = ",time_dim_id
if n_elements(time_dim_id) eq 0 then begin 
	t_dim=1 
	time_nbr=1
endif else begin 
	t_dim=time_dim_id 
	time_nbr=dim_lst(time_dim_id).size
;	print,"DBG: time_nbr=",time_nbr
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin NetCDF commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nc_id=ncdf_open(fl_in)
; now get the valid dimensions for the chosen variable.
; We retrieve the horizontal dimensions every time we reload 
; a variable because the region may have been truncated in the
; meantime. Get them here also so that there are values to use
; at widget instantiation.

ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,'lon',lon

; Get the vertical coordinate arrays. 
; alev and blev arrays are always available in cond-produced files,
; but are named differently in hsum-produced files.
result=ncdf_varid(nc_id,lev_dim_nm)
if result ne -1 then begin
	ncdf_varget,nc_id,lev_dim_nm,lev ; # % NCDF_VARID: Variable not found  "".
endif else begin; end if lev is a netcdf variable
result=ncdf_varid(nc_id,'mlev') ; # % NCDF_VARID: Variable not found  "mlev".
if result ne -1 then ncdf_varget,nc_id,'mlev',lev else lev = 0.0
endelse

; Due to ECMWF idiosyncracies, ECMWF variables are defined on 
; interfaces, not levels, so swap the aint data into the lev 
; array for ECMWF.
result=strpos(fl_in,'ecmwf')
if result ne -1 then begin
result=ncdf_varid(nc_id,lev_dim_nm)
if result ne -1 then begin
result=strpos(fl_in,'omega')
if result eq -1 then begin
;	ncdf_attget,nc_id,'alev',alev,/GLOBAL
;	ncdf_attget,nc_id,'amid',amid,/GLOBAL
;	ncdf_attget,nc_id,'bmid',bmid,/GLOBAL
;	ncdf_attget,nc_id,'aint',aint,/GLOBAL
;	ncdf_attget,nc_id,'bint',bint,/GLOBAL
;	level_pressure=amid*1000+bmid*1013
;	interface_pressure=aint*1000+bint*1013
	var_id=ncdf_varid(nc_id,'hyam')
	ncdf_varget,nc_id,var_id,hyam
	var_id=ncdf_varid(nc_id,'hybm')
	ncdf_varget,nc_id,var_id,hybm
	var_id=ncdf_varid(nc_id,'hyai')
	ncdf_varget,nc_id,var_id,hyai
	var_id=ncdf_varid(nc_id,'hybi')
	ncdf_varget,nc_id,var_id,hybi
	level_pressure=hyam*1000+hybm*1013
	interface_pressure=hyai*1000+hybi*1013
	lev=interface_pressure(1:lev_nbr)
;for level=0,lev_nbr-1 do begin
;print,'level(',level,') = ',lev(level),', level_pressure = ', $
;	level_pressure(level),' mb, interface pressure = ', $
;	interface_pressure(level)
;endfor
endif
endif
endif

; Get time dimension
result=ncdf_varid(nc_id,'time')
;print,"DBG: checking for variable time, got result = ",result
if result ne -1 then begin
	dim_id=ncdf_dimid(nc_id,'time')
	ncdf_diminq,nc_id,dim_id,foo,time_nbr
;	print,"DBG: got time dimension, time_nbr= ",time_nbr
; BEE's ccm2nc program sometimes creates a time variable but
; gives it a dimension size of zero.
	if time_nbr eq 0 then begin
		time_nbr=1
		time=indgen(time_nbr)
	endif else begin
		ncdf_varget,nc_id,result,time
		time_nbr=n_elements(time)
;		print,"DBG: getting time dimension again, time_nbr= ",time_nbr
	endelse
endif else begin
	result=ncdf_dimid(nc_id,'time')
	if result ne -1 then ncdf_diminq,nc_id,result,foo,time_nbr else time_nbr=1
	if time_nbr eq 0 then time_nbr=1
	time=indgen(time_nbr)
;	print,"DBG: getting time dimension in else loop, time_nbr= ",time_nbr
endelse
ncdf_close,nc_id
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of NetCDF commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; If widgets do not exist, draw them, and if they do exist, destroy them and/or update the values
; NB: Old widget IDs stay in common memory after crashes!
;print,'dbg: before x_dim exists 0'

; NB: 1998/10/19 This is where IBP crashes when run from gss1
; Since these two widgets do not do anything useful, the workaround
; is to just comment them out and fix it later.
;if n_elements(Wx_dim) ne 0 then widget_xst=widget_info(Wx_dim,/valid_id) else widget_xst=0
;print,'dbg: before x_dim destroy 1'
;if widget_xst eq 1 then widget_control,Wx_dim,/destroy
;
;print,'dbg: before x_dim create 2'
;Wx_dim=widget_list( $
;		Wdim_base, $
;		value=dim_nm_list, $
;		ysize=min([5,dim_nbr]), $
;		resource_name='x_dim', $
;		uvalue='Wx_dim')
;print,'dbg: before x_dim set_list_select 3'
;widget_control,Wx_dim,set_list_select=x_dim_nm_list_loc
;
;print,'dbg: before y_dim exists 5'
;if n_elements(Wy_dim) ne 0 then widget_xst=widget_info(Wy_dim,/valid_id) else widget_xst=0
;print,'dbg: before y_dim destroy 6'
;if widget_xst eq 1 then widget_control,Wy_dim,/destroy
;print,'dbg: before y_dim create 7'
;Wy_dim=widget_list( $
;		Wdim_base, $
;		value=dim_nm_list, $
;		ysize=min([5,dim_nbr]), $
;		resource_name='y_dim', $
;		uvalue='Wy_dim')
;print,'dbg: before y_dim set_list_select 8'
;widget_control,Wy_dim,set_list_select=y_dim_nm_list_loc

if n_elements(Wlat_list) ne 0 then widget_xst=widget_info(Wlat_list,/valid_id) else widget_xst=0
if widget_xst eq 1 then begin
widget_control,Wlat_list,set_value=strtrim(string(lat),2)
endif else begin
Wlat_list=widget_list( $
		Wlat_list_base, $
		value=strtrim(string(lat),2), $	
		ysize=5, $
		resource_name='lat_list', $
		uvalue='Wlat_list')
endelse

if n_elements(Wlon_list) ne 0 then widget_xst=widget_info(Wlon_list,/valid_id) else widget_xst=0
if widget_xst eq 1 then begin
widget_control,Wlon_list,set_value=strtrim(string(lon),2)
endif else begin
Wlon_list=widget_list( $
		Wlon_list_base, $
		value=strtrim(string(lon),2), $	
		ysize=5, $
		resource_name='lon_list', $
		uvalue='Wlon_list')
endelse

if vrt_slc ge lev_nbr then vrt_slc=lev_nbr-1
vrt_slc_lbl=strarr(lev_nbr)
for level=0,lev_nbr-1 do begin
	vrt_slc_lbl(level)=auto_sng(round(lev(level)),0)+' mb'
endfor; env loop over levels

if time_slc ge time_nbr then time_slc=time_nbr-1
time_lbl=strarr(time_nbr)
for time_idx=0,time_nbr-1 do begin
	time_lbl(time_idx)=auto_sng(time(time_idx),3)
endfor; end loop over times

if n_elements(Wvrt_slc) ne 0 then widget_xst=widget_info(Wvrt_slc,/valid_id) else widget_xst=0
if widget_xst eq 1 then widget_control,Wvrt_slc,/destroy
Wvrt_slc=widget_list( $
		Wvrt_slc_base, $
		value=vrt_slc_lbl, $
		ysize=min([5,lev_nbr]), $
		resource_name='vrt_slc', $
		uvalue='Wvrt_slc')
widget_control,Wvrt_slc,set_list_select=vrt_slc

if n_elements(Wtime_slc) ne 0 then widget_xst=widget_info(Wtime_slc,/valid_id) else widget_xst=0
if widget_xst eq 1 then widget_control,Wtime_slc,/destroy
;print,"DBG: about to set up time widget time_nbr=",time_nbr
Wtime_slc=widget_list( $
		Wtime_slc_base, $
		value=time_lbl, $
		ysize=min([5,time_nbr]), $
		resource_name='time_slc', $
		uvalue='Wtime_slc')
widget_control,Wtime_slc,set_list_select=time_slc

;print,'dbg: end of rfr_dim_lst()'
end; end rfr_dim_lst()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of rfr_dim_lst()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin read_fl()
; Strategy: First read in the available fields data and 
; Pass this information to the corresponding widget lists.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro read_fl
@ibp_cmn.com

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin NetCDF commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nc_id=ncdf_open(fl_in)
cdf_info={dim_nbr:0L,fld_nbr:0L,nbr_gatts:0L,rec_dim:0L}
foo=ncdf_inquire(nc_id)
cdf_info.dim_nbr=foo.ndims
cdf_info.fld_nbr=foo.nvars
cdf_info.nbr_gatts=foo.ngatts
cdf_info.rec_dim=foo.recdim
dim_nbr=cdf_info.dim_nbr
fld_nbr=cdf_info.fld_nbr

; Based on the number of dims and fields, allocate space for the structures
dim_sct={ $
	name:'', $
	id:0L, $
	size:0L}
dim_lst=replicate(dim_sct,dim_nbr)

fld_sct={ $
	name:'', $
	id:0L, $
	data_type:'', $
	dim_nbr:0L, $
	nbr_atts:0L, $
;	att_nm:strarr(10), $
;	att_val:strarr(10), $
	dim:lonarr(4)}

fld_lst=replicate(fld_sct,fld_nbr)
for dim_id=0,dim_nbr-1 do begin
	ncdf_diminq,nc_id,dim_id,name,size
	dim_lst(dim_id).name=name
	dim_lst(dim_id).id=dim_id
	dim_lst(dim_id).size=size
endfor

for dim_id=0,dim_nbr-1 do begin
	print,'dim('+auto_sng(dim_id,0)+') = '+dim_lst(dim_id).name+ $
	' = size of '+auto_sng(dim_lst(dim_id).size,0)
endfor

for fld_id=0,fld_nbr-1 do begin
	foo=ncdf_varinq(nc_id,fld_id)
	fld_lst(fld_id).name=foo.name
	fld_lst(fld_id).id=fld_id
	fld_lst(fld_id).data_type=foo.datatype
	fld_lst(fld_id).dim_nbr=foo.ndims
	fld_lst(fld_id).nbr_atts=foo.natts
	if (foo.ndims gt 0) then fld_lst(fld_id).dim(0:foo.ndims-1)=foo.dim
	dim_sng=''
	for idx=0,fld_lst(fld_id).dim_nbr-1 do begin
		dim_sng=dim_sng+auto_sng(fld_lst(fld_id).dim(idx),0)+' = '+dim_lst(fld_lst(fld_id).dim(idx)).name+', '
	endfor
	print,'fld('+auto_sng(fld_id,0)+') = '+fld_lst(fld_id).name+', '+fld_lst(fld_id).data_type+', ndims = '+auto_sng(fld_lst(fld_id).dim_nbr,0)+', '+dim_sng
;	for att_num=0,fld_lst(fld_id).nbr_atts-1 do begin $
;		fld_lst(fld_id).att_nm(att_num)=ncdf_attname(nc_id,fld_id,att_num)
;		ncdf_attget,nc_id,fld_id,fld_lst(fld_id).att_nm(att_num),value
;		fld_lst(fld_id).att_val(att_num)=string(value(*))
;	endfor
endfor
;
; say bye-bye
;
ncdf_close,nc_id
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of NetCDF commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Do not sort dimensions structure because that screws up just about everything
print,'Before weeding out extraneous fields and alphabetizing the rest, fld_nbr = ',auto_sng(fld_nbr,0)
swap=fld_lst
swap_good_idx=where(swap.dim_nbr gt 0,fld_nbr)
swap=swap(swap_good_idx)
swap=swap(sort(swap.name)) ; Alphabetize fld_lst
fld_lst=swap
swap=0.0
print,'After...fld_nbr = ',auto_sng(fld_nbr,0)

;cmd='touch ' + fl_in
;spawn,cmd

decrypt_fl_nm

end; end read_fl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of read_fl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin drw_pretty_picture()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro drw_pretty_picture

@ibp_cmn.com

; A call to clr_rfr() is wise before drawing any color plots
; because clr_rfr() will set the black and white colors to 
; the correct values according to the output device. Recall that
; X and postscript need to swap color 1 and color 0. The clr_mk()
; routine leaves color 0 and color 1 alone.

; Make sure the data in the current data array matches all the 
; selected dimensions since the rfr_dim_lst call, or else draw()
; will generate errors when attempting to draw data that hasn't 
; been read-in.

widget_control,/hourglass
clr_rfr,clr_tbl,tbl_fl,usr_dfn_clr_tbl

if !d.name ne 'X' then begin
	x_map_cbar_mrg=[7.5,3]	;[10,3] is [left,right] default
	y_map_cbar_mrg=[11,7]	;[4,2] is [bottom,top] default
endif else begin
	x_map_cbar_mrg=[7.5,3]	;[10,3] is [left,right] default
	y_map_cbar_mrg=[8,5]		;[4,2] is [bottom,top] default
endelse

; Figure-dependent procedures follow
if fgr_lst(fgr_idx).fnc_nm eq 'lon_lat_fll' then begin
	make_lon_lat
	cbar_setup
endif; endif lon_lat_fll
if fgr_lst(fgr_idx).fnc_nm eq 'lon_lat_img' then begin
	make_lon_lat
	cbar_setup
endif; endif lon_lat_img
if fgr_lst(fgr_idx).fnc_nm eq 'hst' then begin
	make_hst
endif; endif hst
if fgr_lst(fgr_idx).fnc_nm eq 'lat_lev' then begin
	make_lat_lev
endif; endif lat_lev
if fgr_lst(fgr_idx).fnc_nm eq 'lon_lev' then begin
	make_lon_lev
endif; endif lon_lev
if fgr_lst(fgr_idx).fnc_nm eq 'pol_img' then begin
	make_lon_lat
	cbar_setup
endif; endif pol_img
if fgr_lst(fgr_idx).fnc_nm eq 'time_lat' then begin
	make_time_lat
endif; endif time_lat
if fgr_lst(fgr_idx).fnc_nm eq 'scat' then begin
	make_scat
endif; endif scat
if fgr_lst(fgr_idx).fnc_nm eq 'znl_avg' then begin
	make_znl_avg
endif; endif znl_avg

call_procedure,fgr_lst(fgr_idx).fnc_nm
command_sng=fgr_lst(fgr_idx).english

end; end drw_pretty_picture()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of drw_pretty_picture()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin decrypt_fl_nm()
; Strategy: Much useful titling information can be stored in the
; filename. Here is where it is first translated.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro decrypt_fl_nm
@ibp_cmn.com

; Decipher the filename and decide what we're looking at
; NB: These source types MUST be known before the read_data()
; routine is called if there is to be any source-specific
; data extraction done. Since ERBE latitudes require special
; handling this is indeed currently the case.

src_nm=''
result=strpos(fl_in,'dstccm')
if result ne -1 then begin
	src_nm='dstccm'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'dstmch')
if result ne -1 then begin
	src_nm='dstmch'
	uscr_psn=strpos(fl_in,'_')
	if uscr_psn gt result then src_nm=strmid(fl_in,uscr_psn-8,8)
	goto,end_of_src_nm
endif
result=strpos(fl_in,'dstmzr')
if result ne -1 then begin
	src_nm='dstmzr'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'dst')
if result ne -1 then begin
	src_nm='dst'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'lsm')
if result ne -1 then begin
	src_nm='lsm'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'amip5')
if result ne -1 then begin
	src_nm='omega_amip'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'388')
if result ne -1 then begin
	src_nm='388'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'422')
if result ne -1 then begin
	src_nm='422'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'ecmwf')
if result ne -1 then begin
	src_nm='ecmwf'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'erbe')
if result ne -1 then begin
	src_nm='erbe'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'isccp')
if result ne -1 then begin
	src_nm='isccp'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'legates')
if result ne -1 then begin
	src_nm='legates'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'omega')
if result ne -1 then begin
	result=strpos(fl_in,'spcp')
	if result eq -1 then begin
		src_nm='omega'
		goto,end_of_src_nm
	endif
endif
result=strpos(fl_in,'plx22')
if result ne -1 then begin
	src_nm='plx22'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'spcp')
if result ne -1 then begin
	src_nm='spcp'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'ssmi')
if result ne -1 then begin
	src_nm='ssmi'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'tvbds')
if result ne -1 then begin
	src_nm='tvbds'
	goto,end_of_src_nm
endif
result=strpos(fl_in,'obs')
if result ne -1 then begin
	result=strpos(fl_in,'388')
	if result ne -1 then src_nm='388_obs'
	result=strpos(fl_in,'spcp')
	if result ne -1 then src_nm='spcp_obs'
	result=strpos(fl_in,'omega')
	if result ne -1 then src_nm='omega_obs'
endif 
end_of_src_nm: foo=1
;
end; end decrypt_fl_nm()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of decrypt_fl_nm()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin read_data()
; Strategy: Read in the field data and find out if it is pre-formatted.
; No data is altered in this procedure, just re-arranged zonally.
; Note that the longitude and latitude should always be re-read with the
; data so that they are both padded and regionalized in synchrony.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro read_data
@ibp_cmn.com
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin NetCDF commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nc_id=ncdf_open(fl_in)

; now get the data array
ncdf_varget,nc_id,fld_nm,data

; for now, also reload the horizontal dimensions in case they've
; been truncated. this is really a kludge, the permanent solution
; would be to always store original, untruncated versions of all
; the dimensions in a given file in memory at all times and swap 
; them into the user copy as necessary. then dimensions would only
; need to be read at the read_fl call, and they would only need
; to be redisplayed at the rfr_dim_lst call. 

ncdf_varget,nc_id,'lat',lat
ncdf_varget,nc_id,'lon',lon
lat_nbr=n_elements(lat)
lon_nbr=n_elements(lon)

; Get the vertical coordinate arrays. 
; alev and blev arrays are always available in cond-produced files,
; but are named differently in hsum-produced files.

result=ncdf_varid(nc_id,lev_dim_nm)
if result ne -1 then begin
	ncdf_varget,nc_id,lev_dim_nm,lev
endif else begin; end if lev is a netcdf variable
result=ncdf_varid(nc_id,'mlev')
if result ne -1 then ncdf_varget,nc_id,'mlev',lev else lev = 0.0
endelse

; Get Gaussian weights or other meridional weights
lat_wgt=(2.0/lat_nbr)+0.0*lat ; Default is all weights are equal and sum to 2.0
result=ncdf_varid(nc_id,'gw')
if result ne -1 then ncdf_varget,nc_id,'gw',lat_wgt else print,'WARNING: latitude dimension does not have gw variable to use for weights, using equal weights instead' 

; Due to ECMWF idiosyncracies, ECMWF variables are defined on 
; interfaces, not levels, so swap the aint data into the lev 
; array for ECMWF.

result=strpos(fl_in,'ecmwf')
if result ne -1 then begin
result=ncdf_varid(nc_id,lev_dim_nm)
if result ne -1 then begin
result=strpos(fl_in,'omega')
if result eq -1 then begin
	var_id=ncdf_varid(nc_id,'hyam')
	ncdf_varget,nc_id,var_id,hyam
	var_id=ncdf_varid(nc_id,'hybm')
	ncdf_varget,nc_id,var_id,hybm
	var_id=ncdf_varid(nc_id,'hyai')
	ncdf_varget,nc_id,var_id,hyai
	var_id=ncdf_varid(nc_id,'hybi')
	ncdf_varget,nc_id,var_id,hybi
	level_pressure=hyam*1000+hybm*1013
	interface_pressure=hyai*1000+hybi*1013
	lev=interface_pressure(1:lev_nbr)
;for level=0,lev_nbr-1 do begin
;print,'level(',level,') = ',lev(level),', level_pressure = ', $
;	level_pressure(level),' mb, interface pressure = ', $
;	interface_pressure(level)
;endfor
endif
endif
endif

; say bye-bye
ncdf_close,nc_id
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of NetCDF commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin data reformatting
; Put data into IDL compatible longitude order.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

data_nbr=n_elements(data)
dim_nbr=(size(data))(0)
print,'size(data) = ',size(data)

; For multi-level fields extract the slice which is of interest
; CAUTION: remember vrt_slc is offset by one to account for C
; indexing of IDL data arrays relative to FORTRAN numbering of levels in
; ECMWF and CCM, e.g., vrt_slc=0 corresponds to the bottom level.
; if there are more than two dimensions to the data then print out the
; pressure coordinates

if dim_nbr lt 2 then begin
	print,'Not enough dimensions in data. Bailing'
	return
endif

if dim_nbr eq 2 then begin
	slice_sng=''
	print,fld_nm,' is two-dimensional'
endif

if dim_nbr eq 3 then begin

; Figure out if data dimensions are lon,lat,time or lon,lev,lat or lon,lat,lev

if dim_lst(fld_lst(fld_idx).dim(2)).name eq 'time' then begin
	print,'Extracting time ',time_slc,' which is at ',time(time_slc)
	data=data(*,*,time_slc)
	data=reform(data)
 	slice_sng=''
endif else if dim_lst(fld_lst(fld_idx).dim(1)).name eq 'lev' then begin
	print,'Extracting level ',vrt_slc,' which is at ',lev(vrt_slc),' mb'
	data=data(*,vrt_slc,*)
	data=reform(data)
 	slice_sng='!5'+auto_sng(round(lev(vrt_slc)))+' mb'
endif else if dim_lst(fld_lst(fld_idx).dim(2)).name eq 'lev' then begin
	print,'Extracting level ',vrt_slc,' which is at ',lev(vrt_slc),' mb'
	data=data(*,*,vrt_slc)
	data=reform(data)
 	slice_sng='!5'+auto_sng(round(lev(vrt_slc)))+' mb'
endif else begin
	print,'Unparseable 3-D field in read_data()'
endelse
endif; endif dim_nbr eq 3

if dim_nbr eq 4 then begin
; Find level and time indices
	lev_idx=-1
	for dim_idx=0,dim_nbr-1 do begin
		if dim_lst(fld_lst(fld_idx).dim(dim_idx)).name eq 'time' then time_idx=dim_idx
		if dim_lst(fld_lst(fld_idx).dim(dim_idx)).name eq 'lev' then lev_idx=dim_idx
	endfor
	if time_idx eq -1 then print,'ERROR: dimension time not found'
	if lev_idx eq -1 then print,'ERROR: dimension lev not found'
	print,'Level dimension is index = ',lev_idx
	print,'Extracting level ',vrt_slc,' which is at ',lev(vrt_slc)
	if lev_idx eq 1 then begin
		data=data(*,vrt_slc,*,*)
	endif else if lev_idx eq 2 then begin
		data=data(*,*,vrt_slc,*)
	endif else begin
		print,'Unprogrammed level slice'
	endelse ; endif lev_idx=1
	data=reform(data)
	slice_sng='!5'+auto_sng(round(lev(vrt_slc)))+' mb'
	print,'Time dimension is index = ',time_idx
	print,'Extracting time ',time_slc,' which is at ',time(time_slc)
	if time_idx eq 0 then begin
		data=data(time_slc,*,*)
	endif else if time_idx eq 1 then begin
		data=data(*,time_slc,*)
	endif else if time_idx eq 2 then begin
		data=data(*,*,time_slc)
	endif else if time_idx eq 3 then begin
		data=data(*,*,time_slc)
	endif else begin
		print,'Unprogrammed time slice'
	endelse ; endif time_idx
	data=reform(data)
endif; endif 4 dimensions

; Fix CCM data to be IDL map-compatible with -180.0 < longitude < 180.0
if dim_lst(x_dim).name eq 'lon' then begin
if lon(lon_nbr-1) gt 90.0 and lon(0) ge 0.0 then begin
	if lon(0) eq 0.0 then nbr_shift=lon_nbr/2-1 else nbr_shift=lon_nbr/2
	lon=shift(lon,nbr_shift)
	west_lon_idx=where(lon gt 180.0,count)
	if count gt 0 then lon(west_lon_idx)=-360.0+lon(west_lon_idx)
	data=shift(data,nbr_shift,0)
	widget_control,Wlon_list,set_value=strtrim(string(lon),2)
endif; endif lon seems to be in the right range
endif; endif x_dim eq 'lon'

; ERBE files used to have latitude stored in descending order from highest (90N)
; to lowest (90S). This creates problems in the truncation routines. 
; Weirdly enough, data seem to be stored in correct order, so there's no 
; need to reverse it. Here's the kludge:

if dim_lst(y_dim).name eq 'lat' then begin
if ((lat(0) gt lat(1)) and (src_nm eq 'erbe')) then begin
	lat=reverse(lat)
	widget_control,Wlat_list,set_value=strtrim(string(lat),2)
endif; endif lat(0) > lat(1)
endif; endif y_dim eq 'lat'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End data reformatting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Find if selected field occurs in database of pre-formatted fields
; then take field and unit strings from there, otherwise use default.
pre_fmt=False
pre_fld_idx=0
while ((not pre_fmt) and (pre_fld_idx lt pre_fld_nbr)) do begin
	if (fld_nm eq pre_fld_lst(pre_fld_idx).name) then begin
		pre_fmt=True
		pre_fld=pre_fld_idx
	endif
	pre_fld_idx=pre_fld_idx+1
endwhile
if pre_fmt then print,'Acquired '+pre_fld_lst(pre_fld).english $
else print,'Acquired ',fld_nm

end; end read_data()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of read_data()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin bld_sngs()
; Must have the vertical level data to accomplish this if the field
; is multi-level
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro bld_sngs
@ibp_cmn.com

; Set defaults
year_sng=''
month_sng=''
src_sng=src_nm 

case src_nm of
'388': src_sng='!5CCM2'
'422': src_sng='!5CCM2 AMIP'
'dstccm': src_sng='!5CCM3'
'dstmch': src_sng='!5MATCH'
'dstmzr': src_sng='!5MOZART'
'ecmwf': src_sng='!5ECMWF'
'erbe': src_sng='!5ERBE'
'isccp': src_sng='!5ISCCP'
'legates': src_sng='!5Legates'
'lsm': src_sng='!5LSM'
'omega': src_sng='!7X!5!I5!N'
'omega_amip': begin
	src_sng='!7X!5!I5!N AMIP'
	result=strpos(fl_in,'omega')
	if result ne -1 then src_sng=src_sng+' - !7X!5!I5!N'
	end
'plx22': src_sng='!5CCM2'
'spcp': begin
	src_sng='!5SPCP'
	result=strpos(fl_in,'spcp_')
	if result ne -1 then src_sng='!5SPCP '+strmid(fl_in,result+5,2)
	result=strpos(fl_in,'omega')
	if result ne -1 then src_sng=src_sng+' - !7X!5!I5!N'
	end
'ssmi': src_sng='!5SSM/I'
'tvbds': src_sng='!5CCM2'
'388_obs': src_sng='!5CCM2'
'spcp_obs': src_sng='!5SPCP'
'omega_obs': src_sng='!5CCM2 !7X!5!I5!N'
else: 
endcase
;
result=strpos(fl_in,'_85_')
if result ne -1 then year_sng='!51985'
result=strpos(fl_in,'_86_')
if result ne -1 then year_sng='!51986'
result=strpos(fl_in,'_87_')
if result ne -1 then year_sng='!51987'
result=strpos(fl_in,'_88_')
if result ne -1 then year_sng='!51988'
result=strpos(fl_in,'_89_')
if result ne -1 then year_sng='!51989'
result=strpos(fl_in,'_8588_')
if result ne -1 then year_sng='!51985-1988'
result=strpos(fl_in,'_8589_')
if result ne -1 then year_sng='!51985-1989'
result=strpos(fl_in,'_8688_')
if result ne -1 then year_sng='!51986-1988'
result=strpos(fl_in,'_0120_')
if result ne -1 then year_sng='!520 Yr. Ensemble'
result=strpos(fl_in,'_0105_')
if result ne -1 then year_sng='!55 Yr. Ensemble'
result=strpos(fl_in,'_2080_')
if result ne -1 then year_sng='!51920-1980'
result=strpos(fl_in,'_8994_')
if result ne -1 then year_sng='!51989-1994'
result=strpos(fl_in,'_9094_')
if result ne -1 then year_sng='!51990-1994'
result=strpos(fl_in,'_8387_')
if result ne -1 then year_sng='!51983-1987'
result=strpos(fl_in,'_8488_')
if result ne -1 then year_sng='!51984-1988'
result=strpos(fl_in,'_87m85_')
if result ne -1 then year_sng='!51987-1985'
result=strpos(fl_in,'_87m89_')
if result ne -1 then year_sng='!51987-1989'
result=strpos(fl_in,'_8589x87_')
if result ne -1 then year_sng='!585,86,88,89'

result=strpos(fl_in,'_01')
if result ne -1 then month_sng='!5Jan'
result=strpos(fl_in,'_02')
if result ne -1 then month_sng='!5Feb'
result=strpos(fl_in,'_03')
if result ne -1 then month_sng='!5Mar'
result=strpos(fl_in,'_04')
if result ne -1 then month_sng='!5Apr'
result=strpos(fl_in,'_05')
if result ne -1 then month_sng='!5May'
result=strpos(fl_in,'_06')
if result ne -1 then month_sng='!5Jun'
result=strpos(fl_in,'_07')
if result ne -1 then month_sng='!5Jul'
result=strpos(fl_in,'_08')
if result ne -1 then month_sng='!5Aug'
result=strpos(fl_in,'_09')
if result ne -1 then month_sng='!5Sep'
result=strpos(fl_in,'_10')
if result ne -1 then month_sng='!5Oct'
result=strpos(fl_in,'_11')
if result ne -1 then month_sng='!5Nov'
result=strpos(fl_in,'_12')
if result ne -1 then month_sng='!5Dec'

; This routine must take place after read_data() so that 
; pre_fmt is defined.
if pre_fmt then begin
	abb_sng=pre_fld_lst(pre_fld).abbrev
	english_sng=pre_fld_lst(pre_fld).english
	fld_sng=pre_fld_lst(pre_fld).string
	sym_sng=pre_fld_lst(pre_fld).symbol
	unit_sng=pre_fld_lst(pre_fld).unit

	result=strpos(src_nm,'obs')
	if result ne -1 then begin
		src_sng=src_sng+' - '+pre_fld_lst(pre_fld).diff
	endif

endif else begin
	abb_sng=fld_nm
	english_sng=''
	fld_sng=fld_nm
	sym_sng=fld_nm
	unit_sng=''
endelse

ttl_sng=slice_sng+' '+fld_sng
if month_sng ne '' then ttl_sng=ttl_sng+', '+month_sng
sub_ttl_sng=src_sng+' '+year_sng

lat_min_sng=auto_sng(lat_min,1)
lat_max_sng=auto_sng(lat_max,1)
lat_sng='!5'+lat_min_sng+'!E!12_!5!N < !7u!5 < '+ $
		lat_max_sng+'!E!12_!5!N'

lon_min_sng=auto_sng(lon_min,1)
lon_max_sng=auto_sng(lon_max,1)
lon_sng='!5'+lon_min_sng+'!E!12_!5!N < !7k!5 < '+ $
		lon_max_sng+'!E!12_!5!N'

; WARNING: i have no idea what's going wrong here
;print,'lev_min = ',lev_min
;print,'size(lev_min) = ',size(lev_min)
;if n_elements(lev_min) eq 0 then lev_min=0.0
;if n_elements(lev_max) eq 0 then lev_max=1000.0
;if (size(lev_min))(0) eq 7 then lev_min=0.0
;if (size(lev_max))(0) eq 7 then lev_max=1000.0
;lev_min_sng=auto_sng(round(lev_min),0)
;lev_max_sng=auto_sng(round(lev_max),0)
;lev_sng='!5'+lev_min_sng+' mb < !7g!5 < '+ $
;		lev_max_sng+' mb'

end; end bld_sngs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of bld_sngs()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin longitude padding
; Make an extra column of data so that the right most longitude = 
; the leftmost longitude.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro longitude_pad
@ibp_cmn.com
;
if (min(lon) eq -180.0 and max(lon) lt 180.0) then begin
	print,'Padding longitude...'
	lon_nbr=lon_nbr+1
	lon = [lon,180.0]
	data = [data,data(0,*)]
endif
data_nbr=n_elements(data)
;
end; end longitude_pad()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End longitude padding
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin diagnose_data()
; Compute statistics on the data: min, max, avg.
; No action in this routine will alter the data.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro diagnose_data
@ibp_cmn.com

; Find maxima after masked values have been removed
good_data=where(data lt very_big_nbr,count)
if count ne 0 then good_data=data(good_data) else good_data=data
data_max=max(good_data)
bin_max=data_max
data_min=min(data)
bin_min=data_min

; Find average
avg_data=area_avg( $
	data, $
	lat, $
	lon, $
	lat_nbr, $
	lon_nbr, $
	lat_max_idx, $ ; index into lat array
	lat_min_idx, $ ; index into lat array
	lat_wgt=lat_wgt, $
	msk_data=msk_data, $
	msk_mode=msk_mode, $
	msk_val=msk_val, $
	vld_avg_thr_pct=vld_avg_thr_pct, $
	very_big_nbr=very_big_nbr)

; Print diagnostics on color scaling

print,'maximum all data= ',max(data)
print,'maximum good data= ',data_max
print,'minimum good data= ',data_min
print,'average data= ',avg_data
print,'size(data)= ',size(data)
print,'number good data= ',count

widget_control,Wscl,set_value=scl
widget_control,Wdata_min,set_value=data_min
widget_control,Wdata_max,set_value=data_max
widget_control,Wavg_data,set_value=avg_data

end; end diagnose_data()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of diagnose_data()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin form_cntr_lvl()
; Set the contour levels
; This routine does not change any data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro form_cntr_lvl
@ibp_cmn.com

cbar_fmt=replicate('UNSET',1)
if not usr_cntr_jd(fgr_idx) then cntr_lvl=0
if not usr_cntr_jd(fgr_idx) then cntr_lvl_nbr=0

; Set defaults for contour levels then override with any user specified values
if (cntr_lvl_min_jd(fgr_idx) lt very_big_nbr) then cntr_lvl_min=cntr_lvl_min_jd(fgr_idx)
if (cntr_ntv_jd(fgr_idx) lt very_big_nbr) then cntr_ntv=cntr_ntv_jd(fgr_idx)
if (cntr_lvl_nbr_jd(fgr_idx) lt very_big_nbr) then cntr_lvl_nbr=cntr_lvl_nbr_jd(fgr_idx)

; Find if the selected field is occurs in the database of pre-formatted fields
; then take the field and unit string from there, otherwise use a default.
if pre_fmt and (not usr_cntr_jd(fgr_idx)) then begin
	if pre_fld_lst(pre_fld).cntr_lvl_min lt very_big_nbr then begin
	cntr_lvl_min=pre_fld_lst(pre_fld).cntr_lvl_min
	cntr_lvl_nbr=pre_fld_lst(pre_fld).cntr_lvl_nbr
	cntr_ntv=pre_fld_lst(pre_fld).cntr_ntv
	endif; end if contour defaults are available
endif; endif pre_fmt

if ((not pre_fmt) or (auto_scl_jd(fgr_idx)) or (cntr_lvl_nbr eq 0) and (not usr_cntr_jd(fgr_idx))) then begin
	cntr_lvl_nbr=20
;	grd_ntv=10.0
	rsn_prj=(data_max-data_min)/cntr_lvl_nbr
	if rsn_prj eq 0.0 then grd_ntv = 1.0 $
	else if rsn_prj le 0.1 then grd_ntv = 10.0^round(alog10(rsn_prj)) $
	else if rsn_prj le 1.0 then grd_ntv = 1.0 $
	else if rsn_prj le 5.0 then grd_ntv = 1.0 $
	else if rsn_prj le 10.0 then grd_ntv = 5.0 $
	else if rsn_prj le 100.0 then grd_ntv = 10.0 $
	else grd_ntv = 10.0^round(alog10(rsn_prj))

	cntr_lvl_max=(data_max+grd_ntv) - abs(data_max mod grd_ntv)
	if (data_min lt 0.0) then $
		cntr_lvl_min=(data_min-grd_ntv) + abs(data_min mod grd_ntv) $
	else $
		cntr_lvl_min=data_min - abs(data_min mod grd_ntv) 

	if (cntr_lvl_min lt 0.0 and data_min ge 0.0) then cntr_lvl_min = 0.0
	cntr_ntv=(cntr_lvl_max-cntr_lvl_min)/(cntr_lvl_nbr-1)
	bin_sz=cntr_ntv

endif ; end if performing auto-scaling

cntr_lvl_min_jd(fgr_idx)=cntr_lvl_min
cntr_ntv_jd(fgr_idx)=cntr_ntv
cntr_lvl_nbr_jd(fgr_idx)=cntr_lvl_nbr

; Update the contour widgets, using Wcntr_lvl_min to see if the widgets exist.
widget_xst=0
if n_elements(Wcntr_lvl_min) ne 0 then begin
	widget_xst=widget_info(Wcntr_lvl_min,/valid_id)
endif
if widget_xst then begin
	widget_control,Wcntr_lvl_min,set_value=cntr_lvl_min_jd(fgr_idx)
	widget_control,Wcntr_ntv,set_value=cntr_ntv_jd(fgr_idx)
	widget_control,Wcntr_lvl_nbr,set_value=cntr_lvl_nbr_jd(fgr_idx)
endif; end updating widgets

end; end form_cntr_lvl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End form_cntr_lvl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin comp_cntr_lvl()
; Set the contour levels.
; This routine does not change any data.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro comp_cntr_lvl
@ibp_cmn.com

if n_elements(cntr_lvl) le 1 then cntr_lvl=cntr_lvl_min+cntr_ntv*findgen(cntr_lvl_nbr)

cntr_lvl_min=min(cntr_lvl)
cntr_lvl_max=max(cntr_lvl)

;print,'minimum color/contour level = ',cntr_lvl_min
;print,'maximum color/contour level = ',cntr_lvl_max

cntr_lbl_mk,cntr_lvl,cntr_lvl_min,cntr_ntv,cntr_lvl_nbr,unit_pfx_sng, $
	cntr_fll_idx,cntr_lvl_lbl,cntr_ln_sty,cntr_thk,cntr_which_lbl

; If widget does not exist, do nothing, but if it does exist, update the values.
; NB: Old widget IDs stay in common memory after crashes!
cntr_ntt=cntr_lvl_lbl

if n_elements(Wcntr_lvl) ne 0 then begin
	widget_xst=widget_info(Wcntr_lvl,/valid_id)
	if widget_xst eq 1 then begin
		widget_control,Wcntr_lvl,/destroy
		Wcntr_lvl=cw_bselector( $
				Wcntr_lvl_base, $
				strtrim(string(cntr_lvl),2), $
				label_left='Conts', $
				/return_index, $
;				resource_name='cntr_lvl', $
				uvalue='Wcntr_lvl')
	endif
endif

end; end comp_cntr_lvl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End comp_cntr_lvl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin scale_data()
; Scales the data to pre_fmt unit can actually change the data.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro scale_data
@ibp_cmn.com

scl=1.0
good_data_idx=where(data lt very_big_nbr)
if pre_fmt then begin
	scl=pre_fld_lst(pre_fld).scale
	data(good_data_idx)=data(good_data_idx)*scl
endif; endif

end; end scale_data()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End scale_data()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin reset_fld()
; Resets elements of the common blocks to the default state
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro reset_fld
@ibp_cmn.com

avg_data=0.0
cbar_fmt=replicate('UNSET',1)
cntr_lvl=0
cntr_lvl_lbl='0'
data_max=0.0
data_min=0.0
pre_fmt=False 
scale_fct=1.0e0

end; end reset_fld()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End reset_fld()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin massage_data()
; Fiddle with the data to make it look sharp on both displays and 
; hard-copy.
; This routine can change the data and should be checked very carefully.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro massage_data
@ibp_cmn.com
;
if data_min lt cntr_lvl_min then begin 
	print,'WARNING: some good data lies outside color/contour range'
	print,'       : automatically filling that data with minimal color'
	fudge_data=where(data lt cntr_lvl_min,count)
	if data_min gt 0 then begin
		if count ne 0 then fudge_data=where(data(fudge_data) gt 0,count)
		if count ne 0 then data(fudge_data)=cntr_lvl_min+(cntr_lvl_max-cntr_lvl_min)/250.0
	endif
endif
;
;CAUTION!!! this part is intended to prevent using black for data too
;close to the lower boudary, it's an untested hack intended for horizontal
;slice plots where being off by 1 color level out of 256 is unnoticeable, 
;but it's not intended for use with histograms or other data analysis because
;it corrupts data in the first color range.
;print,'d.table_size = ',!d.table_size
;fudge_data=where(data lt cntr_lvl_min+(cntr_lvl_max-cntr_lvl_min)/250,count)
;if count ne 0 then data(fudge_data)=cntr_lvl_min+(cntr_lvl_max-cntr_lvl_min)/250.0
;
; Get rid of any masked out areas by setting them equal to white
if data_max gt very_big_nbr then begin
	masked_values=where(data gt very_big_nbr)
	if !d.name eq 'X' then $
		data(masked_values)=0 $
	else $
		data(masked_values)=very_big_nbr 
endif
;
if data_max gt cntr_lvl_max then begin 
	print,'WARNING: some data lies outside color/contour range'
	print,'       : automatically filling that data with maximal color'
	fudge_data=where(data gt cntr_lvl_max,count)
	if count ne 0 then fudge_data=where(data(fudge_data) lt very_big_nbr,count)
	if count ne 0 then data(fudge_data)=cntr_lvl_max-(cntr_lvl_max-cntr_lvl_min)/100.0
endif
;
print,'new minimum data= ',min(data)
print,'new maximum data= ',max(data)
;
widget_control,Wdata_min,set_value=data_min
widget_control,Wdata_max,set_value=data_max
;
end; end massage_data()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End massage_data()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Colorbar Setup
; This routine should go after data is in final state
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cbar_setup
@ibp_cmn.com

; Since colors 0,1, and !d.table_size-1 are usually filled in with black and white
; and that is what missing/blocked/largerthanmaxcontour/smallerthanmincontour
; values will show up as, we just use !d.table_size-3 colors on color bar.

if fgr_lst(fgr_idx).fnc_nm eq 'lon_lat_fll' then begin
	cbar_clr_nbr=cntr_lvl_nbr-1
	cbar_lgn=cntr_lvl_lbl
endif; endif contour map

if fgr_lst(fgr_idx).fnc_nm eq 'lat_lev' then begin
	cbar_clr_nbr=cntr_lvl_nbr-1
	cbar_lgn=cntr_lvl_lbl
endif; endif contour map

if fgr_lst(fgr_idx).fnc_nm eq 'lon_lev' then begin
	cbar_clr_nbr=cntr_lvl_nbr-1
	cbar_lgn=cntr_lvl_lbl
endif; endif contour map

if fgr_lst(fgr_idx).fnc_nm eq 'lon_lat_img' or $
	fgr_lst(fgr_idx).fnc_nm eq 'time_lat' or $ 
	fgr_lst(fgr_idx).fnc_nm eq 'pol_img' then begin
	cbar_clr_nbr=(!d.table_size-3)
	cbar_lgn=strarr(cbar_clr_nbr+1)

; find the color index associated with each contour level
	cbar_lgn(0)=cntr_lvl_lbl(0)
	for level=1,cntr_lvl_nbr-2 do begin
		color_idx=(!d.table_size-1)* $
			(cntr_lvl(level)-cntr_lvl_min)/ $
			(cntr_lvl_max-cntr_lvl_min)
		color_idx=round(color_idx)
		cbar_lgn(color_idx)=cntr_lvl_lbl(level)
	endfor
	cbar_lgn(cbar_clr_nbr)=cntr_lvl_lbl(cntr_lvl_nbr-1)

endif; endif image map

cbar_idx=indgen(cbar_clr_nbr)+2
cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
if !d.name ne 'X' then cbar_chr_sz=1.3 else cbar_chr_sz=1.0

if (strlen(unit_jd(fgr_idx)) le 0) then unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
cbar_unit=unit_jd(fgr_idx)
if !d.name ne 'X' then cbar_lbl_sz=1.5 else cbar_lbl_sz=1.0

end; end cbar_setup()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Colorbar Setup
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Region Truncation
; This routine truncates two dimensional data to the given bounds
; Data are rearranged in order of increasing longitude so that the leftmost column.,
; lon_min, must be greater than lon_max if the region crosses the dateline,
; but map_lmt_lon_min must always be numerically less than map_lmt_lon_max.
; diagnose_data() should be called after this routine to recompute the
; data min, max, avg, etc.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro rgn_trunc
@ibp_cmn.com

; Current method requires supplying desired lon_min and lon_max
; as values actually in longitude array. 
; These statements define ???_???_idx for ensuing operations.

;print,'lon = ',lon
;print,'abs(lon-lon_min) = ',abs(lon-lon_min)
foo=min(abs(lon-lon_min),lon_min_idx)
foo=min(abs(lon-lon_max),lon_max_idx)
foo=min(abs(lat-lat_min),lat_min_idx)
foo=min(abs(lat-lat_max),lat_max_idx)

print,'Beginning Region Truncation...'
print,'lon_nbr = ',lon_nbr,' min(lon) = ',min(lon),' max(lon) = ',max(lon)
print,'lon_min = ',lon_min,' lon_min_idx = ',lon_min_idx,' = ',lon(lon_min_idx)
print,'lon_max = ',lon_max,' lon_max_idx = ',lon_max_idx,' = ',lon(lon_max_idx)
print,''
print,'lat_nbr = ',lat_nbr,' min(lat) = ',min(lat),' max(lat) = ',max(lat)
print,'lat_min = ',lat_min,' lat_min_idx = ',lat_min_idx,' = ',lat(lat_min_idx)
print,'lat_max = ',lat_max,' lat_max_idx = ',lat_max_idx,' = ',lat(lat_max_idx)

; Is there to be any zonal regionalization?
if ((lon_min_idx ne 0) or (lon_max_idx ne lon_nbr-1)) then begin

; If zonal range crosses date line.....
if lon_min gt lon_max then begin
; Find current index of desired left-most longitude
	shift_amt=lon_nbr-lon_min_idx
	data=shift(data,shift_amt,0)
	lon=shift(lon,shift_amt)

; Truncate data to right of desired right-most longitude
	data=data(0:lon_max_idx+shift_amt,*)
	lon=lon(0:lon_max_idx+shift_amt)

	lon_ctr=180.0-0.5*(lon_min+lon_max)
	map_lmt_dgr(1)=lon_min
	map_lmt_dgr(3)=360.0+lon_max

endif; endif data region crossed date line

; If zonal range does not cross date line.....
if lon_min lt lon_max then begin
; Truncate data outside desired longitude region
	data=data(lon_min_idx:lon_max_idx,*)
	lon=lon(lon_min_idx:lon_max_idx)

	lon_ctr=0.5*(lon_min+lon_max)
	map_lmt_dgr(1)=lon_min
	map_lmt_dgr(3)=lon_max

endif; endif data region did not cross date line

endif; endif longitude needed regionalization

; Is there to be any meridional regionalization?
if ((lat_min_idx ne 0) or (lat_max_idx ne lat_nbr-1)) then begin
; Truncate data outside desired latitude region
	data=data(*,lat_min_idx:lat_max_idx)
	lat=lat(lat_min_idx:lat_max_idx)

	lat_ctr=0.5*(lat_min+lat_max)
	map_lmt_dgr(0)=lat_min
	map_lmt_dgr(2)=lat_max

endif; endif latitude needed truncating

lat_nbr=n_elements(lat)
lon_nbr=n_elements(lon)
data_nbr=n_elements(data)

widget_control,Wlat_list,set_value=strtrim(string(lat),2)
widget_control,Wlon_list,set_value=strtrim(string(lon),2)

good_idx=where(data lt very_big_nbr,good_data_nbr)
bad_idx=where(data gt very_big_nbr,bad_data_nbr)
print,'Post rgn_trunc....'
print,'lat_nbr = ',lat_nbr
print,'lon_nbr = ',lon_nbr
print,'data_nbr = ',data_nbr
if bad_data_nbr gt 0 then begin
	dbg_lvl=0
	print,'good_data_nbr = ',good_data_nbr
	print,'bad_data_nbr = ',bad_data_nbr
	if dbg_lvl gt 2 then begin
		for dat_idx=0,bad_data_nbr-1 do begin
			print,'data(',bad_idx(dat_idx),') = ',data(bad_idx(dat_idx))
		endfor ; end loop over dat_idx
	endif ; endif dbg_lvl gt 2
endif ; endif bad_data_nbr ne 0

; boundary arguments to map_set are a bit weird:
; When data crosses date line then lon_min and lon_max should
; be specified such that lon_min is numerically less than lon_max no
; matter what. This means it's o.k. to go outside range -180 < lon < 180
; for these arguments.
;map_lmt_dgr=[ $ ;[latmin,lonmin,latmax,lonmax] drawing box for map_set
;			lat_min, $
;			map_lmt_lon_min, $
;			lat_max, $
;			map_lmt_lon_max] 
;
end; end rgn_trunc()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Region Truncation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Zvector 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro zvector, $
	u_data, $
	v_data, $
	mag_data, $
	lon, $
	lat, $
	length=length, $
	color=color, $
	spacing=spacing
;
max_abs_u=max(abs(u_data))
max_abs_v=max(abs(v_data))
;
abs_mag_data=abs(mag_data)
max_abs_mag_data=max(abs_mag_data)
;
lon_nbr=n_elements(lon)
lat_nbr=n_elements(lat)
;
; this is true for even #'s of lats. and lons.
nbr_x_arrows=fix(lon_nbr/spacing)
nbr_y_arrows=fix(lat_nbr/spacing)
;
x_dgr_per_cell=(lon(spacing*nbr_x_arrows-1)-lon(0))/(nbr_x_arrows-1)
y_dgr_per_cell=(lat(spacing*nbr_y_arrows-1)-lat(0))/(nbr_y_arrows-1)
;
arrow_x_length=x_dgr_per_cell*length*u_data/max_abs_u
arrow_y_length=y_dgr_per_cell*length*v_data/max_abs_v
;this fixes the aspect ratio so that the vertical length of the arrow
;corresponds to the same wind speed as the same horizontal length, even
;if the lat:lon aspect ratio of the whole map is not 1:1
arrow_y_length=arrow_y_length*0.5
;
print,'lon(0), lon(lon_nbr-1) = ',lon(0), lon(lon_nbr-1)
print,'lat(0), lat(lat_nbr-1) = ',lat(0), lat(lat_nbr-1)
print,'x_dgr_per_cell = ',x_dgr_per_cell,' degrees lon.'
print,'y_dgr_per_cell = ',y_dgr_per_cell,' degrees lat.'
print,'max(arrow_x_length) = ',max(abs(arrow_x_length)),' degrees long.'
print,'max(arrow_y_length) = ',max(abs(arrow_y_length)),' degrees lat.' 
;
for lon_idx=0,lon_nbr-1,spacing do begin
	arrow_tail_x=lon(lon_idx)
	for lat_idx=0,lat_nbr-1,spacing do begin
		arrow_head_x=arrow_tail_x+arrow_x_length(lon_idx,lat_idx)
		arrow_tail_y=lat(lat_idx)
		arrow_head_y=arrow_tail_y+arrow_y_length(lon_idx,lat_idx)
		if (arrow_head_x gt -180.0) and (arrow_head_x lt 180.0) then $
		arrow, $
			arrow_tail_x, $
			arrow_tail_y, $
			arrow_head_x, $
			arrow_head_y, $
			color=color, $
;			hsize=!d.x_size / 64., $	;head size
			hsize=!d.x_size*abs_mag_data(lon_idx,lat_idx)/ $
				(100.0*max_abs_mag_data), $;head size
			hthick=1.0, $			;head thickness
			thick=1.0, $			;body thickness
;			/SOLID, $
			/DATA
	endfor						; end lat_idx
endfor							; end lon_idx

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Zvector 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Tick Format
; Purpose: Common axes and title drawing procedure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function tck_fmt_get,axis,index,value
@ibp_cmn.com

; Longitude axis tick values are sometimes returned relative to center longitude of map and so may need to be corrected
; 19991230: Is this statement true any longer?
if ((axis eq 0) and (lon_ctr ne 0) and (fgr_lst(fgr_idx).fnc_nm ne 'lat_lev') and fgr_lst(fgr_idx).fnc_nm ne 'lon_lev') then begin
; 19991230: Commented out following line, which gives erroneous results for regions not centered on Greenwich
;	value=value+lon_ctr
	if value gt 180 then value=value-360
endif

nint_value=round(value)
abs_nint_value=abs(nint_value)

if abs_nint_value ge 100 then fmt_sng='(I3)' else $
if abs_nint_value ge 10  then fmt_sng='(I2)' else $
fmt_sng='(I1)' 

if fgr_lst(fgr_idx).fnc_nm eq 'lat_lev' then begin	;latitude axis
	if nint_value eq 0 then return,'!50!E!12_!N!5'
	if nint_value gt 0 then return,'!5'+ $
		string(format=fmt_sng,abs_nint_value)+ $
		'!E!12_!N!5N'
	if nint_value lt 0 then return,'!5'+ $
		string(format=fmt_sng,abs_nint_value)+ $
		'!E!12_!N!5S'
endif

if fgr_lst(fgr_idx).fnc_nm eq 'lon_lev' then begin	;longitude axis
	if nint_value eq 0 then return,'!50!E!12_!N!5'
	if nint_value gt 0 then return,'!5'+ $
		string(format=fmt_sng,abs_nint_value)+ $
		'!E!12_!N!5E'
	if nint_value lt 0 then return,'!5'+ $
		string(format=fmt_sng,abs_nint_value)+ $
		'!E!12_!N!5W'
endif

if axis eq 0 then begin	;longitude axis
	if nint_value eq 0 then return,'!50!E!12_!N!5'
	if nint_value gt 0 then return,'!5'+ $
		string(format=fmt_sng,abs_nint_value)+ $
		'!E!12_!N!5E'
	if nint_value lt 0 then return,'!5'+ $
		string(format=fmt_sng,abs_nint_value)+ $
		'!E!12_!N!5W'
endif

if axis eq 1 then begin	;latitude axis
	if nint_value eq 0 then return,'!50!E!12_!N!5'
	if nint_value gt 0 then return,'!5'+ $
		string(format=fmt_sng,abs_nint_value)+ $
		'!E!12_!N!5N'
	if nint_value lt 0 then return,'!5'+ $
		string(format=fmt_sng,abs_nint_value)+ $
		'!E!12_!N!5S'
endif

return,tck_lbl
end; end tck_fmt_get()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Tick Format
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Draw Axes
; Define the axes drawing procedure for contour and image maps.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro drw_axz
@ibp_cmn.com

if !d.name ne 'X' then tick_charsize=1.5 else tick_charsize=1.0
tick_charthick=1.5

; tick marks on the y-axes tend to get stretched out so make them smaller
x_tck_lnggth=-0.02	;default is 0.02, negative means ticks point outwards
y_tck_lnggth=-0.01	;default is 0.02, negative means ticks point outwards
axz_thk=1.0

lat_top=map_lmt_dgr(2)
lat_btm=map_lmt_dgr(0)

lon_lft=map_lmt_dgr(1)
lon_rgt=map_lmt_dgr(3)

lon_spn=lon_rgt-lon_lft
lat_spn=lat_top-lat_btm

if lon_spn ge 150 then begin
	lon_lbl_int=30.0 
	lon_tck_mnr_nbr=3
endif else $
if lon_spn ge 100 then begin
	lon_lbl_int=20.0 
	lon_tck_mnr_nbr=4
endif else $
if lon_spn ge 30 then begin
	lon_lbl_int=10.0 
	lon_tck_mnr_nbr=5
endif else begin
	lon_lbl_int=5.
	lon_tck_mnr_nbr=5
endelse

if lat_spn ge 120 then begin
	lat_lbl_int=30.0
	lat_tck_mnr_nbr=3
endif else $
if lat_spn ge 90 then begin
	lat_lbl_int=20.0
	lat_tck_mnr_nbr=4
endif else $
if lat_spn ge 60 then begin
	lat_lbl_int=10.0
	lat_tck_mnr_nbr=5
endif else begin
	lat_lbl_int=5.
	lat_tck_mnr_nbr=5
endelse

lon_lbl_nbr=fix(lon_spn/lon_lbl_int)
lat_lbl_nbr=fix(lat_spn/lat_lbl_int)

axis, $
	xaxis=0, $
	xstyle=1, $
	xthick=axz_thk, $
	xtick_get=lon_tck, $
	xticks=lon_lbl_nbr, $
	xminor=lon_tck_mnr_nbr, $
	ticklen=x_tck_lnggth, $
	xtickformat='tck_fmt_get', $
	charthick=tick_charthick, $
	charsize=tick_charsize, $
	/nodata

axis, $
	yaxis=0, $
	ystyle=1, $
	ythick=axz_thk, $
	ytick_get=lat_tck, $
	yticks=lat_lbl_nbr, $
	yminor=lat_tck_mnr_nbr, $
	ticklen=y_tck_lnggth, $
	ytickformat='tck_fmt_get', $
	charthick=tick_charthick, $
	charsize=tick_charsize, $
	/nodata

; NB: using a null string in xtickname (in order to have blank labels)
; does not work because it is over-ridden by some internal labeling 
; routine.

axis, $
	xaxis=1, $
	xstyle=1, $
	xthick=axz_thk, $
	xtick_get=lon_tck, $
	xticks=lon_lbl_nbr, $
	xminor=lon_tck_mnr_nbr, $
	ticklen=x_tck_lnggth, $
	xtickformat='null_fmt', $
;	xtickname=strarr(lon_lbl_nbr), $
	charthick=tick_charthick, $
	charsize=tick_charsize, $
	/nodata

axis, $
	yaxis=1, $
	ystyle=1, $
	ythick=axz_thk, $
	ytick_get=lat_tck, $
	yticks=lat_lbl_nbr, $
	yminor=lat_tck_mnr_nbr, $
	ticklen=y_tck_lnggth, $
	ytickformat='null_fmt', $
;	ytickname=strarr(lat_lbl_nbr), $
	charthick=tick_charthick, $
	charsize=tick_charsize, $
	/nodata

;plots, $
;	[plt_rgn_nrm(0), $
;	plt_rgn_nrm(2), $
;	plt_rgn_nrm(2)], $
;	[plt_rgn_nrm(3), $
;	plt_rgn_nrm(3), $
;	plt_rgn_nrm(1)], $
;	thick=axz_thk, $
;	/normal, $
;	linestyle=0

dbg_lvl=0
if dbg_lvl gt 0 then begin
	print,'ibp_fgr.pro:drw_axz() reports '
	print,'lon_tck = ',lon_tck
	print,'lat_tck = ',lat_tck
	print,'lon_ctr = ',lon_ctr
	print,'lat_ctr = ',lat_ctr
endif; endif dbg

end; end drw_axz()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Draw Axes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Draw Title
; Define the title drawing procedure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro drw_ttl,ttl_sng,sub_ttl_sng

if !d.name ne 'X' then ttl_sz=2.0 else ttl_sz=1.25
if !d.name ne 'X' then sub_ttl_sz=1.5 else sub_ttl_sz=0.75

xyouts,0.5,0.95,ttl_sng,size=ttl_sz,alignment=0.5,/NORMAL
xyouts,0.5,0.90,sub_ttl_sng,size=sub_ttl_sz,alignment=0.5,/NORMAL

end; end drw_ttl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Draw Title
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
