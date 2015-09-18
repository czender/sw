; $Id$

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin znl_avg() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro znl_avg_event,event

@ibp_cmn.com

event_type=tag_names(event,/structure)
;print,event_type,event

; In case event is not handled by the below
command_sng=''

; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value

if (usr_value eq 'Wznl_avg_gph') then begin
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
		inverse_znl_avg
	endif
endif; endif Wznl_avg_gph

if command_sng ne '' then print,command_sng

end; end znl_avg() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End znl_avg() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin make_znl_avg() 
; If the current device is X then this routine will create the 
; widgets and install the event handler for the Zonal Average graphs if
; they don't already exist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_znl_avg

@ibp_cmn.com

if !d.name eq 'X' then begin

; If the window is already up then OK else start it up.
; be careful because old widget ID's stay in common memory 
; after crashes!
widget_exists=1
if n_elements(Wznl_avg_gph) ne 0 then begin
	widget_exists=widget_info(Wznl_avg_gph,/valid_id)
endif

if ((widget_exists) and (n_elements(Wznl_avg_gph) ne 0)) then begin
curr_sz=indgen(2)
widget_control,Wznl_avg_base,tlb_get_size=curr_sz

; If the window has been resized, then destroy the widget hierarchy
; so that the widget can be recreated at the new size.
wset,znl_avg_wdw

if ((znl_avg_wdw_x_sz ne curr_sz(0)) or (znl_avg_wdw_y_sz ne curr_sz(1))) then begin
	znl_avg_wdw_x_sz=curr_sz(0)
	znl_avg_wdw_y_sz=curr_sz(1)
	widget_control,Wznl_avg_base,/destroy
	widget_exists=0
endif

endif

if ((widget_exists ne 1) or (n_elements(Wznl_avg_gph) eq 0)) then begin
; Create the zonal average widget hierarchy
Wznl_avg_base=widget_base( $
	group_leader=Wtop_lvl, $
	title='Zonal Average')
Wznl_avg_gph=widget_draw( $
	Wznl_avg_base, $
		xsize=znl_avg_wdw_x_sz, $
		ysize=znl_avg_wdw_y_sz, $
	/button_events, $
	uvalue='Wznl_avg_gph')

widget_control,Wznl_avg_base,/realize
widget_control,get_value=znl_avg_wdw,Wznl_avg_gph

; Call the event manager for the new base widgets
xmanager, $
	'make_znl_avg', $
	Wznl_avg_base, $
	group_leader=Wtop_lvl, $
	event_handler='znl_avg_event'		

endif; end creating widgets

widget_control,Wznl_avg_base,iconify=0 ; de-iconify when drawing
wset,znl_avg_wdw

endif; endif device is X

end; end make_znl_avg() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End make_znl_avg() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure Zonal Average
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro znl_avg
@ibp_cmn.com

erase

if data_nbr eq sv_data_nbr then begin
	print,'Zonally averaging two fields'
	TWO_FLD_LST=True
endif else begin
	print,'Zonally averaging one field'
	TWO_FLD_LST=False
endelse

; I'll have to do this inside a latitude loop here because the results
; of the where() function are indices into a 1D array, i.e., individual
; latitude information is lost. also must be careful when screening out 
; masked values that i don't graph any latitude circles with less than
; a minimum number of valid points. if all the data is in the valid 
; range then the simple statement z_avg=total(data,1)/lon_nbr will work.

z_avg=fltarr(lat_nbr)
sv_z_avg=fltarr(lat_nbr)
good_data_nbr=0.
sv_good_data_nbr=0.

for lat_idx=0,lat_nbr-1 do begin
	good_idx=where(data(*,lat_idx) lt very_big_nbr,count)
	vld_cnt_pct=100.0*count/lon_nbr
	if vld_cnt_pct lt vld_avg_thr_pct then begin
		print,'WARNING: Only '+auto_sng(vld_cnt_pct,1)+'% of lons valid in zonal average of '+fld_nm+' in lat('+auto_sng(lat_idx,0)+') = '+auto_sng(lat(lat_idx),1)+' degrees.'
	endif
	if count ne 0 then begin
		good_data=data(good_idx,lat_idx) 
		z_avg(lat_idx)=total(good_data)/count
		good_data_nbr=good_data_nbr+count
	endif else begin
		z_avg(lat_idx)=1.0e36
	endelse
endfor

if TWO_FLD_LST then begin
for lat_idx=0,lat_nbr-1 do begin
	good_idx=where(sv_data(*,lat_idx) lt very_big_nbr,count)
	vld_cnt_pct=100.*count/lon_nbr
	if vld_cnt_pct lt vld_avg_thr_pct then begin
		print,'WARNING: Only '+auto_sng(vld_cnt_pct,1)+'% of lons valid in zonal average of '+sv_fld_nm+' in lat('+auto_sng(lat_idx,0)+') = '+auto_sng(lat(lat_idx),1)+' degrees.'
	endif
	if count ne 0 then begin
		good_data=sv_data(good_idx,lat_idx) 
		sv_z_avg(lat_idx)=total(good_data)/count
		sv_good_data_nbr=sv_good_data_nbr+count
	endif else begin
		sv_z_avg(lat_idx)=1.0e36
	endelse
endfor
endif; endif TWO_FLD_LST

good_lat_idx=where(z_avg lt very_big_nbr,count)
good_lat_idx=where(sv_z_avg lt very_big_nbr,sv_count)
if (count eq 0) or (sv_count eq 0) then begin
	print,'There needs to be at least one latitude circle without a blocked point'
	goto,end_of_procedure
endif

;print,"fld_nm = ",fld_nm
;print,"data_nbr = ",data_nbr
;print,"good_data_nbr = ",good_data_nbr
;print,"z_avg = ",z_avg
;if TWO_FLD_LST then begin
;	print,"sv_fld_nm = ",sv_fld_nm
;	print,"sv_data_nbr = ",sv_data_nbr
;	print,"sv_good_data_nbr = ",sv_good_data_nbr
;	print,"sv_z_avg = ",sv_z_avg
;endif; endif TWO_FLD_LST

if TWO_FLD_LST then begin
	if sv_abb_sng eq abb_sng then begin
	if (strlen(ttl_jd(fgr_idx)) le 0) then ttl_jd(fgr_idx)=sv_abb_sng
	if (strlen(xttl_jd(fgr_idx)) le 0) then xttl_jd(fgr_idx)=slice_sng+' '+sym_sng+' ('+unit_sng+')'
	endif else begin
	if (strlen(ttl_jd(fgr_idx)) le 0) then ttl_jd(fgr_idx)=sv_abb_sng+' !8vs.!5 '+abb_sng
	if (strlen(xttl_jd(fgr_idx)) le 0) then xttl_jd(fgr_idx)=slice_sng+' '+sym_sng+', '+sv_slc_sng+' '+sv_sym_sng+' ('+unit_sng+')'
	endelse
	abc_min=min([min(z_avg),min(sv_z_avg)])
	abc_max=max( $
	[max(z_avg(where(z_avg lt very_big_nbr))), $
	max(sv_z_avg(where(sv_z_avg lt very_big_nbr)))])
endif else begin; endif there are two fields
	if (strlen(ttl_jd(fgr_idx)) le 0) then ttl_jd(fgr_idx)=abb_sng
	if (strlen(xttl_jd(fgr_idx)) le 0) then xttl_jd(fgr_idx)=slice_sng+' '+sym_sng+' ('+unit_sng+')'
	abc_min=min(z_avg)
	abc_max=max(z_avg(where(z_avg lt very_big_nbr)))
endelse; endif there is only one field

good_lat_idx=where(z_avg eq 0.0,count)
if count eq lat_nbr then begin
	print,'WARNING: Field seems to be identically zero.'
	abc_max=1.
endif
if (strlen(unit_jd(fgr_idx)) le 0) then unit_jd(fgr_idx)=unit_sng
if (strlen(yttl_jd(fgr_idx)) le 0) then yttl_jd(fgr_idx) ='!5Latitude (!5!E!12_!N!5)'

ord_min=min(lat)
ord_max=max(lat)
ord_min=(ord_min-5.0) + abs(ord_min mod 5.0)
ord_max=(ord_max+5.0) - abs(ord_max mod 5.0)

rng_x_min=abc_min
max_rng_x=abc_max
rng_y_min=ord_min
rng_y_max=ord_max
                 
; Set defaults for x,y ranges then override with any user specified values

if auto_scl_jd(fgr_idx) then begin
	rng_x_min_jd(fgr_idx)=rng_x_min
	max_rng_x_jd(fgr_idx)=max_rng_x
	rng_y_min_jd(fgr_idx)=rng_y_min
	rng_y_max_jd(fgr_idx)=max_rng_x
endif else begin
	if (rng_x_min_jd(fgr_idx) lt very_big_nbr) then rng_x_min=rng_x_min_jd(fgr_idx)
	if (max_rng_x_jd(fgr_idx) lt very_big_nbr) then max_rng_x=max_rng_x_jd(fgr_idx)
	if (rng_y_min_jd(fgr_idx) lt very_big_nbr) then rng_y_min=rng_y_min_jd(fgr_idx)
	if (rng_y_max_jd(fgr_idx) lt very_big_nbr) then rng_y_max=rng_y_max_jd(fgr_idx)
endelse

if TWO_FLD_LST then line_style=2 else line_style=0
if !d.name ne 'X' then chr_sz=2.0 else chr_sz=1.0

plot, $
	z_avg, $
	lat, $
	tit=ttl_jd(fgr_idx), $
	xtit=xttl_jd(fgr_idx), $
	ytit=yttl_jd(fgr_idx), $
	max_value=very_big_nbr, $
	xstyle=1, $
	ystyle=13, $
	xrange=[rng_x_min,max_rng_x], $
	yrange=[rng_y_min,rng_y_max], $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[3,3], $ ;[10,3] is [left,right] default
	ymargin=[3,2], $  ;[4,2] is [bottom,top] default
	linestyle=line_style

; Save the axes ranges for the next time through
rng_x_min_jd(fgr_idx)=!x.crange(0)
max_rng_x_jd(fgr_idx)=!x.crange(1)
rng_y_min_jd(fgr_idx)=!y.crange(0)
rng_y_max_jd(fgr_idx)=!y.crange(1)

lat_axz_drw, $
	lat_top=rng_y_max, $
	lat_btm=rng_y_min, $
	axz_vrt=1

if TWO_FLD_LST then begin
	oplot,	$
		sv_z_avg, $
		lat, $
		max_value=very_big_nbr, $
		thick=3.0, $
		linestyle=0
endif; endif there are two fields

if TWO_FLD_LST then begin
	ln_lgn_x1=.3
	ln_lgn_dx=0.07
	ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
	txt_lgn_x=ln_lgn_x2+0.01
	if !d.name ne 'X' then txt_lgn_sz=1.75 else txt_lgn_sz=1.0
	lgn_y_top=.7
	lgn_dy=0.1
	lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
		lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
	plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
	plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=2,thick=3.0,/NORMAL
;
	xyouts,txt_lgn_x,lgn_y(0),sv_src_sng+' '+sv_sym_sng,size=txt_lgn_sz,/NORMAL
	xyouts,txt_lgn_x,lgn_y(1),src_sng+' '+sym_sng,size=txt_lgn_sz,/NORMAL

endif; endif there are two fields

if (info_jd(fgr_idx) gt 0) then begin
if TWO_FLD_LST then begin
	xyouts,0.98,.025,sv_sub_ttl_sng+' '+sv_sym_sng+' !8vs.!5 '+sub_ttl_sng+' '+sym_sng+' in '+rgn_lst(rgn_idx).string+': '+lon_sng+', '+lat_sng,alignment=0.0,orientation=90.0,size=dbg_text_sz,/NORMAL
	info_sng='Used '+auto_sng(sv_good_data_nbr,0)+'/'+auto_sng(sv_data_nbr,0)+', '+auto_sng(good_data_nbr,0)+'/'+auto_sng(data_nbr,0)+' points.'
	info_sng=info_sng+' Data Avg. = '+auto_sng(sv_avg_data,1)+', '+auto_sng(avg_data,1)+'; Min. = '+auto_sng(sv_data_min,1)+', '+auto_sng(data_min,1)
	info_sng=info_sng+'; Max. = '+auto_sng(sv_data_max,1)+', '+auto_sng(data_max,1)+' '+unit_sng+', respectively.'
	xyouts,1.0,.025,info_sng,alignment=0.0,orientation=90.0,size=dbg_text_sz,/NORMAL
endif else begin; endif there are two fields
	xyouts,0.98,.025,sub_ttl_sng+' '+sym_sng+' in '+rgn_lst(rgn_idx).string+': '+lon_sng+', '+lat_sng,alignment=0.0,orientation=90.0,size=dbg_text_sz,/NORMAL
	info_sng='Used '+auto_sng(good_data_nbr,0)+'/'+auto_sng(data_nbr,0)+' points.'
	info_sng=info_sng+' Data Avg. = '+auto_sng(avg_data,1)+'; Min. = '+auto_sng(data_min,1)+'; Max. = '+auto_sng(data_max,1)+' '+unit_sng+'.'
	xyouts,1.0,.025,info_sng,alignment=0.0,orientation=90.0,size=dbg_text_sz,/NORMAL
endelse; endif there is only one field

sys_time=systime(0) 
xyouts,1.0,0.0,usr_sng+'@'+hostname_sng+' '+sys_time,size=dbg_text_sz,orientation=0.0,alignment=1.0,/NORMAL
endif; endif info

end_of_procedure: foo=1
end; end znl_avg()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure Zonal Average
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
