; $Id$

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin make_cstm() 
; If the current device is X then this routine will create the 
; widgets and install the event handler for the XY Range panel if
; they don't already exist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_cstm

@ibp_cmn.com

if !d.name eq 'X' then begin

; If the window is already up then OK else start it up.
; be careful because old widget ID's stay in common memory 
; after crashes!

widget_xst=1
if n_elements(Wcstm_base) ne 0 then begin
	widget_xst=widget_info(Wcstm_base,/valid_id)
endif

if ((widget_xst ne 1) or (n_elements(Wcstm_base) eq 0)) then begin

; Create the custom widget hierarchy
case fgr_lst(fgr_idx).fnc_nm of
	'lon_lat_fll': begin
	end; end lon_lat_fll
	'lon_lat_img': begin
	end; end lon_lat_img
	'hst': begin
	end; end hst
        'lat_lev': begin
	end; end lat_lev
        'lon_lev': begin
	end; end lon_lev
	'pol_img': begin
	end; end lon_lat_img
        'time_lat': begin
	end; end time_lat
	'scat': begin
	end; end scat
	'znl_avg': begin
	end; end znl_avg
	else: begin
	print,'Figure not defined'
	end; end else
endcase

cstm_ttl='Customize '+fgr_lst(fgr_idx).english

Wcstm_base=widget_base( $
	group_leader=Wtop_lvl, $
	title=cstm_ttl, $
	resource_name='ibp', $
	/column)

Wcstm_top=widget_base($
	Wcstm_base, $
	/column)

Wcstm_first_pnl=widget_base($
	Wcstm_base, $
	/row)

Wcstm_second_pnl=widget_base($
	Wcstm_base, $
	/row)

Wcstm_third_pnl=widget_base($
	Wcstm_base, $
	/row)

Wcstm_btm=widget_base($
	Wcstm_base, $
	/row)

Wttl=cw_field( $
	Wcstm_top, $
	title='       title:', $
	value=ttl_jd(fgr_idx), $
	xsize= 40, $
	/all_events, $
	uvalue='Wttl', $
	/row)

Wx_ttl=cw_field( $
	Wcstm_top, $
	title='      xtitle:', $
	value=xttl_jd(fgr_idx), $
	xsize= 40, $
	/all_events, $
	uvalue='Wx_ttl', $
	/row)

Wy_ttl=cw_field( $
	Wcstm_top, $
	title='      ytitle:', $
	value=yttl_jd(fgr_idx), $
	xsize= 40, $
	/all_events, $
	uvalue='Wy_ttl', $
	/row)

Wunit=cw_field( $
	Wcstm_top, $
	title='       unit:', $
	value=unit_jd(fgr_idx), $
	xsize= 40, $
	/all_events, $
	uvalue='Wunit', $
	/row)

Wmisc=cw_field( $
	Wcstm_top, $
	title='        misc:', $
	value=misc_jd(fgr_idx), $
	xsize= 40, $
	/all_events, $
	uvalue='Wmisc', $
	/row)

Winfo_base=widget_base($
	Wcstm_first_pnl, $
	/exclusive, $
	/row)

Winfo=widget_button( $
	Winfo_base, $
	value='Info', $
	uvalue='Winfo')

Wauto_scl_base=widget_base($
	Wcstm_first_pnl, $
	/exclusive, $
	/row)

Wauto_scl=widget_button( $
	Wauto_scl_base, $
	value='Auto_Scl', $
	uvalue='Wauto_scl')

Wrng_x_min=cw_field( $
	Wcstm_third_pnl, $
	title='X Min ', $
	value=rng_x_min_jd(fgr_idx), $
	xsize=11, $
	/all_events, $
	/row, $
	/float, $
	uvalue='Wrng_x_min')

Wmax_rng_x=cw_field( $
	Wcstm_third_pnl, $
	title='X Max ', $
	value=max_rng_x_jd(fgr_idx), $
	xsize=11, $
	/all_events, $
	/row, $
	/float, $
	uvalue='Wmax_rng_x')

Wrng_y_min=cw_field( $
	Wcstm_third_pnl, $
	title='Y Min ', $
	value=rng_y_min_jd(fgr_idx), $
	xsize=11, $
	/all_events, $
	/row, $
	/float, $
	uvalue='Wrng_y_min')

Wrng_y_max=cw_field( $
	Wcstm_third_pnl, $
	title='Y Max ', $
	value=rng_y_max_jd(fgr_idx), $
	xsize=11, $
	/all_events, $
	/row, $
	/float, $
	uvalue='Wrng_y_max')

Wcstm_done=widget_button( $
	Wcstm_btm, $
	value='           Done          ', $
	resource_name='done', $
	uvalue='Wcstm_done')

Wcstm_dfl=widget_button( $
	Wcstm_btm, $
	value='     Load Defaults        ', $
	uvalue='Wcstm_dfl')

Wcstm_draw=widget_button( $
	Wcstm_btm, $
	value='     Draw        ', $
	resource_name='draw', $
	uvalue='Wcstm_draw')

if (fgr_lst(fgr_idx).fnc_nm eq 'lon_lat_fll') or $
(fgr_lst(fgr_idx).fnc_nm eq 'lon_lat_img') or $
(fgr_lst(fgr_idx).fnc_nm eq 'pol_img') or $
(fgr_lst(fgr_idx).fnc_nm eq 'lat_lev') or $
(fgr_lst(fgr_idx).fnc_nm eq 'lon_lev') or $
(fgr_lst(fgr_idx).fnc_nm eq 'time_lat') then begin

Wbackground_base=widget_base($
	Wcstm_first_pnl, $
	/exclusive, $
	/row)

Wbackground=widget_button( $
	Wbackground_base, $
	value='Background', $
	uvalue='Wbackground')

Wcntr_ovl_base=widget_base($
	Wcstm_first_pnl, $
	/exclusive, $
	/row)

Wcntr_ovl=widget_button( $
	Wcntr_ovl_base, $
	value='Overlay contours', $
	uvalue='Wcntr_ovl')

Wusr_cntr_base=widget_base($
	Wcstm_first_pnl, $
	/exclusive, $
	/row)

Wcntr_usr=widget_button( $
	Wusr_cntr_base, $
	value='User-specified levels', $
	uvalue='Wcntr_usr')

	Wcntr_lvl_nbr=cw_field( $
		Wcstm_second_pnl, $
		title='#', $
		value=cntr_lvl_nbr, $
		xsize=2, $
		/return_events, $
		/row, $
		/integer, $
;		resource_name='cntr_lvl_nbr', $
		uvalue='Wcntr_lvl_nbr')

	Wcntr_lvl_min=cw_field( $
		title='Min', $
		Wcstm_second_pnl, $
		value=cntr_lvl_min, $
		xsize=11, $
		/float, $
		/return_events, $
		/row, $
;		resource_name='cntr_lvl_min', $
		uvalue='Wcntr_lvl_min')

	Wcntr_ntv=cw_field( $
		Wcstm_second_pnl, $
		title='Int', $
		value=cntr_ntv, $
		xsize=11.0, $
		/return_events, $
		/row, $
		/float, $
;		resource_name='cntr_ntv', $
		uvalue='Wcntr_ntv')

	Wcntr_lvl_base=widget_base( $
		Wcstm_second_pnl, $
;		resource_name='cntr_lvl_base', $
		uvalue='Wcntr_lvl_base')

	Wcntr_lvl=cw_bselector( $
		Wcntr_lvl_base, $
		strtrim(string(cntr_lvl),2), $
		label_left='Conts', $
		/return_index, $
;		resource_name='cntr_lvl', $
		uvalue='Wcntr_lvl')

		widget_control,Wcstm_base,/realize
		widget_control,Wauto_scl,set_button=auto_scl_jd(fgr_idx)

	Wusr_cntr_lvl=cw_field( $
		Wcstm_second_pnl, $
		title='User lev', $
		value=1.0e36, $
		xsize=11.0, $
		/all_events, $
		/row, $
		/float, $
;		resource_name='usr_cntr_lvl', $
		uvalue='Wusr_cntr_lvl')

widget_control,Wauto_scl,set_button=auto_scl_jd(fgr_idx)
widget_control,Wbackground,set_button=background_jd(fgr_idx)
widget_control,Wcntr_ovl,set_button=cntr_ovl_jd(fgr_idx)
widget_control,Wcntr_usr,set_button=usr_cntr_jd(fgr_idx)

endif; endif figure is a contour drawing

case fgr_lst(fgr_idx).fnc_nm of
; 19991230 fxm: I have no idea what y_ttl_rgt is used for
	'hst': begin
	Wy_ttl_rgt=cw_field( $
		Wcstm_top, $
		title='right ytitle:', $
		value=ryttl_jd(fgr_idx), $
		xsize= 40, $
		/all_events, $
		uvalue='Wy_ttl_rgt', $
		/row)

	Wbin_sz=cw_field($
		Wcstm_second_pnl, $
		title='Bin Size:', $
		value=bin_sz, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
		uvalue='Wbin_sz')

	Wbin_min=cw_field( $
		Wcstm_second_pnl, $
		title='Min Bin:', $
		value=bin_min, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
		uvalue='Wbin_min')

	Wbin_max=cw_field( $
		Wcstm_second_pnl, $
		title='Max Bin:', $
		value=bin_max, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
		uvalue='Wbin_max')

	widget_control,Wcstm_base,/realize
        end; end hst
        'lat_lev': begin
            widget_control,Wcstm_base,/realize
        end; end lat_lev
        'lon_lev': begin
            widget_control,Wcstm_base,/realize
        end; end lon_lev
        'time_lat': begin
            widget_control,Wcstm_base,/realize
        end; end time_lat
	'scat': begin
            widget_control,Wcstm_base,/realize
	end; end scat
	'znl_avg': begin
            widget_control,Wcstm_base,/realize
	end; end znl_avg
	else:
    endcase

widget_control,Winfo,set_button=info_jd(fgr_idx)
xmanager, $
	'make_cstm', $
	Wcstm_base, $
	group_leader=Wtop_lvl, $
	event_handler='cstm_event'		

endif; end creating widgets

widget_control,Wcstm_base,iconify=0 ; de-iconify when drawing

endif; endif device is X

end; end make_cstm() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End make_cstm() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin cstm_event() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cstm_event,event
@ibp_cmn.com

event_type=tag_names(event,/structure)

print,'event_type = ',event_type
print,event_type+' tags names = ',tag_names(event)
print,'event = ',event

; In case event is not handled by the below
command_sng=''

; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value

print,'user value=Customize ',usr_value

case usr_value of
'Wauto_scl': begin
	auto_scl_jd(fgr_idx)=(auto_scl_jd(fgr_idx)+1) mod 2
	form_cntr_lvl
	comp_cntr_lvl
end; end Winfo

'Wbackground': begin
	background_jd(fgr_idx)=(background_jd(fgr_idx)+1) mod 2
end; end Wbackground

'Wbin_sz': begin
	bin_sz=event.value
end; end Wbin_sz

'Wcntr_ntv': begin
	cbar_fmt='UNSET'
	cntr_lvl=0
	cntr_ntv_jd(fgr_idx)=event.value
	cntr_ntv=event.value
	comp_cntr_lvl
	command_sng='cntr_ntv = event.value = '+string(event.value)
end; end Wcntr_ntv

'Wcntr_lvl': begin

; Copy the value of the user-specified widget into this particular level	
	if usr_cntr_jd(fgr_idx) then begin
	widget_control,get_value=usr_cntr_lvl,Wusr_cntr_lvl
	cntr_lvl(event.index)=usr_cntr_lvl
	widget_control,Wcntr_lvl,/destroy
	Wcntr_lvl=cw_bselector( $
		Wcntr_lvl_base, $
		strtrim(string(cntr_lvl),2), $
		label_left='Conts', $
		/return_index, $
;		resource_name='cntr_lvl', $
		uvalue='Wcntr_lvl')
	command_sng='cntr_lvl('+auto_sng(event.index,0)+') = '+auto_sng(cntr_lvl(event.index),3)
	cntr_lvl_min=min(cntr_lvl)
	cntr_lvl_min_jd(fgr_idx)=cntr_lvl_min
	widget_control,Wcntr_lvl_min,set_value=cntr_lvl_min_jd(fgr_idx)
	comp_cntr_lvl
	endif; endi if usr_cntr_jd(fgr_idx)
end; end Wcntr_lvl

'Wcntr_ovl': begin
	cntr_ovl_jd(fgr_idx)=(cntr_ovl_jd(fgr_idx)+1) mod 2
end; end Winfo

'Wcstm_dfl': begin
;	rng_x_min_jd(fgr_idx)=very_big_nbr
;	widget_control,Wrng_x_min,set_value=rng_x_min_jd(fgr_idx)
;	max_rng_x_jd(fgr_idx)=very_big_nbr
;	widget_control,Wmax_rng_x,set_value=max_rng_x_jd(fgr_idx)
;	rng_y_min_jd(fgr_idx)=very_big_nbr
;	widget_control,Wrng_y_min,set_value=rng_y_min_jd(fgr_idx)
;	rng_y_max_jd(fgr_idx)=very_big_nbr
;	widget_control,Wrng_y_max,set_value=rng_y_max_jd(fgr_idx)
	case fgr_lst(fgr_idx).fnc_nm of
	'lon_lat_fll': begin
		ttl_jd(fgr_idx)=slice_sng+' '+fld_sng+', '+month_sng
		xttl_jd(fgr_idx)=src_sng+' '+year_sng
		form_cntr_lvl
		comp_cntr_lvl
		unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
	end; end lon_lat_fll
	'lon_lat_img': begin
		ttl_jd(fgr_idx)=slice_sng+' '+fld_sng+', '+month_sng
		xttl_jd(fgr_idx)=src_sng+' '+year_sng
		form_cntr_lvl
		comp_cntr_lvl
		unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
	end; end lon_lat_img
	'pol_img': begin
		ttl_jd(fgr_idx)=slice_sng+' '+fld_sng+', '+month_sng
		xttl_jd(fgr_idx)=src_sng+' '+year_sng
		form_cntr_lvl
		comp_cntr_lvl
		unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
	end; end pol_img
	'lat_lev': begin
		ttl_jd(fgr_idx)=src_sng+' '+fld_sng+', '+month_sng
		result=strpos(fl_in,'spcp_')
		if result ne -1 then result=strpos(fl_in,'omega')
		if result ne -1 then ttl_jd(fgr_idx)='!7D!5'+sym_sng+' '+abb_sng+' ['+src_sng+'] '+month_sng
		result=strpos(fl_in,'amip')
		if result ne -1 then result=strpos(fl_in,'omega')
		if result ne -1 then ttl_jd(fgr_idx)='!7D!5'+sym_sng+' '+abb_sng+' ['+src_sng+'] '+month_sng
		xttl_jd(fgr_idx)='Latitude !7u!5 (!E!12_!5!N)'
		yttl_jd(fgr_idx)='Hybrid Pressure !7g!5 (mb)'
		form_cntr_lvl
		comp_cntr_lvl
		unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
	end; end lat_lev
	'lon_lev': begin
		ttl_jd(fgr_idx)=src_sng+' '+fld_sng+', '+month_sng
		result=strpos(fl_in,'spcp_')
		if result ne -1 then result=strpos(fl_in,'omega')
		if result ne -1 then ttl_jd(fgr_idx)='!7D!5'+sym_sng+' '+abb_sng+' ['+src_sng+'] '+month_sng
		result=strpos(fl_in,'amip')
		if result ne -1 then result=strpos(fl_in,'omega')
		if result ne -1 then ttl_jd(fgr_idx)='!7D!5'+sym_sng+' '+abb_sng+' ['+src_sng+'] '+month_sng
		xttl_jd(fgr_idx)='Longitude !7k!5 (!E!12_!5!E)'
		yttl_jd(fgr_idx)='Hybrid Pressure !7g!5 (mb)'
		form_cntr_lvl
		comp_cntr_lvl
		unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
	end; end lon_lev
	'time_lat': begin
		ttl_jd(fgr_idx)=src_sng+' '+fld_sng+', '+month_sng
		xttl_jd(fgr_idx)='Time'
		yttl_jd(fgr_idx)='Latitude !7u!5 (!E!12_!5!N)'
		form_cntr_lvl
		comp_cntr_lvl
		unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
	end; end time_lat
	'hst': begin
		ttl_jd(fgr_idx)=abb_sng+', '+month_sng
		xttl_jd(fgr_idx)=slice_sng+' '+sym_sng+' ('+unit_jd(fgr_idx)+')'
		yttl_jd(fgr_idx)='Frequency (%)'
		ryttl_jd(fgr_idx)='!5Cumulative Frequency (%)'
		widget_control,Wy_ttl_rgt,set_value=ryttl_jd(fgr_idx)
		unit_jd(fgr_idx)=unit_sng
	end; end hst
	'scat': begin
		ttl_jd(fgr_idx)=sv_abb_sng+' !8vs.!5 '+abb_sng+', '+month_sng+' '+year_sng
		misc_jd(fgr_idx)=src_sng+' '+rgn_lst(rgn_idx).string
		xttl_jd(fgr_idx)=sv_slc_sng+' '+sv_sym_sng+' ('+sv_unit_sng+')'
		yttl_jd(fgr_idx)=slice_sng+' '+sym_sng+' ('+unit_sng+')'
		unit_jd(fgr_idx)=unit_sng
	end; end scat
	'znl_avg': begin
		if data_nbr eq sv_data_nbr then begin
		if sv_abb_sng eq abb_sng then begin
			ttl_jd(fgr_idx)=sv_abb_sng+', '+month_sng+' '+year_sng
			xttl_jd(fgr_idx)=slice_sng+' '+sym_sng+' ('+unit_sng+')'
		endif else begin
			ttl_jd(fgr_idx)=sv_abb_sng+' !8vs.!5 '+abb_sng+', '+month_sng+' '+year_sng
			xttl_jd(fgr_idx)=slice_sng+' '+sym_sng+', '+sv_slc_sng+' '+sv_sym_sng+' ('+unit_sng+')'
		endelse
		endif else begin
			ttl_jd(fgr_idx)=abb_sng+', '+month_sng+' '+year_sng
			xttl_jd(fgr_idx)=slice_sng+' '+sym_sng+' ('+unit_sng+')'
		endelse
		yttl_jd(fgr_idx)='!5Latitude (!5!E!12_!N!5)'
		unit_jd(fgr_idx)=unit_sng
	end; end znl_avg
endcase; endcase over fgr type
widget_control,Wttl,set_value=ttl_jd(fgr_idx)
widget_control,Wx_ttl,set_value=xttl_jd(fgr_idx)
widget_control,Wy_ttl,set_value=yttl_jd(fgr_idx)
widget_control,Wmisc,set_value=misc_jd(fgr_idx)
widget_control,Wunit,set_value=unit_jd(fgr_idx)
end; end Wcstm_dfl

'Wcstm_draw': begin

	drw_pretty_picture	

end; end Wcstm_draw

'Wcstm_done': begin
        widget_control,Wcstm_base,/destroy
end; end Wcstm_done

'Winfo': begin
	info_jd(fgr_idx)=(info_jd(fgr_idx)+1) mod 2
end; end Winfo

'Wmisc': begin
	misc_jd(fgr_idx)=event.value
end; end Wmisc

'Wbin_max': begin
	bin_max=event.value
end; end Wbin_max

'Wmax_rng_x': begin
	max_rng_x_jd(fgr_idx)=event.value
end; end Wmax_rng_x

'Wrng_y_max': begin
	rng_y_max_jd(fgr_idx)=event.value
end; end Wrng_y_max

'Wbin_min': begin
	bin_min=event.value
end; end Wbin_min

'Wcntr_lvl_min': begin
	cbar_fmt='UNSET'
	cntr_lvl=0
	cntr_lvl_min_jd(fgr_idx)=event.value
	cntr_lvl_min=event.value
	comp_cntr_lvl
	command_sng='cntr_lvl_min = event.value = '+string(event.value)
end; end Wcntr_lvl_min

'Wrng_x_min': begin
	rng_x_min_jd(fgr_idx)=event.value
end; end Wrng_x_min

'Wrng_y_min': begin
	rng_y_min_jd(fgr_idx)=event.value
end; end Wrng_y_min

'Wcntr_lvl_nbr': begin
	cbar_fmt='UNSET'
	cntr_lvl=0
	cntr_lvl_nbr_jd(fgr_idx)=event.value
	cntr_lvl_nbr=event.value
	comp_cntr_lvl
	command_sng='cntr_lvl_nbr = event.value = '+string(event.value)
end; end Wcntr_lvl_nbr

'Wy_ttl_rgt': begin
	ryttl_jd(fgr_idx)=event.value
end; end Wy_ttl_rgt

'Wttl': begin
	ttl_jd(fgr_idx)=event.value
end; end Wttl

'Wunit': begin
	unit_jd(fgr_idx)=event.value
end; end Wunit

'Wcntr_usr': begin
	usr_cntr_jd(fgr_idx)=(usr_cntr_jd(fgr_idx)+1) mod 2
end; end Winfo

'Wusr_cntr_lvl': begin
	command_sng='usr_cntr_lvl = event.value = '+string(event.value)
end; end Wusr_cntr_lvl

'Wx_ttl': begin
	xttl_jd(fgr_idx)=event.value
end; end Wx_ttl

'Wy_ttl': begin
	yttl_jd(fgr_idx)=event.value
end; end Wy_ttl

endcase; endcase over event type
end; end cstm_event()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End cstm_event() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
