; $Id$

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin scat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro scat_event,event

@ibp_cmn.com

event_type=tag_names(event,/structure)
;print,event_type,event

; In case event is not handled by the below
command_sng=''

; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value

if (usr_value eq 'Wscat_gph') then begin
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
		inverse_scat
	endif
endif; endif Wscat_gph

if command_sng ne '' then print,command_sng

end; end scat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End scat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin make_scat() 
; If the current device is X then this routine will create the 
; widgets and install the event handler for the scatterplot graphs if
; they don't already exist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_scat

@ibp_cmn.com

if !d.name eq 'X' then begin

; If the window is already up then OK else start it up.
; be careful because old widget ID's stay in common memory 
; after crashes!

widget_xst=1
if n_elements(Wscat_gph) ne 0 then begin
	widget_xst=widget_info(Wscat_gph,/valid_id)
endif

if ((widget_xst) and (n_elements(Wscat_gph) ne 0)) then begin
curr_sz=indgen(2)
widget_control,Wscat_base,tlb_get_size=curr_sz

; If the window has been resized, then destroy the widget hierarchy
; so that the widget can be recreated at the new size.

wset,scat_wdw

if ((scat_wdw_x_sz ne curr_sz(0)) or (scat_wdw_y_sz ne curr_sz(1))) then begin
	scat_wdw_x_sz=curr_sz(0)
	scat_wdw_y_sz=curr_sz(1)
	widget_control,Wscat_base,/destroy
	widget_xst=0
endif

endif

if ((widget_xst ne 1) or (n_elements(Wscat_gph) eq 0)) then begin

; Create the scatterplot widget hierarchy

Wscat_base=widget_base( $
	group_leader=Wtop_lvl, $
	title='Scatter Plot')
Wscat_gph=widget_draw( $
	Wscat_base, $
		xsize=scat_wdw_x_sz, $
		ysize=scat_wdw_y_sz, $
	/button_events, $
	uvalue='Wscat_gph')

widget_control,Wscat_base,/realize
widget_control,get_value=scat_wdw,Wscat_gph

; Call the event manager for the new base widgets
xmanager, $
	'make_scat', $
	Wscat_base, $
	group_leader=Wtop_lvl, $
	event_handler='scat_event'		

endif; end creating widgets

widget_control,Wscat_base,iconify=0 ; de-iconify when drawing
wset,scat_wdw

endif; endif device is X

end; end make_scat() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End make_scat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure Scatterplot 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro scat
@ibp_cmn.com

erase

; Make sure there are two equally sized fields selected

if n_elements(sv_data) eq 0 then begin
	print,'Must save a field first'
	goto, end_of_procedure
endif
if n_elements(sv_data) ne n_elements(data) then begin
	print,'Fields are of differing sizes, check your data...'
	goto, end_of_procedure
endif

; find the maxima after the masked values have been removed

good_idx=where(data lt very_big_nbr,count)
good_data=data(good_idx)
good_sv_data=sv_data(good_idx)
good_idx=where(good_sv_data lt very_big_nbr,count)
good_data=good_data(good_idx)
good_sv_data=good_sv_data(good_idx)
good_data_nbr=n_elements(good_data)

min_good_data=min(good_data)
max_good_data=max(good_data)
min_good_sv_data=min(good_sv_data)
max_good_sv_data=max(good_sv_data)

rng_x_min=min_good_sv_data
max_rng_x=max_good_sv_data
rng_y_min=min_good_data
rng_y_max=max_good_data

if (strlen(unit_jd(fgr_idx)) le 0) then unit_jd(fgr_idx)=unit_sng
if (strlen(ttl_jd(fgr_idx)) le 0) then ttl_jd(fgr_idx)=sv_abb_sng+' !8vs.!5 '+abb_sng
if (strlen(misc_jd(fgr_idx)) le 0) then misc_jd(fgr_idx)=''
if (strlen(xttl_jd(fgr_idx)) le 0) then xttl_jd(fgr_idx)=sv_slc_sng+' '+sv_sym_sng+' ('+sv_unit_sng+')'
if (strlen(yttl_jd(fgr_idx)) le 0) then yttl_jd(fgr_idx) =slice_sng+' '+sym_sng+' ('+unit_sng+')'

; Set defaults for x,y ranges then override with any user specified values

if (rng_x_min_jd(fgr_idx) lt very_big_nbr) then rng_x_min=rng_x_min_jd(fgr_idx)
if (max_rng_x_jd(fgr_idx) lt very_big_nbr) then max_rng_x=max_rng_x_jd(fgr_idx)
if (rng_y_min_jd(fgr_idx) lt very_big_nbr) then rng_y_min=rng_y_min_jd(fgr_idx)
if (rng_y_max_jd(fgr_idx) lt very_big_nbr) then rng_y_max=rng_y_max_jd(fgr_idx)

;print,very_big_nbr
;print,good_sv_data
;print,good_data
if !d.name ne 'X' then chr_sz=2.0 else chr_sz=1.0

plot, $
	good_sv_data, $
	good_data, $
	max_value=very_big_nbr, $
;	tit=ttl_jd(fgr_idx), $
	xtit=xttl_jd(fgr_idx), $
	ytit=yttl_jd(fgr_idx), $
	xstyle=1, $
	ystyle=1, $
	xrange=[rng_x_min,max_rng_x], $
	yrange=[rng_y_min,rng_y_max], $
	thick=2.0, $
	charsize=chr_sz, $
	charthick=2.0, $
	xmargin=[7,2], $ ;[10,3] is [left,right] default
;	ymargin=[3,2], $  ;[4,2] is [bottom,top] default
	ymargin=[3,3], $  ;[4,2] is [bottom,top] default
	psym=1

; Save the axes ranges for the next time through
rng_x_min_jd(fgr_idx)=!x.crange(0)
max_rng_x_jd(fgr_idx)=!x.crange(1)
rng_y_min_jd(fgr_idx)=!y.crange(0)
rng_y_max_jd(fgr_idx)=!y.crange(1)

; Perform a statistical analysis of the correlation between the two datasets.
fit_coeffs=poly_fit( $
		good_sv_data, $
		good_data, $
		1, $
		fitted_ords)

if (info_jd(fgr_idx) gt 0) then begin
oplot, $
	[min_good_sv_data, $
	max_good_sv_data], $
	[fit_coeffs(0)+min_good_sv_data*fit_coeffs(1), $
	fit_coeffs(0)+max_good_sv_data*fit_coeffs(1)], $
	thick=3.0, $
	linestyle=0	

	xyouts,0.98,.025,sv_sub_ttl_sng+' '+sv_sym_sng+' !8vs.!5 '+sub_ttl_sng+' '+sym_sng+' in '+rgn_lst(rgn_idx).string+': '+lon_sng+', '+lat_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
	sys_time=systime(0) 
	xyouts,1.0,0.0,usr_sng+'@'+hostname_sng+' '+sys_time,size=dbg_txt_sz,orientation=0.0,alignment=1.0,/NORMAL
endif; endif info
xyouts,.55,0.96,ttl_jd(fgr_idx),alignment=0.5,orientation=0.0,charthick=2.0,size=3.0,/NORMAL
xyouts,.55,0.91,misc_jd(fgr_idx),alignment=0.5,orientation=0.0,charthick=2.0,size=2.0,/NORMAL

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin linear (Pearson) regression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
matrix=fltarr(2,good_data_nbr)
matrix(0,*)=good_sv_data
matrix(1,*)=good_data
pearson_mtx=correl_matrix(matrix)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End linear (Pearson) regression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Non-parametric (Spearman) correlation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
matrix=fltarr(2,good_data_nbr)
matrix(0,*)=good_sv_data
matrix(1,*)=good_data
spearman,matrix,spearman_mtx
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Non-parametric (Spearman) correlation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

b_sng=auto_sng(abs(fit_coeffs(0)),2)
m_sng=auto_sng(fit_coeffs(1),2)
if fit_coeffs(0) ge 0. then sign_sng=' + ' else sign_sng=' - '
fit_sng=sym_sng+' = '+m_sng+' '+sv_sym_sng+sign_sng+b_sng
pearson_sng=string(format='(F5.2)',pearson_mtx(0,1))
pearson_sng='!8r!5 = '+pearson_sng
spearman_sng=string(format='(F5.2)',spearman_mtx(0,1))
spearman_sng='!8r!Is!N!5 = '+spearman_sng

if (info_jd(fgr_idx) gt 0) then begin
	xyouts,1.0,.025,'Least squares fit: '+fit_sng+'; Linear '+pearson_sng+', Non-parametric '+spearman_sng,size=dbg_txt_sz,orientation=90.0,alignment=0.0,/NORMAL
	;xyouts,0.2,0.75,fit_sng,size=2.0,orientation=90.0,alignment=0.0,/NORMAL
	;xyouts,0.2,0.75,pearson_sng,size=2.0,orientation=90.0,alignment=0.0,/NORMAL
	;xyouts,0.2,0.70,spearman_sng,size=2.0,orientation=90.0,alignment=0.0,/NORMAL
endif; endif info

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Full multiple regression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
matrix=fltarr(1,good_data_nbr)
matrix(0,*)=good_sv_data
weights=fltarr(good_data_nbr)
weights(*)=1
m_term = regress( $
		matrix, $
		reform(good_data,good_data_nbr), $
		weights, $
		/relative_weight, $
		fitted_ords, $
		b_term, $
		sigma, $
		f_test, $
		pearson, $
		multi_pearson, $
		chi_squared)

print,'b_term = ',b_term
print,'m_term = ',m_term
print,'sigma = ',sigma
print,'f_test = ',f_test
print,'pearson = ',pearson
print,'multi_pearson = ',multi_pearson
print,'chi_squared = ',chi_squared
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Full multiple regression
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end_of_procedure: foo=1
end; end scat()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure Scatterplot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Inverse Scatterplot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro inverse_scat
@ibp_cmn.com

x_min=ptr_psn(0)
y_min=ptr_psn(3)
x_max=ptr_psn(2)
y_max=ptr_psn(1)

; Search for all points lying within the bounded region.
; The simplest way seems to be to convert the coords to lat/lon indices.
x_idx=where((sv_data ge x_min) and (sv_data le x_max),count)
if count eq 0 then return
y_idx=where((data ge y_min) and (data le y_max),count)
if count eq 0 then return

nbr_x_idx=n_elements(x_idx)
nbr_y_idx=n_elements(y_idx)

indices=0L; bogus hack
for i=0,nbr_x_idx-1 do begin
	foo=where(y_idx eq x_idx(i),count)
	if count eq 1 then indices=[indices,x_idx(i)]
endfor
nbr_scat_data=n_elements(indices)-1
if nbr_scat_data gt 0 then indices=indices(1:nbr_scat_data) else return

;print,'(x_min,y_max) = ',x_min,y_max
;print,'(x_max,y_min) = ',x_max,y_min
;print,'indices = ',indices

; indices is now an array of longword indices into the original data
; we can now invert it to find the latitude and longitude of the data.
lat_idx=indgen(nbr_scat_data)
lon_idx=indgen(nbr_scat_data)
for i=0,nbr_scat_data-1 do begin
	lon_idx(i)=indices(i) mod lon_nbr
	lat_idx(i)=fix(indices(i)/lon_nbr)

	print,'lon = '+string(format='(F7.2)',lon(lon_idx(i)))+ $
		', lat = '+string(format='(F7.2)',lat(lat_idx(i)))+ $
		', '+sv_fld_nm+' = '+ $
		string(format='(F7.2)',sv_data(lon_idx(i),lat_idx(i)))+ $
		', '+fld_nm+' = '+string(format='(F7.2)', $
		data(lon_idx(i),lat_idx(i)))
endfor

wset,lon_lat_wdw
outline_points,nbr_scat_data,lon_idx,lat_idx

end; end inverse_scat()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Inverse Scatterplot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Outline Points
; This routine assumes, for the first time that i'm aware of in this
; program, that the longitudes are regularly spaced.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro outline_points,nbr_scat_data,lon_idx,lat_idx
@ibp_cmn.com

print,'nbr_scat_data = ',nbr_scat_data
print,'lon_idx = ',lon_idx
print,'lat_idx = ',lat_idx

; The map coordinates should be reset because IDL only seems to remember
; the coordinates of the last window graphed in, which would be the
; histogram or scatterplot in this place.
map_set, $
	lat_ctr, $		;[dgr] Latitude at center of map
	lon_ctr, $		;[dgr] Longitude at center of map
	0.0, $			;[dgr] North rotation angle
	limit=map_lmt_dgr, $;[latmin,lonmin,latmax,lonmax] drawing box
	/cylindrical, $		;cylindrical projection
	xmargin=x_map_cbar_mrg, $
	ymargin=y_map_cbar_mrg, $
	/noerase, $
	/noborder 		;don't draw a border around the grid

lon_spacing=lon(1)-lon(0)
lon_spacingd2=lon_spacing/2.
for i=0,nbr_scat_data-1 do begin
; Find the bounding coordinates of the box containing the data point
	box_lon_ctr=lon(lon_idx(i))
	box_lat_ctr=lat(lat_idx(i))
	if lat_idx(i) ne 0 then $
		lat_spacing=lat(lat_idx(i))-lat(lat_idx(i)-1) else $
		lat_spacing=lat(lat_idx(i))-lat_min
	lat_spacingd2=lat_spacing/2.

; Here it's assumed that the data point is centered within the box to be
; drawn. This is not correct but close enough for a first cut.

	plots, $
		[box_lon_ctr-lon_spacingd2, $ ;lft_box
		box_lon_ctr-lon_spacingd2, $ ;lft_box
		box_lon_ctr+lon_spacingd2, $ ;rgt_box
		box_lon_ctr+lon_spacingd2, $ ;rgt_box
		box_lon_ctr-lon_spacingd2], $ ;lft_box
		[box_lat_ctr-lat_spacingd2, $ ;bot_box
		box_lat_ctr+lat_spacingd2, $ ;top_box
		box_lat_ctr+lat_spacingd2, $ ;top_box
		box_lat_ctr-lat_spacingd2, $ ;bot_box
		box_lat_ctr-lat_spacingd2], $ ;bot_box
		thick=2.0, $
		/data, $
		linestyle=0	

endfor

end; end outline_points()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Outline Points
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
