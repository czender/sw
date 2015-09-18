; $Id$

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin hst() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hst_event,event

@ibp_cmn.com

event_type=tag_names(event,/structure)
;print,event_type,event

; In case event is not handled by the below
command_sng=''

; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value

if (usr_value eq 'Whst_gph') then begin
	data_crd=convert_coord( $
		event.x, $
		event.y, $
		/device, $
		/to_data)
	if ((event.type eq 0) and (event.press eq 1)) then begin
; left button was pressed
		inverse_hst,data_crd
	endif
endif; endif Whst_gph

if command_sng ne '' then print,command_sng

end; end hst() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End hst() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin make_hst() 
; If the current device is X then this routine will create the 
; widgets and install the event handler for the histogram graphs if
; they don't already exist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_hst

@ibp_cmn.com

if !d.name eq 'X' then begin

; If the window is already up then OK else start it up.
; be careful because old widget ID's stay in common memory 
; after crashes!
widget_xst=1
if n_elements(Whst_gph) ne 0 then begin
	widget_xst=widget_info(Whst_gph,/valid_id)
endif

if ((widget_xst) and (n_elements(Whst_gph) ne 0)) then begin
curr_sz=indgen(2)
widget_control,Whst_base,tlb_get_size=curr_sz

; If the window has been resized, then destroy the widget hierarchy
; so that the widget can be recreated at the new size.

wset,hst_wdw

if ((hst_wdw_x_sz ne curr_sz(0)) or (hst_wdw_y_sz ne curr_sz(1))) then begin
	hst_wdw_x_sz=curr_sz(0)
	hst_wdw_y_sz=curr_sz(1)
	widget_control,Whst_base,/destroy
	widget_xst=0
endif

endif

if ((widget_xst ne 1) or (n_elements(Whst_gph) eq 0)) then begin
; Create the histogram widget hierarchy
Whst_base=widget_base( $
	group_leader=Wtop_lvl, $
	title='Histogram')
Whst_gph=widget_draw( $
	Whst_base, $
		xsize=hst_wdw_x_sz, $
		ysize=hst_wdw_y_sz, $
	/button_events, $
	uvalue='Whst_gph')

widget_control,Whst_base,/realize
widget_control,get_value=hst_wdw,Whst_gph

; Call the event manager
xmanager, $
	'make_hst', $
	Whst_base, $
	group_leader=Wtop_lvl, $
	event_handler='hst_event'		

endif; end creating widgets

widget_control,Whst_base,iconify=0 ; de-iconify when drawing
wset,hst_wdw

endif; endif device is X

end; end make_hst() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End make_hst() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Fgr Histogram
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hst
@ibp_cmn.com

erase

; find the maxima after the masked values have been removed
good_data=where(data lt very_big_nbr,count)
if count ne -1 then good_data=data(good_data) else good_data=data
max_good_data=max(good_data)
good_data_nbr=n_elements(good_data)

if n_elements(bin_sz) eq 0 then bin_sz=cntr_ntv
if n_elements(bin_max) eq 0 then bin_max=max_good_data
if n_elements(bin_min) eq 0 then bin_min=data_min

if bin_sz eq 0. then begin
	print,'Histogram bin_sz not allowed to be zero'
	goto,end_of_procedure
endif
if bin_max eq bin_min then begin
	print,'Histogram bin_max = bin_min = ',bin_max,' is not allowed'
	goto,end_of_procedure
endif

if bin_min ne data_min then begin
	outside=where(data lt bin_min,count)
	if count ne 0 then print,'points smaller than minimum bin: ',count
endif
if bin_max ne data_max then begin
	outside=where(data gt bin_max,count)
	if count ne 0 then print,'points larger than maximum bin: ',count
endif

hst=histogram(    $
		data, $
		binsize=bin_sz, $
		omax=hst_max, $
		omin=hst_min, $
		max=bin_max, $
		min=bin_min, $
		reverse_indices=inverse_data)

; Scale the histogram to frequency
frequency_axis=100.*hst/float(good_data_nbr)

print,'hst_min = ',hst_min
print,'hst_max = ',hst_max

bin_nbr=(size(hst))(1)
abc_axis=bin_sz*(findgen(bin_nbr)+.5)+hst_min
print,'There are ',bin_nbr,' bins covering from ',hst_min,' to ',hst_min+bin_nbr*bin_sz

cum_freq_axis=fltarr(bin_nbr)
cum_freq_axis(0)=frequency_axis(0)
for i=1,bin_nbr-1 do begin
	cum_freq_axis(i)=cum_freq_axis(i-1)+frequency_axis(i)
endfor

if (strlen(unit_jd(fgr_idx)) le 0) then unit_jd(fgr_idx)=unit_sng
if (strlen(ttl_jd(fgr_idx)) le 0) then ttl_jd(fgr_idx)=abb_sng
if (strlen(xttl_jd(fgr_idx)) le 0) then xttl_jd(fgr_idx)=slice_sng+' '+sym_sng+ '('+unit_jd(fgr_idx)+')'
if (strlen(yttl_jd(fgr_idx)) le 0) then yttl_jd(fgr_idx) ='Frequency (%)'
if (strlen(ryttl_jd(fgr_idx)) le 0) then ryttl_jd(fgr_idx) ='!5Cumulative Frequency (%)'

; Set defaults for x,y ranges then override with any user specified values
rng_x_min=data_min
max_rng_x=data_max
rng_y_min=min(frequency_axis)
rng_y_max=max(frequency_axis)

; Set defaults for x,y ranges then override with any user specified values
if (rng_x_min_jd(fgr_idx) lt very_big_nbr) then rng_x_min=rng_x_min_jd(fgr_idx)
if (max_rng_x_jd(fgr_idx) lt very_big_nbr) then max_rng_x=max_rng_x_jd(fgr_idx)
if (rng_y_min_jd(fgr_idx) lt very_big_nbr) then rng_y_min=rng_y_min_jd(fgr_idx)
if (rng_y_max_jd(fgr_idx) lt very_big_nbr) then rng_y_max=rng_y_max_jd(fgr_idx)
if !d.name ne 'X' then chr_sz=2.0 else chr_sz=1.0

plot, $
	abc_axis, $
	frequency_axis, $
	tit=ttl_jd(fgr_idx), $
	xtit=xttl_jd(fgr_idx), $
	ytit=yttl_jd(fgr_idx), $
	xrange=[rng_x_min,max_rng_x], $
	yrange=[rng_y_min,rng_y_max], $
	xstyle=0, $
	ystyle=8, $
	thick=3.0, $
	charsize=chr_sz, $
	xmargin=[7,7.5], $  ;[10,3] is [left,right] default
	ymargin=[3,2], $  ;[4,2] is [bottom,top] default
	psym=10

; Save the axes ranges for the next time through
rng_x_min_jd(fgr_idx)=!x.crange(0)
max_rng_x_jd(fgr_idx)=!x.crange(1)
rng_y_min_jd(fgr_idx)=!y.crange(0)
rng_y_max_jd(fgr_idx)=!y.crange(1)
if !d.name ne 'X' then chr_sz=2.0 else chr_sz=1.0

plots, $
	[hst_min, $
	hst_min, $
	abc_axis(0)], $
	[0.0, $
	frequency_axis(0), $
	frequency_axis(0)], $
	linestyle=0, $
	thick=3.0, $
	/data

plots, $
	[hst_min+bin_nbr*bin_sz, $
	hst_min+bin_nbr*bin_sz, $
	abc_axis(bin_nbr-1)], $
	[0.0, $
	frequency_axis(bin_nbr-1), $
	frequency_axis(bin_nbr-1)], $
	linestyle=0, $
	thick=3.0, $
	/data

plot, $
	abc_axis, $
	cum_freq_axis, $
	xstyle=4, $
	ystyle=4, $
	/ynozero, $
	thick=3.0, $
	charsize=chr_sz, $
	linestyle=2, $
	xmargin=[7,7.5], $ ;[10,3] is [left,right] default
	ymargin=[3,2], $  ;[4,2] is [bottom,top] default
	/noerase

axis, $
	yaxis=1, $
	ystyle=1, $
	ytitle= ryttl_jd(fgr_idx), $
;	ytitle='!5Cumulative Frequency (%)', $
	charsize=chr_sz

possible_elements_sng=auto_sng(data_nbr,0)
plotted_elements_sng=auto_sng(round(total(hst)),0)
bin_sz_sng=auto_sng(bin_sz,1)
bin_sz_sng='!5Bin size = '+bin_sz_sng+' '+unit_jd(fgr_idx)+'.'
avg_sng=auto_sng(avg_data,1)
avg_sng='!5Average '+sym_sng+' = '+avg_sng+' '+unit_jd(fgr_idx)+'.'
mdn_sng=auto_sng(median(data),1)
mdn_sng='!5Median '+sym_sng+' = '+mdn_sng+' '+unit_jd(fgr_idx)+'.'

print,''
print,'lat_max = '+auto_sng(lat_max,1)+', lat_min = '+auto_sng(lat_min,1)+', lat_spn = '+auto_sng(abs(lat_max-lat_min),1)
yprint,'lon_max = '+auto_sng(lon_max,1)+', lon_min = '+auto_sng(lon_min,1)+', lon_spn = '+auto_sng(abs(lon_max-lon_min),1)
print,'good_data_nbr = '+auto_sng(good_data_nbr,0)+', avg_data = '+auto_sng(avg_data,2)

if (info_jd(fgr_idx) gt 0) then begin
	xyouts,0.98,.025,sub_ttl_sng+' '+sym_sng+' in '+rgn_lst(rgn_idx).string+': '+lon_sng+', '+lat_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
;
	xyouts,1.0,.025,'Used '+plotted_elements_sng+'/'+possible_elements_sng+' points.  '+bin_sz_sng+'  '+avg_sng+'  '+mdn_sng,alignment=0.0,orientation=90.0,size=dbg_txt_sz,/NORMAL
endif; endif info

sys_time=systime(0) 
xyouts,1.0,0.0,usr_sng+'@'+hostname_sng+' '+sys_time,size=dbg_txt_sz,orientation=0.0,alignment=1.0,/NORMAL

end_of_procedure: foo=1
end; end hst()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure Histogram
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Inverse Histogram
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro inverse_hst,data_crd
@ibp_cmn.com

ptr_x=data_crd(0)
ptr_y=data_crd(1)

; Find the bin indicated by the x cursor value. Assume linear bins.
bin=fix((ptr_x-hst_min)/bin_sz)

; Was the mouse clicked over a bin or not?
if ((bin lt 0) or (bin gt bin_nbr-1)) then begin
	print,'Non-existent bin'
	return
endif

; Now that we know the bin selected we can print out the data in that bin
if inverse_data(bin) ne inverse_data(bin+1) then begin
nbr_bin_data=n_elements(data(inverse_data(inverse_data(bin):inverse_data(bin+1)-1)))
print,'There are ',nbr_bin_data,' data points in bin ',bin

lat_idx=indgen(nbr_bin_data)
lon_idx=indgen(nbr_bin_data)
for datum=0,nbr_bin_data-1 do begin
	datum_idx=inverse_data(inverse_data(bin)+datum)
	lon_idx(datum)=datum_idx mod lon_nbr
	lat_idx(datum)=fix(datum_idx/lon_nbr)

	print,'lon = '+string(format='(F7.2)',lon(lon_idx(datum)))+ $
		', lat = '+string(format='(F7.2)',lat(lat_idx(datum)))+ $
		', '+fld_nm+' = '+ $
		string(format='(F7.2)',data(lon_idx(datum),lat_idx(datum)))
;
;	print,'datum = ',datum,', datum_idx = ',datum_idx,', data = ',data(datum_idx)
;	print,'data(',lon_idx,',',lat_idx,') = ',data(lon_idx,lat_idx)
;	print,'lat = ',lat(lat_idx),' degrees, lon = ',lon(lon_idx),' degrees.'
endfor

endif else return; end if there are data in this selected bin

wset,lon_lat_wdw
outline_points,nbr_bin_data,lon_idx,lat_idx

end; end inverse_hst()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Inverse Histogram
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
