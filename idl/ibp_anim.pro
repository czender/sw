;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin make_anim() 
; If the current device is X then this routine will create the 
; widgets and install the event handler for the XY Range panel if
; they don't already exist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_anim

@ibp_cmn.com

if !d.name eq 'X' then begin

; If the window is already up then OK else start it up.
; be careful because old widget ID's stay in common memory 
; after crashes!
widget_xst=1
if n_elements(Wanim_base) ne 0 then begin
	widget_xst=widget_info(Wanim_base,/valid_id)
endif

if ((widget_xst ne 1) or (n_elements(Wanim_base) eq 0)) then begin

; If we're here then we can assume the animation is starting from scratch

anim_idx=0

; Create anim widget hierarchy
Wanim_base=widget_base( $
	group_leader=Wtop_lvl, $
	title='Animation', $
	/row)

Wanim_nbr=cw_field( $
		Wanim_base, $
		title='# Frames ', $
		value=anim_nbr, $
		xsize=5, $
		/return_events, $
		/row, $
		/integer, $
;		resource_name='anim_nbr', $
		uvalue='Wanim_nbr')

Wanim_load=widget_button( $
		Wanim_base, $
		value='Load Image', $
;		resource_name='anim_load', $
		uvalue='Wanim_load')

Wanim_lev=widget_button( $
		Wanim_base, $
		value='Auto Anim Lev', $
;		resource_name='anim_lev', $
		uvalue='Wanim_lev')

Wanim_time=widget_button( $
		Wanim_base, $
		value='Auto Anim Time', $
;		resource_name='anim_time', $
		uvalue='Wanim_time')

Wanim_play=widget_button( $
		Wanim_base, $
		value='Play', $
;		resource_name='anim_play', $
		uvalue='Wanim_play')

Wanim_done=widget_button( $
		Wanim_base, $
		value='Done', $
;		resource_name='anim_done', $
		uvalue='Wanim_done')

widget_control,Wanim_base,/realize

xmanager, $
	'make_anim', $
	Wanim_base, $
	group_leader=Wtop_lvl, $
	event_handler='anim_event'		

Wprojector_base=widget_base( $
	group_leader=Wanim_base, $
	title='Projector', $
	/row)

Wprojector=cw_animate( $
		Wprojector_base, $
		lon_lat_wdw.x_size, $
		lon_lat_wdw_y_sz, $
		anim_nbr, $
;		resource_name='projector', $
		uvalue='Wprojector') 

widget_control,Wprojector_base,/realize

xmanager, $
	'make_anim', $
	Wprojector_base, $
	group_leader=Wanim_base, $
	event_handler='anim_event'		

endif; end creating widgets

widget_control,Wanim_base,iconify=0 ; de-iconify when drawing

endif; endif device is X

end; end make_anim() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End make_anim() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin anim() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro anim_event,event

@ibp_cmn.com

event_type=tag_names(event,/structure)

print,'event_type = ',event_type
print,event_type+' tags names = ',tag_names(event)
print,'event = ',event

; In case event is not handled by the below
command_sng=''

; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value

if (usr_value eq 'Wanim_nbr') then begin
	anim_nbr=event.value
	command_sng='anim_nbr = event.value = '+string(event.value)
endif; endif Wanim_nbr

if (usr_value eq 'Wanim_lev') then begin

; Only animate in levels if there is more than one level.
foo=min(abs(lev-lev_min),lev_min_idx)
foo=min(abs(lev-lev_max),lev_max_idx)
loop_step=1
if lev_min_idx gt lev_max_idx then loop_step=-1

if abs(lev_min_idx-lev_max_idx) gt 0 then begin

; The documentation for widget_list events seems to be incorrect.
; Also, and more inexplicable, is why naming the button event
; structure widget_button causes the top_lvl event loop to
; crash whenever buttons are pressed, and similarly for lists.
; e.g.,
; % Array dimensions must be greater than 0.
; % Execution halted at XMANAGER </opt/idl3.6/lib/prog/xmanager.pro( 457)> (WIDGET_EVENT).
;
; Do i not name pull down structures cw_pdmenu_s?
; Why must i create an identical structure to widget_button,
; but call it something else like button_event, in order to 
; stop the crashes?
;
foo={button_event,id:0L,top:0L,handler:0L,select:0L}
foo={list_event,id:0L,top:0L,handler:0L,index:0L,clicks:0L}
read_event={button_event,Wread,Wtop_lvl,Wtop_lvl,1}
drw_event={button_event,Wdraw,Wtop_lvl,Wtop_lvl,1}
anim_load_event={button_event,Wanim_load,Wanim_base,Wanim_base,1}
print,'lev_min_idx = ',lev_min_idx
print,'lev_max_idx = ',lev_max_idx
tmp_vrt_slc=vrt_slc
for vrt_slc=lev_min_idx,lev_max_idx,loop_step do begin
	print,'loading vrt_slc = ',vrt_slc
	vrt_slc_event={list_event,Wvrt_slc,Wtop_lvl,Wtop_lvl,vrt_slc,1}

; It appears the /nocopy keyword can foul things up so don't use it
	widget_control,Wvrt_slc,send_event=vrt_slc_event
	widget_control,Wread,send_event=read_event
	widget_control,Wdraw,send_event=drw_event
;	widget_control,Wanim_load,send_event=anim_load_event
;	wait,5
	print,'loading current image'
	cw_animate_load,Wprojector,frame=anim_idx,/track,window=lon_lat_wdw
	anim_idx=anim_idx+1
endfor; env loop over levels
vrt_slc=tmp_vrt_slc
command_sng='animate levels'
endif else begin
command_streing='Need a multi-level field to animate'
endelse
endif; endif Wanim_lev

if (usr_value eq 'Wanim_load') then begin
cw_animate_load, $
	Wprojector, $
	frame=anim_idx, $
	/track, $
	window=lon_lat_wdw
anim_idx=anim_idx+1
command_sng='cw_animate_load'
endif; endif Wanim_load

if (usr_value eq 'Wanim_time') then begin
	command_sng='animate in time'
endif; endif Wanim_time

if (usr_value eq 'Wanim_play') then begin
	rate=10
	cw_animate_run,Wprojector,rate,nframes=anim_idx
	command_sng='cw_animate_run'
endif; endif Wanim_play

if (usr_value eq 'Wprojector') then begin
	widget_control,Wprojector_base,/destroy
	command_sng='widget_control,Wprojector_base,/destroy'
endif; endif Wprojector

if (usr_value eq 'Wanim_done') then begin
	widget_control,Wanim_base,/destroy
	command_sng='widget_control,Wanim_base,/destroy'
endif; endif Wrng_y_max

if command_sng ne '' then print,command_sng

end; end anim() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End anim() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




