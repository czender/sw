; $Id$

; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro

; Purpose: Routines to visualize global datasets.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin lon_lat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lon_lat_event,event
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
if (usr_value eq 'Wlon_lat_gph') then begin
;
; We need the position for all mouse events
;
	data_crd=convert_coord( $
				event.x, $
				event.y, $
				/device, $
				/to_data)
;
	if ((event.type eq 0) and (event.press eq 1)) then begin
;
; left button was pressed
;
		foo=min(abs(lon-data_crd(0)),lon_idx)
		foo=min(abs(lat-data_crd(1)),lat_idx)
		widget_control,Wlon,set_value=lon(lon_idx)
		widget_control,Wlat,set_value=lat(lat_idx)
		widget_control,Wvalue,set_value=data(lon_idx,lat_idx)
	endif
	if ((event.type eq 0) and (event.press eq 2)) then begin
;
; middle button was pressed
;
		ptr_psn(0)=data_crd(0)
		ptr_psn(1)=data_crd(1)
	endif
	if ((event.type eq 1) and (event.release eq 2)) then begin
;
; middle button was released
;
		ptr_psn(2)=data_crd(0)
		ptr_psn(3)=data_crd(1)
		lon_min=ptr_psn(0)
		lat_min=ptr_psn(3)
		lon_max=ptr_psn(2)
		lat_max=ptr_psn(1)
		widget_control,Wlat_min,set_value=lat_min
		widget_control,Wlat_max,set_value=lat_max
		widget_control,Wlon_min,set_value=lon_min
		widget_control,Wlon_max,set_value=lon_max
	endif
endif; endif Wlon_lat_gph

if command_sng ne '' then print,command_sng

end; end lon_lat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End lon_lat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin make_lon_lat() 
; If the current device is X then this routine will create the 
; widgets and install the event handler for the Latitude Level graphs if
; they don't already exist.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_lon_lat

@ibp_cmn.com

if !d.name eq 'X' then begin

; If the window is already up then OK else start it up.
; be careful because old widget ID's stay in common memory 
; after crashes!

widget_xst=1
if n_elements(Wlon_lat_gph) ne 0 then begin
	widget_xst=widget_info(Wlon_lat_gph,/valid_id)
endif

if ((widget_xst) and (n_elements(Wlon_lat_gph) ne 0)) then begin
curr_sz=indgen(2)
widget_control,Wlon_lat_base,tlb_get_size=curr_sz

; If the window has been resized, then destroy the widget hierarchy
; so that the widget can be recreated at the new size.

wset,lon_lat_wdw

if ((lon_lat_wdw_x_sz ne curr_sz(0)) or (lon_lat_wdw_y_sz ne curr_sz(1))) then begin
	lon_lat_wdw_x_sz=curr_sz(0)
	lon_lat_wdw_y_sz=curr_sz(1)
	widget_control,Wlon_lat_base,/destroy
	widget_xst=0
endif

endif

if ((widget_xst ne 1) or (n_elements(Wlon_lat_gph) eq 0)) then begin

; Create the Latitude Level widget hierarchy

Wlon_lat_base=widget_base( $
	group_leader=Wtop_lvl, $
	title='Longitude Latitude')
Wlon_lat_gph=widget_draw( $
	Wlon_lat_base, $
		xsize=lon_lat_wdw_x_sz, $
		ysize=lon_lat_wdw_y_sz, $
	/button_events, $
	uvalue='Wlon_lat_gph')

widget_control,Wlon_lat_base,/realize
widget_control,get_value=lon_lat_wdw,Wlon_lat_gph

; Call the event manager for the new base widgets
xmanager, $
	'make_lon_lat', $
	Wlon_lat_base, $
	group_leader=Wtop_lvl, $
	event_handler='lon_lat_event'		

endif; end creating widgets

widget_control,Wlon_lat_base,iconify=0 ; de-iconify when drawing
wset,lon_lat_wdw

endif; endif device is X

end; end make_lon_lat() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End make_lon_lat() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Contour Overlay
; It's useful to have this routine stored separately in case it desired
; to plot contours of one field over the image of another.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro contour_ovl
@ibp_cmn.com

contour,data, $
	lon, $			;for maps, first dimen is longitude
	lat, $			;for maps, second dimen is latitude
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $	
	c_annotation=cntr_ntt, $	
;	/closed, $				
	/overplot, $		;don't erase image before contouring
	/follow, $		;smooth contours
	max_value=very_big_nbr, $ ;
	levels=cntr_lvl	;specify my own contour levels

end; end cntr_ovl()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Contour Overlay
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure Contour Map
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lon_lat_fll
@ibp_cmn.com

erase

if (strlen(unit_jd(fgr_idx)) le 0) then unit_jd(fgr_idx)=unit_pfx_sng+unit_sng
if (strlen(ttl_jd(fgr_idx)) le 0) then ttl_jd(fgr_idx)=ttl_sng
if (strlen(xttl_jd(fgr_idx)) le 0) then xttl_jd(fgr_idx)=sub_ttl_sng

map_set, $
	0.0, $			;[dgr] Latitude at center of map
;	lat_ctr, $		;[dgr] Latitude at center of map
	lon_ctr, $		;[dgr] Longitude at center of map
	0.0, $			;[dgr] North rotation angle
	limit=map_lmt_dgr, $;[latmin,lonmin,latmax,lonmax] drawing box
	/cylindrical, $		;cylindrical projection
	xmargin=x_map_cbar_mrg, $
	ymargin=y_map_cbar_mrg, $
	/noborder 		;don't draw a border around the grid

; Convert edge position of map to normal coordinates
img_psn_nrm=convert_coord( $
	[map_lmt_dgr(1), $ ;lon_min
	map_lmt_dgr(3)], $ ;lon_max
	[map_lmt_dgr(0), $ ;lat_min
	map_lmt_dgr(2)], $ ;lat_max
	/data, $
	/to_normal)

; Get rid of z information which was returned, convert to x-y vector
plt_rgn_nrm=[ $
	img_psn_nrm(0), $ ;lon_min
	img_psn_nrm(1), $ ;lat_min
	img_psn_nrm(3), $ ;lon_max
	img_psn_nrm(4)] ;lat_max

cbar_psn=[ $
	plt_rgn_nrm(0), $ ;lon_min
	.06, $
	plt_rgn_nrm(2), $ ;lon_max
	.10] 

clr_mk,cbar_clr_nbr,color_order,pll_idx,background_jd(fgr_idx)

tvlct,r_curr,g_curr,b_curr

clr_bar_drw, $
	bar_chr_sz=cbar_chr_sz, $
	bar_fnt=cbar_fnt, $
	bar_idx=cbar_idx, $
	bar_lbl_sz=cbar_lbl_sz, $
	bar_lgn=cbar_lgn, $
	bar_psn=cbar_psn, $
	bar_txt_clr=cbar_txt_clr, $
	bar_unit=cbar_unit, $
	bar_clr_nbr=cbar_clr_nbr

contour, $
	data, $
	lon, $
	lat, $
	level=cntr_lvl, $                 
	xrange= $
		[map_lmt_dgr(1), $ ;lon_min
		map_lmt_dgr(3)], $ ;lon_max
	yrange= $
		[map_lmt_dgr(0), $ ;lat_min
		map_lmt_dgr(2)], $ ;lat_max
	c_color=cntr_fll_idx, $		
;	c_labels=cntr_which_lbl, $	
;	c_annotation=cntr_ntt, $	
	position=plt_rgn_nrm, $
	xstyle=13, $
	ystyle=13, $
	max_value=very_big_nbr, $ ;
;	/closed, $				
	/fill, $				
	/noerase

if cntr_ovl_jd(fgr_idx) then $
contour, $
	data, $
	lon, $			;for maps, first dimen is longitude
	lat, $			;for maps, second dimen is latitude
	level=cntr_lvl, $                 
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $	
	c_annotation=cntr_ntt, $	
	c_charsize=0.83*!p.charsize, $
;	xstyle=13, $
;	ystyle=13, $
	max_value=very_big_nbr, $ ;
;	/closed, $				
	/follow, $				
	/overplot

; Set boundaries for Cylindrical projections
; NB: Until IDL 4.0, these commands were not needed because map_set() set x.crange and y.crange, which axes routines need to do correct labeling
; Without following two lines, x.crange and y.crange will be [0.0,1.0] and axis labels will be wrong
!x.range=[map_lmt_dgr(1),map_lmt_dgr(3)] ; lon_min,lon_max ; Longitude increases on x-axis
!y.range=[map_lmt_dgr(0),map_lmt_dgr(2)] ; lat_min,lat_max ; Latitude increases on y-axis

map_set, $
	0.0, $			;[dgr] Latitude at center of map
;	lat_ctr, $		;[dgr] Latitude at center of map
	lon_ctr, $		;[dgr] Longitude at center of map
	0.0, $			;[dgr] North rotation angle
	latdel=10.0, $		;latitude grid spacing
	londel=10.0, $		;longitude grid spacing
	limit=map_lmt_dgr, $;[latmin,lonmin,latmax,lonmax] drawing box
;	tit=ttl_sng, $
	/continent, $		;plot the continental boudaries
	/cylindrical, $		;cylindrical projection
;	/label, $		;label the parallels and meridians
;	latalign=0.5, $		;.5 centers the labels above/below the text
;	lonalign=0.5, $		;.5 centers the labels next to the text
;	lonlab=[ $		;latitudes at which to place longitude labels
;		map_lmt_dgr(0), $
;		map_lmt_dgr(2)], $
;	latlab=[ $		;longitudes at which to place latitude labels
;		map_lmt_dgr(1), $
;		map_lmt_dgr(3)], $
	xmargin=x_map_cbar_mrg, $
	ymargin=y_map_cbar_mrg, $
	/noborder, $		;don't draw a border around the grid
	/noerase, $		;don't erase the current image before drawing
;	/usa, $			;draw the usa state boundaries 
	/grid			;draw the grid of parallels and meridians

; Employ common axes and title drawing procedures
drw_axz

drw_ttl,ttl_jd(fgr_idx),xttl_jd(fgr_idx)

end; end lon_lat_fll()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Contour Map
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure Image Map
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lon_lat_img
@ibp_cmn.com

erase

map_set, $
	0.0, $			;[dgr] Latitude at center of map
;	lat_ctr, $		;[dgr] Latitude at center of map
	lon_ctr, $		;[dgr] Longitude at center of map
	0.0, $			;[dgr] North rotation angle
	limit=map_lmt_dgr, $;[latmin,lonmin,latmax,lonmax] drawing box
	/cylindrical, $		;cylindrical projection
	xmargin=x_map_cbar_mrg, $
	ymargin=y_map_cbar_mrg, $
	/noborder 		;don't draw a border around the grid

; convert the position of the edges of the map to normal coordinates
img_psn_nrm=convert_coord( $
	[map_lmt_dgr(1), $ ;lon_min
	map_lmt_dgr(3)], $ ;lon_max
	[map_lmt_dgr(0), $ ;lat_min
	map_lmt_dgr(2)], $ ;lat_max
	/data, $
	/to_normal)

;get rid of the z information which was returned, convert to an x-y vector
plt_rgn_nrm=[ $
	img_psn_nrm(0), $ ;lon_min
	img_psn_nrm(1), $ ;lat_min
	img_psn_nrm(3), $ ;lon_max
	img_psn_nrm(4)] ;lat_max

cbar_psn=[ $
	plt_rgn_nrm(0), $ ;lon_min
	.06, $
	plt_rgn_nrm(2), $ ;lon_max
	.10] 

;the following must probably be used for non-cylindrical maps.
;using map_image has the effect of interpolating between data points 
;which results in a smoother image (and larger file) at printout
;data_image=map_image(data,startx,starty,scale=0.02)
;but this works best for cylindrical maps.
;using the raw data produces a grainier image and smaller file.
data_img=data

arr_sz=size(data_img)
x_arr_sz=arr_sz(1)
y_arr_sz=arr_sz(2)

;convert the x-y normal coord. vector to device coords
img_psn_dvc=fix(plt_rgn_nrm*[ $
	!d.x_size, $
	!d.y_size, $
	!d.x_size, $
	!d.y_size])

;get the width and height in device unit
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
	bar_chr_sz=cbar_chr_sz, $
	bar_fnt=cbar_fnt, $
	bar_idx=cbar_idx, $
	bar_lbl_sz=cbar_lbl_sz, $
	bar_lgn=cbar_lgn, $
	bar_psn=cbar_psn, $
	bar_txt_clr=cbar_txt_clr, $
	bar_unit=cbar_unit, $
	bar_clr_nbr=cbar_clr_nbr

; Set boundaries for Cylindrical projections
; NB: Until IDL 4.0, this commands was not needed because map_set() set x.crange and y.crange, which the axes routines need to do correct labeling
; Without the following two lines, x.crange and y.crange will be [0.0,1.0] and axis labels will be wrong
!x.range=[map_lmt_dgr(1),map_lmt_dgr(3)] ; lon_min,lon_max ; Longitude increases on x-axis
!y.range=[map_lmt_dgr(0),map_lmt_dgr(2)] ; lat_min,lat_max ; Latitude increases on y-axis

map_set, $
	0.0, $			;[dgr] Latitude at center of map
;	lat_ctr, $		;[dgr] Latitude at center of map
	lon_ctr, $		;[dgr] Longitude at center of map
	0.0, $			;[dgr] North rotation angle
	latdel=10.0, $		;latitude grid spacing
	londel=10.0, $		;longitude grid spacing
	limit=map_lmt_dgr, $;[latmin,lonmin,latmax,lonmax] drawing box
;	tit=ttl_sng, $
	/continent, $		;plot the continental boudaries
	/cylindrical, $		;cylindrical projection
;	/label, $		;label the parallels and meridians
;	latalign=0.5, $		;.5 centers the labels above/below the text
;	lonalign=0.5, $		;.5 centers the labels next to the text
;	lonlab=[ $		;latitudes at which to place longitude labels
;		map_lmt_dgr(0), $
;		map_lmt_dgr(2)], $
;	latlab=[ $		;longitudes at which to place latitude labels
;		map_lmt_dgr(1), $
;		map_lmt_dgr(3)], $
	xmargin=x_map_cbar_mrg, $
	ymargin=y_map_cbar_mrg, $
	/noborder, $		;don't draw a border around the grid
;	/usa, $			;draw the usa state boundaries 
;	/grid, $		;draw the grid of parallels and meridians
	/noerase 		;don't erase the current image before drawing

if cntr_ovl_jd(fgr_idx) then $
contour,data, $
	lon, $			;for maps, first dimen is longitude
	lat, $			;for maps, second dimen is latitude
	levels=cntr_lvl, $	;specify my own contour levels
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $	
	c_annotation=cntr_ntt, $	
	c_charsize=0.83*!p.charsize, $
	max_value=very_big_nbr, $ ;
;	/closed, $				
	/overplot, $		;don't erase image before contouring
	/follow			;smooth contours

; employ the common axes and title drawing procedures
drw_axz

drw_ttl,ttl_jd(fgr_idx),xttl_jd(fgr_idx)

end; end lon_lat_img()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Image Map 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure Polar Map
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro pol_img
@ibp_cmn.com

erase

min_lat=min(lat)
max_lat=max(lat)
min_lon=min(lon)
max_lon=max(lon)

if (max_lat gt abs(min_lat)) then begin
	middle_lat=lat_max
	lon_lbl_lat=lat_min-2.5
	lat_lbl_lon=0.
endif else begin
	middle_lat=lat_min
	lon_lbl_lat=lat_max+2.5
	lat_lbl_lon=0.
endelse

map_set, $
	middle_lat, $		;[dgr] Latitude at center of map
	lon_ctr, $		;[dgr] Longitude at center of map
	lon_ctr, $		;[dgr] Longitude at center of map
	limit=map_lmt_dgr, $;[latmin,lonmin,latmax,lonmax] drawing box
	/azi, $			;azimuthal projection
	xmargin=[1.0,1.0], $	;[10,3] is [left,right] default
	ymargin=[12.0,10.0], $	;[4,2] is [bottom,top] default
	/noborder 		;don't draw a border around the grid

cbar_psn=[0.01,0.06,0.99,0.10] 

print,'middle_lat = ',middle_lat
print,'lat_ctr = ',lat_ctr
print,'lon_ctr = ',lon_ctr
print,'map_lmt_dgr = ',map_lmt_dgr
print,'size(data) = ',size(data)
print,'lat_min = ',lat_min
print,'lat_max = ',lat_max
print,'lon_min = ',lon_min
print,'lon_max = ',lon_max
print,'min_lat = ',min_lat
print,'max_lat = ',max_lat
print,'min_lon = ',min_lon
print,'max_lon = ',max_lon
;print,'img_psn_nrm = ',img_psn_nrm
;print,'plt_rgn_nrm = ',plt_rgn_nrm
;print,'lon_lat_wdw_x_sz = ',lon_lat_wdw_x_sz
;print,'lon_lat_wdw_y_sz = ',lon_lat_wdw_y_sz
;print,'img_psn_dvc = ',img_psn_dvc
;print,'dvc_wdt = ',dvc_wdt
;print,'dvc_hgt = ',dvc_hgt

; When scaling images it is important to place the top argument at
; no greater than !d.table_size-1 since there are at most 256 colors
; with PostScript and !d.table_size=256. Since the minimum value is
; always 0, then the colors that will show on the screen/paper are
; color indices 0..255.

; See User's Guide p. 19-15

; here follows and attempt to add an extra longitude column

;	img_sz=size(data)
;	img_data=fltarr(img_sz(1)+2,img_sz(2))
;	img_data(1,0)=data
;	img_data(0,0)=img_data(0,*)
;	img_data(img_sz(1)+1,0)=img_data(0,*)
;	img_lon=[lon_min,lon,lon_max]

	data_scl=bytscl(data,min=cntr_lvl_min,max=cntr_lvl_max, $
;min=data_min, $ ; in order to use data_min, must tie lowest contour to data_min, i.e., cntr_lvl_min must equal data_min
;max=data_max, $ ; in order to use data_max, must tie lowest contour to data_min, i.e., cntr_lvl_max must equal data_max
	                top=(!d.table_size-3))+2

	data_scl_map=map_image(data_scl,startx,starty,xsize,ysize, $
;			latmin=lat_min, $
;			latmax=lat_max, $
;			lonmin=lon_min, $
;			lonmax=lon_max, $
			latmin=min_lat, $
			latmax=max_lat, $
			lonmin=min_lon, $
			lonmax=max_lon, $
			/bilin)

if(!d.name eq 'PS') then tv,data_scl_map,startx,starty,xsize=xsize,ysize=ysize else tv,data_scl_map,startx,starty

clr_bar_drw, $
	bar_chr_sz=cbar_chr_sz, $
	bar_fnt=cbar_fnt, $
	bar_idx=cbar_idx, $
	bar_lbl_sz=cbar_lbl_sz, $
	bar_lgn=cbar_lgn, $
	bar_psn=cbar_psn, $
	bar_txt_clr=cbar_txt_clr, $
	bar_unit=cbar_unit, $
	bar_clr_nbr=cbar_clr_nbr

if !d.name ne 'X' then chr_sz=1.25 else chr_sz=1.0
map_set, $
	middle_lat, $		;[dgr] Latitude at center of map
	lon_ctr, $		;[dgr] Longitude at center of map
	lon_ctr, $		;rotation
	latdel=15., $		;latitude grid spacing
	londel=15., $		;longitude grid spacing
	limit=map_lmt_dgr, $;[latmin,lonmin,latmax,lonmax] drawing box
;	tit=ttl_sng, $
	/continent, $		;plot the continental boudaries
	/azi, $			;azimuthal
	/label, $		;label the parallels and meridians
	latalign=0.5, $		;.5 centers the labels above/below the text
	lonalign=0.5, $		;.5 centers the labels next to the text
	lonlab=lon_lbl_lat, $	;latitudes at which to place longitude labels
	latlab=lat_lbl_lon, $	;longitudes at which to place latitude labels
	xmargin=[1.0,1.0], $	;[10,3] is [left,right] default
	ymargin=[12.0,10.0], $	;[4,2] is [bottom,top] default
	charsize=chr_sz, $
	/noborder, $		;don't draw a border around the grid
	/noerase, $		;don't erase the current image before drawing
;	/usa, $			;draw the usa state boundaries 
	/grid			;draw the grid of parallels and meridians

if cntr_ovl_jd(fgr_idx) then $
contour,data, $
	lon, $			;for maps, first dimen is longitude
	lat, $			;for maps, second dimen is latitude
	levels=cntr_lvl, $	;specify my own contour levels
	c_thick=cntr_thk, $ 		;line_thk
	c_linestyle=cntr_ln_sty, $	;line_styles
	c_labels=cntr_which_lbl, $	
	c_annotation=cntr_ntt, $	
	max_value=very_big_nbr, $ ;
;	/closed, $				
	/overplot, $		;don't erase image before contouring
	/follow			;smooth contours

if !d.name ne 'X' then ttl_sz=2.0 else ttl_sz=1.25
if !d.name ne 'X' then sub_ttl_sz=1.5 else sub_ttl_sz=0.75

xyouts,0.5,0.95,ttl_jd(fgr_idx),size=ttl_sz,alignment=0.5,/NORMAL
xyouts,0.5,0.91,xttl_jd(fgr_idx),size=sub_ttl_sz,alignment=0.5,/NORMAL

;drw_ttl,ttl_jd(fgr_idx),xttl_jd(fgr_idx)

end; end pol_img()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure Polar Map 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

