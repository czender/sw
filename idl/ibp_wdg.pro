; $Id$

; Usage: 
; ibp,file='/data/zender/_aux0_/ccm/isccp.jul.nc'
; ibp,file='/fs/cgd/data0/zender/dst06/dst06_01.nc'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin ibp() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ibp_event,event

@ibp_cmn.com

event_type=tag_names(event,/structure)

print,'event_type = ',event_type
print,event_type+' tags names = ',tag_names(event)
print,'event = ',event

; In case event is not handled by the below
command_sng=''

; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value
print,'usr_value = ',usr_value

if (usr_value eq 'Wquit') then begin
	widget_control,/destroy,event.top
	command_sng=usr_value
endif; endif Wquit

if (usr_value eq 'Wdraw') then begin
drw_pretty_picture	
endif; endif Wdraw

if (usr_value eq 'Wfld') then begin
	reset_fld
	fld_idx=event.index
	fld_nm=fld_lst(fld_idx).name
	rfr_dim_lst
	command_sng='fld = '+fld_lst(fld_idx).name
endif; endif Wfld

if (usr_value eq 'Wread') then begin
	decrypt_fl_nm
	read_data
	scale_data
	rgn_trunc
	diagnose_data
	form_cntr_lvl
	comp_cntr_lvl
	bld_sngs
	command_sng='read fld '+fld_lst(fld_idx).name
endif; endif Wread

if (usr_value eq 'Wfl_in') then begin
	widget_control,/hourglass
	pth_nm_get,fl_in,path,name
	new_fl=dialog_pickfile( $
		/read, $
		path=path, $
		file=fl_in, $
		filter='*.nc')
	if (new_fl ne '') then begin
		fl_in=new_fl
		pth_nm_get,fl_in,path,name
		widget_control,Wfl_in,set_value=name
		read_fl
		decrypt_fl_nm
		foo=where(fld_lst.name eq fld_nm,count)
		if count eq 0 then fld_idx=0 else fld_idx=foo(0)
		if(fld_lst(fld_idx).dim_nbr lt 2) then fld_idx=0
		fld_nm=fld_lst(fld_idx).name
		widget_control,Wfld,set_value=fld_lst.name
		widget_control,Wfld,set_list_select=fld_idx
		rfr_dim_lst
	endif
	command_sng='dialog_pickfile() = '+fl_in
endif; endif Wfl_in

if (usr_value eq 'Wfl_out') then begin
if !d.name eq 'X' then begin
	widget_control,/hourglass
	pth_nm_get,fl_out,path,name
	new_fl=dialog_pickfile( $
		title='Select output file', $
		path=path, $
		file=fl_out, $
		filter='*ps')
	if (new_fl ne '') then begin
		fl_out=new_fl
		pth_nm_get,fl_out,path,name
		widget_control,Wfl_out,set_value=name
	endif
	command_sng='dialog_pickfile() = '+fl_out
endif else begin
	command_sng='Changing output filenames is only allowed when device is X windows'
endelse
endif; endif Wfl_out

if (usr_value eq 'Wpalette') then begin
	pll_idx=event.index
	command_sng='pll_idx = '+pll_idx+' = '+pll_lst(pll_idx).english
endif; endif Wpalette

if (usr_value eq 'Wcolor_order') then begin
	color_order=(color_order+1) mod 2
	clr_rfr,clr_tbl,tbl_fl,usr_dfn_clr_tbl
endif; end Wcolor_order

if (usr_value eq 'Wusr_hook') then begin
	if event.value eq 'Strong ocean precip.' then begin
		precc_cmfmc	
	endif else if event.value eq 'Cloud clusters' then begin
		precc_cmfmc
	endif
	command_sng=event.value
endif; endif Wusr_hook

if (usr_value eq 'Wcstm') then begin
	command_sng='bringing up customization screen'

; Popup the context-appropriate widget screen
	make_cstm	
endif; endif Wcstm

if (usr_value eq 'Wtbl_fl') then begin
        widget_control,/hourglass
        pth_nm_get,tbl_fl,path,name
        new_fl=dialog_pickfile( $
                /read, $
                path=path, $
                file=tbl_fl, $
                filter='*.tbl')
        if (new_fl ne '') then begin
		usr_dfn_clr_tbl=True
                tbl_fl=new_fl
                pth_nm_get,tbl_fl,path,name
		clr_rfr,clr_tbl,tbl_fl,usr_dfn_clr_tbl
                widget_control,Wtbl_fl,set_value=name
        endif
        command_string='dialog_pickfile() = '+tbl_fl
endif; endif Wtbl_fl

if (usr_value eq 'Wsv') then begin
; Make it so that saving the same data twice unsaves the data.
; This is a simpler method than creating an unsave button.
; Changing sv_data to scalar zero, and zeroing a few
; other should suffice.
; The 
; KLUDGE KLUDGE KLUDGE KLUDGE KLUDGE KLUDGE KLUDGE KLUDGE KLUDGE
; is comparing data by checking their total.

	OK_TO_SV_DATA=True
	if total(sv_data) eq total(data) then begin
		sv_data=0
		sv_data_nbr=0
	        sv_data_max=0
	        sv_data_min=0
		command_sng='unsaving  '+sv_fld_nm
		sv_fld_nm='None Saved'
		OK_TO_SV_DATA=False
	endif; endif the totals of data and sv_data are equal

	if OK_TO_SV_DATA then begin
	bld_sngs
	sv_abb_sng=abb_sng
	sv_avg_data=avg_data
	sv_data=data
	sv_fld_nm=fld_nm
	sv_fld_sng=fld_sng
	sv_lat=lat
	sv_lon=lon
        sv_data_max=data_max
        sv_data_min=data_min
	sv_data_nbr=data_nbr
	sv_lat_nbr=lat_nbr
	sv_lon_nbr=lon_nbr
	sv_slc_sng=slice_sng
	sv_src_sng=src_sng
	sv_sub_ttl_sng=sub_ttl_sng
	sv_sym_sng=sym_sng
	sv_unit_sng=unit_sng
	command_sng='saving '+sv_fld_nm
	endif; endif saving data

	widget_control,Wsv,set_value=sv_fld_nm
endif; endif Wsv

if (usr_value eq 'Wdifference') then begin
	difference
	command_sng='differencing '+fld_nm+' - '+sv_fld_nm
endif; endif Wdifference

if (usr_value eq 'Wmask') then begin
	msk_mode=True
	msk_data=data
	msk_fld_nm=fld_nm
	widget_control,Wmask,set_value='Mask: '+msk_fld_nm
	command_sng='mask field = '+msk_fld_nm
endif; endif Wmask

if (usr_value eq 'Wrgn') then begin
	rgn_idx=event.index
	lon_min=rgn_lst(rgn_idx).lon_min
	lat_min=rgn_lst(rgn_idx).lat_min
	lon_max=rgn_lst(rgn_idx).lon_max
	lat_max=rgn_lst(rgn_idx).lat_max
	if lon_min gt lon_max then lon_ctr=180.0-0.5*(lon_min+lon_max)
	if lon_min lt lon_max then lon_ctr=0.5*(lon_min+lon_max)
	lat_ctr=0.5*(lat_min+lat_max)
	map_lmt_dgr=[ $ ;[latmin,lonmin,latmax,lonmax]
			lat_min, $
			lon_min, $
			lat_max, $
			lon_max] 
	widget_control,Wlat_min,set_value=lat_min
	widget_control,Wlat_max,set_value=lat_max
	widget_control,Wlon_min,set_value=lon_min
	widget_control,Wlon_max,set_value=lon_max
	command_sng='region = '+rgn_lst(rgn_idx).english
endif; endif Wrgn

;if (usr_value eq 'Woptions') then begin
;	if event.value eq 'Hi/Lo Labels' then begin
;		hi_lo_plt=not hi_lo_plt
;		high_low_lbls		
;		command_sng='hi_lo_plt = '+auto_sng(hi_lo_plt,0)
;	endif else if event.value eq 'Mask Mode' then begin
;		msk_mode=not msk_mode
;		command_sng='msk_mode = '+auto_sng(msk_mode,0)
;	endif
;endif; endif Woptions

if (usr_value eq 'Wx_dim') then begin
	x_dim=event.index
	command_sng='x_dim = event.index = '+auto_sng(x_dim,0)
endif; endif Wx_dim

if (usr_value eq 'Wy_dim') then begin
	y_dim=event.index
	command_sng='y_dim = event.index = '+auto_sng(y_dim,0)
endif; endif Wy_dim

if (usr_value eq 'Wtime_slc') then begin
	time_slc=event.index
	command_sng='time_slc = event.index = '+auto_sng(time_slc,0)
endif; endif Wtime_slc

if (usr_value eq 'Wvrt_slc') then begin
	vrt_slc=event.index
	command_sng='vrt_slc = event.index = '+auto_sng(vrt_slc,0)
endif; endif Wvrt_slc

if (usr_value eq 'Wfgr') then begin
	if (fgr_lst(event.index).fnc_nm ne '') then fgr_idx=event.index
	command_sng='fgr = '+fgr_lst(fgr_idx).fnc_nm
endif; endif Wfgr

if (usr_value eq 'Wformat') then begin
args=''
fmt_idx=event.index
case fmt_lst(fmt_idx).english of
	'X': begin
		close_ps,fl_nm=fl_out
	end ; end X
	'EPS': begin
; if fl_out name suffix ends in .ps then change it .eps
;		result=strpos(fl_out,'.ps')
;		if result ne -1 then begin
;			fl_out=strmid(fl_out,0,strlen(fl_out)-2)+'eps'
;			pth_nm_get,fl_out,path,name
;			widget_control,Wfl_out,set_value=name
;		endif
		args=args+fmt_lst(fmt_idx).arg
		args=args+chroma_lst(chroma_idx).arg
		args=args+shape_lst(shape_idx).arg
		args=args+',fl_nm='''+fl_out+'''
		args='open_ps'+args
		result=execute(args)
		command_sng=args
	end ; end EPS
	'GIF': begin
	; 20011012: Using GIF now requires license from Unisys!
		image=tvrd()
		pth_nm_get,fl_out,path,name
		raster_fl=path+'/'+fld_nm+'.gif'
		print,'fl_out = ',fl_out
		print,'path = ',path
		print,'name = ',name
		print,'raster_fl = ',raster_fl
		write_gif,raster_fl,image
		fmt_idx=0
		widget_control,Wformat,set_list_select=fmt_idx
		command_sng='Wrote currently displayed image to GIF format file '+raster_fl
	end ; end GIF
	'JPEG': begin
		; true=[1,2,3] required for truecolor images, specifies dimension for color data
		image=tvrd(true=1)
		pth_nm_get,fl_out,path,name
		raster_fl=path+'/'+fld_nm+'.jpg'
		write_jpeg,raster_fl,image,true=1
		fmt_idx=0
		widget_control,Wformat,set_list_select=fmt_idx
		command_sng='Wrote currently displayed image to JPEG/JFIF format file '+raster_fl
	end ; end JPEG
	'PNG': begin
		; true=[1,2,3] required for truecolor images, specifies dimension for color data
		image=tvrd(true=1)
		pth_nm_get,fl_out,path,name
		raster_fl=path+'/'+fld_nm+'.png'
		write_png,raster_fl,image
		fmt_idx=0
		widget_control,Wformat,set_list_select=fmt_idx
		command_sng='Wrote currently displayed image to PNG format file '+raster_fl
	end ; end PNG
	'PPM': begin
		; Warning: write_ppm() only works on 2D rasters, so 3-channel images (True Color displays) do not work with it
		image=tvrd()
		pth_nm_get,fl_out,path,name
		raster_fl=path+'/'+fld_nm+'.ppm'
		write_ppm,raster_fl,image
		fmt_idx=0
		widget_control,Wformat,set_list_select=fmt_idx
		command_sng='Wrote currently displayed image to PPM format file '+raster_fl
	end ; end PPM
	'Postscript': begin
; if fl_out name suffix ends in .eps then change it .ps
;		result=strpos(fl_out,'.eps')
;		if result ne -1 then begin
;			fl_out=strmid(fl_out,0,strlen(fl_out)-3)+'ps'
;			pth_nm_get,fl_out,path,name
;			widget_control,Wfl_out,set_value=name
;		endif
		args=args+chroma_lst(chroma_idx).arg
		args=args+shape_lst(shape_idx).arg
		args=args+',fl_nm='''+fl_out+'''
		args='open_ps'+args
		result=execute(args)
		command_sng=args
	end ; end Postscript
	else:
endcase
endif; endif Wformat

if (usr_value eq 'Wdest') then begin
args=''
dest_idx=event.index
args=args+dest_lst(dest_idx).arg
args='close_ps'+args+',fl_nm='''+fl_out+'''
result=execute(args)
fmt_idx=0
widget_control,Wformat,set_list_select=fmt_idx
widget_control,Wdest,set_list_select=-1
command_sng=args
endif; endif Wdest

if (usr_value eq 'Wshape') then begin
shape_idx=event.index
command_sng='shape = '+shape_lst(shape_idx).english
endif; endif Wshape

if (usr_value eq 'Wchroma') then begin
chroma_idx=event.index
command_sng='chroma = '+chroma_lst(chroma_idx).english
endif; endif Wchroma

if (usr_value eq 'Wclr_tbl') then begin
	usr_dfn_clr_tbl=False
	clr_tbl=event.value
	if clr_tbl ge 40 then begin
		clr_tbl=22
		widget_control,Wclr_tbl,set_value=clr_tbl
	endif
	clr_rfr,clr_tbl,tbl_fl,usr_dfn_clr_tbl
	command_sng='clr_tbl = event.value = '+string(event.value)
endif; endif Wclr_tbl

if (usr_value eq 'Wlat_min') then begin
	lat_min=event.value
	foo=min(abs(lat-lat_max),lat_max_idx)
	widget_control,Wrgn,set_list_select=-1
	command_sng='lat_min = event.value = '+string(event.value)
endif; endif Wlat_min

if (usr_value eq 'Wlat_max') then begin
	lat_max=event.value
	foo=min(abs(lat-lat_max),lat_max_idx)
	widget_control,Wrgn,set_list_select=-1
	command_sng='lat_max = event.value = '+string(event.value)
endif; endif Wlat_max

if (usr_value eq 'Wlon_min') then begin
	lon_min=event.value
	foo=min(abs(lon-lon_min),lon_min_idx)
	widget_control,Wrgn,set_list_select=-1
	command_sng='lon_min = event.value = '+string(event.value)
endif; endif Wlon_min

if (usr_value eq 'Wlon_max') then begin
	lon_max=event.value
	foo=min(abs(lon-lon_max),lon_max_idx)
	widget_control,Wrgn,set_list_select=-1
	command_sng='lon_max = event.value = '+string(event.value)
endif; endif Wlon_max

if (usr_value eq 'Wlev_min') then begin
	lev_min=event.value
	foo=min(abs(lev-lev_min),lev_min_idx)
	command_sng='lev_min = event.value = '+string(event.value)
endif; endif Wlev_min

if (usr_value eq 'Wlev_max') then begin
	lev_max=event.value
	foo=min(abs(lev-lev_max),lev_max_idx)
	command_sng='lev_max = event.value = '+string(event.value)
endif; endif Wlev_max

if (usr_value eq 'Wanim') then begin
	command_sng='firing up animation widgets'

; Popup a widget to allow animation settings.
make_anim

endif; endif Wanim

if (usr_value eq 'Wbreak') then begin
; intentionally crash the program for easier restarts
	foo=non_existant_array(-1)
endif; endif Wbreak

if command_sng ne '' then print,command_sng

end; end ibp() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End ibp() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin main program ibp() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin main() Initialization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ibp,file=file

; ibp,file='/fs/cgd/csm/input/atm/SOx_emissions.nc'

@ibp_cmn.com

CVS_Revision='$Revision$'

;Perform initial housekeeping
if !d.name ne 'X' then device,/close_file
set_plot,'X'
close,/all

;Initialize internal defaults
True=1
False=0

; Default window sizes should be multiples of 16 for easier MPEG encoding
color_order=1
clr_tbl=22
encapsulate=False
hi_lo_plt=False
hst_wdw_x_sz=512
hst_wdw_y_sz=512
lat_lev_wdw_x_sz=512
lat_lev_wdw_y_sz=512
lon_lev_wdw_x_sz=512
lon_lev_wdw_y_sz=512
lon_lat_wdw_x_sz=768
lon_lat_wdw_y_sz=512
msk_mode=False
msk_val=0.
anim_nbr=18
fl_out='./idl.ps'
pll_idx=4
sv_data=0
sv_data_nbr=0
scat_wdw_x_sz=512
scat_wdw_y_sz=512
time_lat_wdw_x_sz=512
time_lat_wdw_y_sz=512
time_slc=0
usr_dfn_clr_tbl=False
vld_avg_thr_pct=1.
vrt_slc=7
very_big_nbr=1.0e20
znl_avg_wdw_x_sz=512
znl_avg_wdw_y_sz=512

; make space for mouse event histories
ptr_psn=fltarr(4)

; Provide enough room to allow printing with or without the colorbar, etc.
x_map_mrg=[7,7]		;[10,3] is [left,right] default
y_map_mrg=[7,7]		;[4,2] is [bottom,top] default

; The following two fill the window and let the device offsets do the scaling
if !d.name ne 'X' then begin
	x_map_cbar_mrg=[7.5,3]	;[10,3] is [left,right] default
	y_map_cbar_mrg=[11,7]	;[4,2] is [bottom,top] default
endif else begin
	x_map_cbar_mrg=[7.5,3]	;[10,3] is [left,right] default
	y_map_cbar_mrg=[8,5]		;[4,2] is [bottom,top] default
endelse

data_sng=getenv('DATA')
usr_sng=getenv('LOGNAME')
hostname_sng=getenv('HOST')
pwd_sng=getenv('PWD')
dbg_txt_sz=0.65
dbg=0
dbg_vrb=True

if usr_sng ne 'zender' then begin
	if n_elements(file) eq 0 then fl_in='/data/zender/data/dst_T42.nc' else fl_in=file
endif else begin
	if n_elements(file) eq 0 then fl_in=data_sng+'/data/dst_T42.nc' else fl_in=file
endelse
tbl_fl=pwd_sng+'/default.tbl'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End main() Initialization
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Data Structures
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@ibp_sct.pro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Data Structures
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Pre-Widget Initialize
; Some data is essential before mapping widgets
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
read_fl
decrypt_fl_nm
reset_fld

if n_elements(rgn_idx) eq 0 then rgn_idx=rgn_nbr-1
lon_min=rgn_lst(rgn_idx).lon_min
lat_min=rgn_lst(rgn_idx).lat_min
lon_max=rgn_lst(rgn_idx).lon_max
lat_max=rgn_lst(rgn_idx).lat_max
lev_min=100.0
lev_max=1000.0
map_lmt_dgr=[ $ ;[latmin,lonmin,latmax,lonmax]
		lat_min, $
		lon_min, $
		lat_max, $
		lon_max] 
if lon_min gt lon_max then lon_ctr=180.-.5*(lon_min+lon_max)
if lon_min lt lon_max then lon_ctr=0.5*(lon_min+lon_max)
lat_ctr=0.5*(lat_min+lat_max)

bin_sz=1.
chroma_idx=1
cntr_lvl=0.
dest_idx=2
device=16
fld_idx=0
fld_nm=fld_lst(fld_idx).name
fgr_idx=0
fmt_idx=0
bin_max=0.
bin_min=0.
shape_idx=1

pth_nm_get,fl_in,path,name
pth_nm_get,fl_out,foo_path,out_nm
pth_nm_get,tbl_fl,foo_path,tbl_nm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Pre-Widget Initialize
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Widgets
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Set up base widgets
Wtop_lvl=widget_base(title='IBP '+strmid(CVS_Revision,11,3),/row,resource_name='ibp')
Wfirst_lft_clm=widget_base(Wtop_lvl,/frame,/column,resource_name='first_lft_clm')
Wfirst_rgt_clm=widget_base(Wtop_lvl,/frame,/column,resource_name='first_rgt_clm')
Wsecond_rgt_clm=widget_base(Wtop_lvl,/frame,/column,resource_name='second_rgt_clm')

Wquit=widget_button( $
		Wfirst_lft_clm, $
		value='Quit', $
		resource_name='quit', $
		uvalue='Wquit')

Wdraw=widget_button( $
		Wfirst_lft_clm, $
		value='Draw', $
		resource_name='draw', $
		uvalue='Wdraw')

Wcstm=widget_button( $
		Wfirst_lft_clm, $
		value='Customize plots', $
		resource_name='custom',$
		uvalue='Wcstm')

Wbreak=widget_button( $
		Wfirst_lft_clm, $
		value='Break', $
		resource_name='break', $
		uvalue='Wbreak')

Wusr_hook=cw_pdmenu( $
		Wfirst_lft_clm, $
		[{cw_pdmenu_s,1,'User Hook'}, $
			{cw_pdmenu_s,0,'Strong ocean precip.'}, $
			{cw_pdmenu_s,2,'Cloud clusters'}], $
		/return_name, $
;		resource_name='usr_hook', $
		uvalue='Wusr_hook')

Wscl=cw_field( $
		Wfirst_lft_clm, $
		title='Scale', $
		value=scl, $
		xsize=13, $
		/row, $
;		resource_name='scl', $
		/noedit)

Wdata_min=cw_field( $
		Wfirst_lft_clm, $
		title='Min', $
		value=data_min, $
		xsize=13, $
		/row, $
;		resource_name='data_min', $
		/noedit)

Wdata_max=cw_field( $
		Wfirst_lft_clm, $
		title='Max', $
		value=data_max, $
		xsize=13, $
		/row, $
;		resource_name='data_max', $
		/noedit)

Wavg_data=cw_field( $
		Wfirst_lft_clm, $
		title='Avg', $
		value=avg_data, $
		xsize=13, $
		/row, $
;		resource_name='avg_data', $
		/noedit)

;Woptions=cw_pdmenu( $
;		Wfirst_lft_clm, $
;		[{cw_pdmenu_s,1,'Options'}, $
;			{cw_pdmenu_s,0,'Hi/Lo Labels'}, $
;			{cw_pdmenu_s,2,'Mask Mode'}], $
;		/return_name, $
;;		resource_name='options', $
;		uvalue='Woptions')

Wvalue=cw_field( $
		Wfirst_lft_clm, $
		title='Val', $
		value=0.0, $
		xsize=13, $
		/row, $
;		resource_name='value', $
		/noedit)

Wlat=cw_field( $
		Wfirst_lft_clm, $
		title='Lat', $
		value=0.0, $
		xsize=5, $
		/row, $
;		resource_name='lat', $
		/noedit)

Wlon=cw_field( $
		Wfirst_lft_clm, $
		title='Lon', $
		value=0.0, $
		xsize=5, $
		/row, $
;		resource_name='lon', $
		/noedit)

Wlev_min=cw_field( $
		Wfirst_lft_clm, $
		title='Min Lev ', $
		value=lev_min, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
;		resource_name='lev_min', $
		uvalue='Wlev_min')

Wlev_max=cw_field( $
		Wfirst_lft_clm, $
		title='Max Lev ', $
		value=lev_max, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
;		resource_name='lev_max', $
		uvalue='Wlev_max')

Wlat_min=cw_field( $
		Wfirst_lft_clm, $
		title='Min Lat ', $
		value=lat_min, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
;		resource_name='lat_min', $
		uvalue='Wlat_min')

Wlat_max=cw_field( $
		Wfirst_lft_clm, $
		title='Max Lat ', $
		value=lat_max, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
;		resource_name='lat_max', $
		uvalue='Wlat_max')

Wlon_min=cw_field( $
		Wfirst_lft_clm, $
		title='Min Lon', $
		value=lon_min, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
;		resource_name='lon_min', $
		uvalue='Wlon_min')

Wlon_max=cw_field( $
		Wfirst_lft_clm, $
		title='Max Lon ', $
		value=lon_max, $
		xsize=11, $
		/all_events, $
		/row, $
		/float, $
;		resource_name='lon_max', $
		uvalue='Wlon_max')

Wanim=widget_button( $
		Wfirst_lft_clm, $
		value='Animate', $
;		resource_name='anim', $
		uvalue='Wanim')

Wfl_in=widget_button( $
		Wfirst_rgt_clm, $
		value=name, $
		resource_name='fl_in', $
		uvalue='Wfl_in')

Wsv=widget_button( $
		Wfirst_rgt_clm, $
		value='None Saved', $
		resource_name='save', $
		uvalue='Wsv')

Wdifference=widget_button( $
		Wfirst_rgt_clm, $
		value='Difference', $
		resource_name='difference', $
		uvalue='Wdifference')

Wmask=widget_button( $
		Wfirst_rgt_clm, $
		value='No Mask', $
		resource_name='mask', $
		uvalue='Wmask')

Wfld=widget_list( $
		Wfirst_rgt_clm, $
		value=fld_lst.name, $
		ysize=5, $
		resource_name='fld', $
		uvalue='Wfld') 

Wread=widget_button( $
		Wfirst_rgt_clm, $
		value='Read', $
		resource_name='read', $
		uvalue='Wread')

Wrgn=widget_list( $
		Wfirst_rgt_clm, $
		value=rgn_lst.english, $		
		ysize=5, $
		resource_name='rgn', $
		uvalue='Wrgn')

Wdim_base=widget_base( $
		Wfirst_rgt_clm, $
		/row, $
;		resource_name='dim_base', $
		uvalue='Wdim_base')

Wtime_slc_base=widget_base( $
		Wfirst_rgt_clm, $
;		resource_name='time_slc_base', $
		uvalue='Wtime_slc_base')

Wlat_list_base=widget_base( $
		Wfirst_rgt_clm, $
;		resource_name='lat_list_base', $
		uvalue='Wlat_list_base')

Wvrt_slc_base=widget_base( $
		Wfirst_rgt_clm, $
;		resource_name='vrt_slc_base', $
		uvalue='Wvrt_slc_base')

Wlon_list_base=widget_base( $
		Wfirst_rgt_clm, $
;		resource_name='lon_list_base', $
		uvalue='Wlon_list_base')

Wfl_out=widget_button( $
		Wsecond_rgt_clm, $
		value=out_nm, $
		resource_name='fl_out', $
		uvalue='Wfl_out')

Wfgr=widget_list( $
		Wsecond_rgt_clm, $
		value=fgr_lst.english, $
		ysize=fgr_nbr, $
		resource_name='fgr', $
		uvalue='Wfgr') 

Wshape=widget_list( $
		Wsecond_rgt_clm, $
		value=shape_lst.english, $
		ysize=shape_nbr, $
		resource_name='shape', $
		uvalue='Wshape') 

Wpalette=widget_list( $
		Wsecond_rgt_clm, $
		value=pll_lst.english, $
		ysize=5, $
		resource_name='palette', $
		uvalue='Wpalette')

Wcolor_order_base=widget_base($
		Wsecond_rgt_clm, $
		/exclusive, $
		/row)

Wcolor_order=widget_button( $
		Wcolor_order_base, $
		value='Order', $
		uvalue='Wcolor_order')

Wclr_tbl=cw_field( $
		Wsecond_rgt_clm, $
		title='Color', $
		value=clr_tbl, $
		xsize=2, $
		/all_events, $
		/integer, $
		/row, $
;		resource_name='clr_tbl', $
		uvalue='Wclr_tbl')

Wtbl_fl=widget_button( $
                Wsecond_rgt_clm, $
                value=tbl_nm, $
                resource_name='tbl_fl', $
                uvalue='Wtbl_fl')

Wchroma=widget_list( $
		Wsecond_rgt_clm, $
		value=chroma_lst.english, $
		ysize=chroma_nbr, $
		resource_name='chroma', $
		uvalue='Wchroma') 

Wformat=widget_list( $
		Wsecond_rgt_clm, $
		value=fmt_lst.english, $
		ysize=fmt_nbr, $
		resource_name='format', $
		uvalue='Wformat') 

Wdest=widget_list( $
		Wsecond_rgt_clm, $
		value=dest_lst.english, $
		ysize=dest_nbr, $
		resource_name='dest', $
		uvalue='Wdest') 

widget_control,Wcolor_order,set_button=color_order
widget_control,Wfld,set_list_select=fld_idx
;widget_control,Wdest,set_list_select=dest_idx
widget_control,Wformat,set_list_select=fmt_idx
widget_control,Wshape,set_list_select=shape_idx
widget_control,Wchroma,set_list_select=chroma_idx
widget_control,Wfgr,set_list_select=fgr_idx
widget_control,Wpalette,set_list_select=pll_idx
Widget_control,Wrgn,set_list_select=rgn_idx
widget_control,Wtop_lvl,/realize
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Widgets
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Post-Widget Initialize
; Now that widget IDs are available, load default field
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
rfr_dim_lst
read_data
scale_data
rgn_trunc
diagnose_data
form_cntr_lvl
comp_cntr_lvl
bld_sngs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Post-Widget Initialize
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Putting up a default graph should be done here
clr_rfr,clr_tbl,tbl_fl,usr_dfn_clr_tbl

; Call event manager
xmanager, $
	'ibp', $
	Wtop_lvl, $
	event_handler='ibp_event'

end; end ibp()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End main program ibp() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


