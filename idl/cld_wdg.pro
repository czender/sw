;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RCS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: cld_wdg.pro,v $
; $Source: /home/zender/cvs/idl/cld_wdg.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; Purpose: All the IDL figures for the shape and size sensitivity paper.
;
; $Log: not supported by cvs2svn $
; Revision 1.9  2000-03-01 02:28:33  zender
; *** empty log message ***
;
; Revision 1.8  2000/01/24 20:03:32  zender
; *** empty log message ***
;
; Revision 1.4  2000/01/01 01:55:48  zender
; *** empty log message ***
;
; Revision 1.3  1999/12/31 20:12:46  zender
; *** empty log message ***
;
; Revision 1.2  1999/12/31 02:09:36  zender
; *** empty log message ***
;
; Revision 1.1  1999/12/31 00:18:17  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/04 23:37:12  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:46  zender
; Imported sources
;
; Revision 2.0  1994/03/08  00:54:07  zender
;  the wonderful widgeting cloud data viewer
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin cld() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cld_event, event
@cld_cmn.com
;
; In case event is not handled by the below
command_sng=''
;
event_type=tag_names(event,/structure)
;print,event_type,event
;
; Make sure all widgets have a usr_value before doing this
widget_control,event.id,get_uvalue=usr_value
;
if (usr_value eq 'top_buttons') then begin
; intentionally crash the program for easier restarts
	foo=non_existant_array(-1)
endif; endif
;
if (usr_value eq 'Wfl_out') then begin
if !d.name eq 'X' then begin
	widget_control,/hourglass
	pth_nm_get,fl_out,path,name
	new_fl=pickfile( $
		/read, $
		path=path, $
		file=fl_out, $
		filter=['*ps'])
	if (new_fl ne '') then begin
		fl_out=new_fl
		pth_nm_get,fl_out,path,name
		widget_control,Wfl_out,set_value=name
	endif
	command_sng='pickfile() = '+fl_out
endif else begin
	command_sng='Changing output filenames is only allowed when device is X windows'
endelse
endif; endif Wfl_out
;
if (usr_value eq 'Wformat') then begin
args=''
fmt_idx=event.index
case fmt_lst(fmt_idx).english of
	'X': begin
		close_ps,fl_nm=fl_out
	end
	'eps': begin
		result=strpos(fl_out,'.ps')
		args=args+fmt_lst(fmt_idx).arg
		args=args+chroma_lst(chroma_idx).arg
		args=args+shape_lst(shape_idx).arg
		args=args+',fl_nm='''+fl_out+'''
		args='open_ps'+args
		result=execute(args)
		command_sng=args
	end
	'gif': begin
		image=tvrd()
		raster_fl='/cgd/home/zender/web/'+src_nm+'/'+fld_nm+'.gif'
;		raster_fl=fl_out
		write_gif,raster_fl,image
		fmt_idx=0
		widget_control,Wformat,set_list_select=fmt_idx
		command_sng='Wrote currently displayed image to '+raster_fl
	end
	'ps': begin
		args=args+chroma_lst(chroma_idx).arg
		args=args+shape_lst(shape_idx).arg
		args=args+',fl_nm='''+fl_out+'''
		args='open_ps'+args
		result=execute(args)
		command_sng=args
	end
	else:
endcase
endif; endif Wformat
;
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
;
if (usr_value eq 'Wshape') then begin
shape_idx=event.index
command_sng='shape = '+shape_lst(shape_idx).english
endif; endif Wshape
;
if (usr_value eq 'Wchroma') then begin
chroma_idx=event.index
command_sng='chroma = '+chroma_lst(chroma_idx).english
endif; endif Wchroma
;
if (usr_value eq 'fl_buttons') then begin
	case event.value of
		'Primary': begin
		filename=pickfile( $
				/read, $
				path='/data/zender/_aux0_/cld', $
				file=primary_fl, $
				filter=['*.nc'])
		if (filename ne '') then primary_fl=filename
		end
		'Secondary': begin
		filename=pickfile( $
				/read, $
				path='/data/zender/data', $
				file=secondary_fl, $
				filter=['*.nc'])
		if (filename ne '') then secondary_fl=filename
		end
		'Tertiary': begin
		filename=pickfile( $
				/read, $
				path='/data/zender/data', $
				file=tertiary_fl, $
				filter=['*.nc'])
		if (filename ne '') then tertiary_fl=filename
		end
		'Quaternary': begin
		filename=pickfile( $
				/read, $
				path='/data/zender/data', $
				file=quaternary_fl, $
				filter=['*.nc'])
		if (filename ne '') then quaternary_fl=filename
		end
		else: 
	endcase	
;
	command_sng=event.value+" = pickfile() = "+filename
endif; endif fl_buttons
;
if (usr_value eq 'palette_button') then begin
;	print,"event.value = ",event.value
;	print,"usr_value = ",usr_value
	pll_idx=event.value
	command_sng="pll_idx=event.value = "+string(event.value)
endif; endif palette_button
;
if (event_type eq 'WIDGET_BUTTON') then begin
	widget_control,event.id,get_value=value
	case value of
		'Done': begin
			widget_control, /destroy, event.top
			command_sng="widget_control, /destroy, event.top"
		end
		else: print,'Unavailable Function'
	endcase	
endif	;end if WIDGET_BUTTON

if (usr_value eq 'fgr_list') then begin
	if (fgr_lst(event.index).fnc_nm ne '') then $
		clr_rfr,clr_tbl,'',0
		call_procedure,fgr_lst(event.index).fnc_nm
	command_sng=fgr_lst(event.index).fnc_nm
endif
;
if command_sng ne '' then print,command_sng
;
end	;end event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End cld() event handler
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin main program cld() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cld, $
	primary_fl=Pprimary_fl, $
	secondary_fl=Psecondary_fl, $
	tertiary_fl=Ptertiary_fl, $
	quaternary_fl=Pquaternary_fl, $
	color_order=Pcolor_order, $
	pll_idx=Ppalette, $
	clr_tbl=Pclr_tbl, $
	pause=pause, $
	ps=ps
;
; The function of the parameters prefixed with "P" is to 
; remap passed parameters that are needed in common blocks into named variables. 
;
;Example usage:
;
;cld,Pprimary_fl='/cgd/data/zender/data/H.T.final.50.nc'
;
@cld_cmn.com

if n_elements(primary_fl) eq 0 then primary_fl='/data/zender/data/H.T.final.50.nc'
if n_elements(secondary_fl) eq 0 then secondary_fl='/data/zender/data/H.T.final.35.nc'
if n_elements(tertiary_fl) eq 0 then tertiary_fl='/data/zender/data/S.T.final.50.nc'
if n_elements(quaternary_fl) eq 0 then quaternary_fl='/data/zender/data/S.T.final.35.nc'
if n_elements(color_order) eq 0 then color_order=1
if n_elements(pll_idx) eq 0 then pll_idx=3
if n_elements(clr_tbl) eq 0 then clr_tbl=34
if n_elements(pause) eq 0 then pause='y';ready to roll
if n_elements(ps) eq 0 then ps='n'; print to .ps file or printer?
if pause eq 'y' then !p.multi=0 else begin
        !p.multi=[0,3,1,0,0];=[0,num_cols,num_rows,0,0]
        pause='n'
endelse
;
sys_time=systime(1) 
if !d.name ne 'X' then device,/close_file
set_plot,'X'
close,/all
dest_idx=1
fmt_idx=0
shape_idx=3
chroma_idx=0
device=16
fl_out=getenv('PWD')+'/idl.ps'
pth_nm_get,fl_out,out_path,out_nm

; Set up the data structure
; Note: the fnc_nm field MUST be <= 15 characters
num_fgr_lst=22
foo={fgr,index:0,english:'',fnc_nm:'',args:''}
fgr_lst=replicate({fgr},num_fgr_lst)

fgr_lst(0)={fgr,0,'Albedo Emissivity of Time','alb_emis_of_time',''}
fgr_lst(1)={fgr,1,'Albedo Emissivity','alb_emis',''}
fgr_lst(2)={fgr,2,'CCN Activated','CCN_activated',''}
fgr_lst(3)={fgr,3,'Fall Speed B&W','fall_speed_bw',''}
fgr_lst(4)={fgr,4,'Heating Distribution B&W','heat_dst_bw',''}
fgr_lst(5)={fgr,5,'Heating Rates','heating_rates',''}
fgr_lst(6)={fgr,6,'Humidity Level','humidity',''}
fgr_lst(7)={fgr,7,'IWC B&W','IWC_bw',''}
fgr_lst(8)={fgr,8,'IWC Color','IWC_color',''}
fgr_lst(9)={fgr,9,'IWP Reff Albedo Emissivity','IWP_reff_alb_em',''}
fgr_lst(10)={fgr,10,'Local IWC','local_IWC',''}
fgr_lst(11)={fgr,11,'Mass Distribution B&W','mass_dst_bw',''}
fgr_lst(12)={fgr,12,'Mass Distribution Color','mass_dst_color',''}
fgr_lst(13)={fgr,13,'Radiative Forcing','rad_forcing',''}
fgr_lst(14)={fgr,14,'Size Distribution B&W','sz_dst_bw',''}
fgr_lst(15)={fgr,15,'Size Distribution Color','sz_dst_color',''}
fgr_lst(16)={fgr,16,'Temperature','temperature',''}
fgr_lst(17)={fgr,17,'Heating Rate Sensitivity','heat_rate_sns',''}
fgr_lst(18)={fgr,18,'Multi Radiative Forcing','multi_rad_forcing',''}
fgr_lst(19)={fgr,19,'Multi Albedo Emissivity','multi_alb_emis',''}
fgr_lst(20)={fgr,20,'Mass Mixing Ratio','mass_mix_ratio',''}
fgr_lst(21)={fgr,21,'Dens Temp Pres','dens_temp_pres',''}

fgr_lst=fgr_lst(sort(fgr_lst.english))
fgr_lst.index=indgen(num_fgr_lst)

; dest_lst: structure to hold valid hardcopy destinations
dest_nbr=6
foo={dest_sct,index:0,english:'',arg:''}
dest_lst=replicate({dest_sct},dest_nbr)
dest_lst(0)={dest_sct,0,'cms',',/cms'}
dest_lst(1)={dest_sct,0,'lp',',/lp'}
dest_lst(2)={dest_sct,0,'ruby',',/ruby'}
dest_lst(3)={dest_sct,0,'tags','/tags'}
dest_lst(4)={dest_sct,0,'tcpr',',/tcpr'}
dest_lst(5)={dest_sct,0,'zeke',',/zeke'}
dest_lst=dest_lst(sort(dest_lst.english))
dest_lst.index=indgen(dest_nbr)
;
; fmt_lst: structure to hold valid output fmt_lst
;
fmt_nbr=4
foo={fmt_sct,index:0,english:'',arg:''}
fmt_lst=replicate({fmt_sct},fmt_nbr)
fmt_lst(0)={fmt_sct,0,'X',''}
fmt_lst(1)={fmt_sct,0,'eps',',/eps'}
fmt_lst(2)={fmt_sct,0,'gif',''}
fmt_lst(3)={fmt_sct,0,'ps',''}
fmt_lst=fmt_lst(sort(fmt_lst.english))
fmt_lst.index=indgen(fmt_nbr)
;
; shape_lst: structure to hold valid hardcopy size options
;
shape_nbr=4
foo={shape_sct,index:0,english:'',arg:''}
shape_lst=replicate({shape_sct},shape_nbr)
shape_lst(0)={shape_sct,0,'half',',/half'}
shape_lst(1)={shape_sct,0,'land',',/land'}
shape_lst(2)={shape_sct,0,'port',',/port'}
shape_lst(3)={shape_sct,0,'square',',/square'}
shape_lst=shape_lst(sort(shape_lst.english))
shape_lst.index=indgen(shape_nbr)
;
; chroma_lst: structure to hold valid color options
;
chroma_nbr=3
foo={chroma_sct,index:0,english:'',arg:''}
chroma_lst=replicate({chroma_sct},chroma_nbr)
chroma_lst(0)={chroma_sct,0,'B&W',''}
chroma_lst(1)={chroma_sct,0,'Color',',/color'}
chroma_lst(2)={chroma_sct,0,'Tcpr',',/tcpr,/color'}
chroma_lst=chroma_lst(sort(chroma_lst.english))
chroma_lst.index=indgen(chroma_nbr)
;
; Set up the base widgets
base=widget_base(title='IDL Cloud Graphs',/row)
left_clm=widget_base(base,/frame,/column)
Wright_clm=widget_base(base,/frame,/column)
;
top_buttons=cw_bgroup( $
		Wright_clm, $
		['Foo','Moo','Goo'], $
		/row, $
		/exclusive, $
		uvalue='top_buttons')
;
done=widget_button( $
		Wright_clm, $
		value='Done', $
		uvalue='done_button')
;
fl_buttons=cw_bgroup( $
		left_clm, $
		['Primary','Secondary','Tertiary','Quaternary'], $
		/row, $
		/return_name, $
		label_left='Files:', $
		uvalue='fl_buttons')
;
palette_button=cw_bselector( $
		left_clm, $
		['Cool','Hue-Sat-Val','Linear Gray','Hue-Lgt-Sat'], $
		label_left='Palette:', $
		set_value=3, $
		/return_index, $
		uvalue='palette_button')
;
;print,"primary_fl = ",primary_fl
;txt_base=widget_base(left_clm,/column,/frame)
;primary_fl=cw_field( $
;		txt_base, $
;		title='Primary File:', $
;		/string, $
;		xsize=50, $
;		ysize=1, $
;		uvalue='primary_fl', $ 
;		value=primary_fl)
;print,"primary_fl = ",primary_fl
;secondary_fl=cw_field( $
;		txt_base, $
;		title='Secondary File:', $
;		/string, $
;		xsize=50, $
;		ysize=1, $
;		uvalue='secondary_fl', $ 
;		value=secondary_fl)
;
fgr_list=widget_list( $
		Wright_clm, $
		value=fgr_lst.english, $
		ysize=15, $
		uvalue='fgr_list') 
;
Wformat=widget_list( $
		Wright_clm, $
		value=fmt_lst.english, $
		ysize=fmt_nbr, $
		uvalue='Wformat') 
;
Wshape_idx=widget_list( $
		Wright_clm, $
		value=shape_lst.english, $
		ysize=shape_nbr, $
		uvalue='Wshape') 
;
Wchroma=widget_list( $
		Wright_clm, $
		value=chroma_lst.english, $
		ysize=chroma_nbr, $
		uvalue='Wchroma') 
;
Wdest=widget_list( $
		Wright_clm, $
		value=dest_lst.english, $
		ysize=dest_nbr, $
		uvalue='Wdest') 
;
Wfl_out=widget_button( $
		Wright_clm, $
		value=fl_out, $
		uvalue='Wfl_out')
;
draw=widget_draw( $
		left_clm, $
		xsize=512, $
		ysize=512)
;
widget_control,base,/realize
widget_control,get_value=window,draw
widget_control,Wformat,set_list_select=fmt_idx
widget_control,Wshape,set_list_select=shape_idx
widget_control,Wchroma,set_list_select=chroma_idx
wset,window
;
clr_rfr,clr_tbl,'',0
;
; Put up a default graph
IWP_reff_alb_em
;
; Call the event manager
xmanager, $
	'cld', $
	base, $
	event_handler='cld_event'
;
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End main program cld() 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
