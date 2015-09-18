; $Id$

pro clr_tbl_dmp, $
	clr_tbl_nbr=clr_tbl_nbr, $ ; color table number
	ctr_fct=ctr_fct		; contrast factor

; Example usage:
; clr_tbl_dmp,clr_tbl_nbr=22,ctr_fct=0.85 ; pistachio
; clr_tbl_dmp,clr_tbl_nbr=3,ctr_fct=0.85 ; red 
; clr_tbl_dmp,clr_tbl_nbr=8,ctr_fct=0.85 ; green
; clr_tbl_dmp,clr_tbl_nbr=1,ctr_fct=0.85 ; blue
; clr_tbl_dmp,clr_tbl_nbr=1,ctr_fct=0.90 ; blue

if n_elements(clr_tbl_nbr) eq 0 then clr_tbl_nbr=22
if n_elements(ctr_fct) eq 0 then ctr_fct=0.

@ibp_clr.com
set_plot,'ps'
loadct,clr_tbl_nbr
r_curr=r_orig
g_curr=g_orig
b_curr=b_orig

r_curr=r_curr+(255-r_curr)*ctr_fct
g_curr=g_curr+(255-g_curr)*ctr_fct
b_curr=b_curr+(255-b_curr)*ctr_fct

r_curr(0)=255
g_curr(0)=255
b_curr(0)=255

r_curr(255)=200
g_curr(255)=200
b_curr(255)=200

dbg=1
if dbg then begin
fl_nm_out='/home/zender/idl/bg.tbl'
prn_unit=1 
openw,prn_unit,fl_nm_out
for idx=0,n_elements(r_curr)-1 do begin
printf,prn_unit,string(format='(3(I3,1x))',r_curr(idx),g_curr(idx),b_curr(idx))
endfor; end loop over idx
close,prn_unit
endif; endif dbg

set_plot,'X'
end; end clr_tbl_dmp()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Refresh Colors routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro clr_rfr,clr_tbl_nbr,tbl_fl,usr_dfn_clr_tbl

@ibp_clr.com

; Enable superior dithering technique for B&W monitors

if !d.table_size eq 2 then device,/floyd	

; Load the color table then replace the background/foreground so text is 
; always black and background is always white

print,'usr_dfn_clr_tbl = ',usr_dfn_clr_tbl

if usr_dfn_clr_tbl then color_init,tbl_fl else loadct,clr_tbl_nbr

r_curr=r_orig
g_curr=g_orig
b_curr=b_orig

idx=n_elements(r_curr)-1
if !d.name eq 'X' then begin
	clr_wht_idx=0
	color_1_value=0
	clr_blk_idx=idx
endif else begin
	clr_blk_idx=0
	color_1_value=255
	clr_wht_idx=idx
endelse

if idx gt 2 then begin
if color_order eq 0 then begin
	r_curr=rotate(r_curr,2)
	g_curr=rotate(g_curr,2)
	b_curr=rotate(b_curr,2)
endif; endif color_order eq 0
endif; endif idx gt 2

white_color_val=255
if strpos(tbl_fl,'bg.tbl') eq -1 then black_color_val=0 else black_color_val=230
r_curr(clr_wht_idx)=white_color_val
g_curr(clr_wht_idx)=white_color_val
b_curr(clr_wht_idx)=white_color_val
r_curr(clr_blk_idx)=black_color_val
g_curr(clr_blk_idx)=black_color_val
b_curr(clr_blk_idx)=black_color_val
r_curr(1)=color_1_value
g_curr(1)=color_1_value
b_curr(1)=color_1_value

tvlct,r_curr,g_curr,b_curr

end; end clr_rfr()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Refresh Colors routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Zecca Color Table Initialization Routine
; This routine loads in home-grown ASCII format color tables in
; the form of a 256-line ordered triplet of RGB values.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro color_init, $
	tbl_fl

@ibp_clr.com

close,1
openr,1,tbl_fl
bt=intarr(3,256)
readf,1,bt
r_usr=reform(bt(0,*))
g_usr=reform(bt(1,*))
b_usr=reform(bt(2,*))
close,1

tbl_places=congrid(indgen(256),!d.table_size)
r_orig=interpolate(r_usr,tbl_places)
g_orig=interpolate(g_usr,tbl_places)
b_orig=interpolate(b_usr,tbl_places)

;tvlct,r_curr,g_curr,b_curr

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Zecca Color Table Routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Colorbar routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Id$
;+
; NAME:		CLR_BAR_DRW
; PURPOSE:  To provide for horizontal and vertical color scales 
;           with associated legends.
; CALLING SEQUENCE:
;       tv,dist(200),150,150
;       loadct,4
;       clr_bar_drw
; INPUTS:
;     	None
; KEYWORD PARAMETERS:
; bar_psn        Vector describing lower left and upper right coordinates
;            of the color bar as x1,y1,x2,y2. Default is NORMAL
;            coordinates.
; bar_clr_nbr   Number of colors in bar. Default is !d.table_size.
; bar_idx     Indices of colors to use. Default is 0,1,2,...
; bar_lgn     Array containing string/int/long/float values to be
;            used to annotate colors on bar. Default is no annotation.
; bar_fmt     Format to be used for int/long/float values in bar_lgn, Default is '(F10.1)'.
; bar_fnt       Font for bar_lgn text.
; bar_txt_clr      Color for bar_lgn text.
; bar_chr_sz   Size for bar_lgn text.
; label	     a label to place above the colorbar(usually unit)
; OUTPUTS:
;	Color bar type legend
; RESTRICTIONS:
;	Must be run on an imaging device
; TESTED ON THESE PLATFORMS:
;	SUN4
; MODIFICATION HISTORY:
;	SAL 9/91 , Precision Visuals, Inc.
;-
pro clr_bar_drw, $
	bar_chr_sz=bar_chr_sz, $
	bar_dvc=bar_dvc, $
        bar_fnt=bar_fnt, $
	bar_idx=bar_idx, $
	bar_lbl_sz=bar_lbl_sz, $
        bar_lgn=bar_lgn, $
	bar_last=bar_last, $
	bar_psn=bar_psn, $
	bar_txt_clr=bar_txt_clr, $
	bar_unit=bar_unit, $
	bar_clr_nbr=bar_clr_nbr, $
	dbg_lvl=dbg_lvl, $ ; Debugging level
	dat_min=dat_min, $ ; Minimum value of data plotted
	dat_max=dat_max ; Maximum value of data plotted

if n_elements(bar_chr_sz) eq 0 then bar_chr_sz=!p.charsize
if n_elements(bar_dvc) eq 0 then bar_dvc=0
if n_elements(bar_fnt) eq 0 then bar_fnt=!p.font
if n_elements(bar_idx) eq 0 then bar_idx=indgen(bar_clr_nbr)
if n_elements(bar_lbl_sz) eq 0 then begin
	if !d.name ne 'X' then bar_lbl_sz=1.5 else bar_lbl_sz=1.0
endif; endif bar_lbl_sz
if n_elements(bar_lgn) eq 0 then bar_lgn=strarr(bar_clr_nbr+1)
if n_elements(bar_psn) eq 0 then bar_psn=[0.1,0.1,0.2,0.9]
if n_elements(bar_txt_clr) eq 0 then bar_txt_clr=!p.color
if n_elements(bar_unit) eq 0 then bar_unit=''
if n_elements(bar_clr_nbr) eq 0 then bar_clr_nbr=!d.table_size
if n_elements(dat_min) eq 0 then dat_min=0.0
if n_elements(dat_max) eq 0 then dat_max=0.0
if n_elements(dbg_lvl) eq 0 then dbg_lvl=0

if dbg_lvl gt 0 then begin
	print,"bar_idx = ",bar_idx
	print,"bar_lgn = ",bar_lgn
	print,"bar_psn = ",bar_psn
	print,"bar_txt_clr = ",bar_txt_clr
	print,"clr_blk_idx = ",clr_blk_idx
	print,"r_curr(bar_txt_clr) = ",r_curr(bar_txt_clr)
	print,"g_curr(bar_txt_clr) = ",g_curr(bar_txt_clr)
	print,"b_curr(bar_txt_clr) = ",b_curr(bar_txt_clr)
	print,"r_curr(0) = ",r_curr(0)
	print,"g_curr(0) = ",g_curr(0)
	print,"b_curr(0) = ",b_curr(0)
	print,"r_curr(1) = ",r_curr(1)
	print,"g_curr(1) = ",g_curr(1)
	print,"b_curr(1) = ",b_curr(1)
	print,"r_curr(2) = ",r_curr(2)
	print,"g_curr(2) = ",g_curr(2)
	print,"b_curr(2) = ",b_curr(2)
	print,"r_curr(3) = ",r_curr(3)
	print,"g_curr(3) = ",g_curr(3)
	print,"b_curr(3) = ",b_curr(3)
endif; endif dbg

if bar_chr_sz eq 0 then bar_chr_sz=1.0
nbr_lgns=n_elements(bar_lgn)

; Convert bar_psn values to normal coordinates if /DEVICE set.
if bar_dvc eq 1 then begin
	bar_psn=float(bar_psn)
	bar_psn(0)=bar_psn(0)/!d.x_vsize
	bar_psn(1)=bar_psn(1)/!d.y_vsize
	bar_psn(2)=bar_psn(2)/!d.x_vsize
	bar_psn(3)=bar_psn(3)/!d.y_vsize
endif
chr_sz_y_nrm=bar_chr_sz*(float(!d.y_ch_size)/!d.y_vsize)
chr_sz_x_nrm=bar_chr_sz*(float(!d.x_ch_size)/!d.x_vsize)

if bar_psn(2)-bar_psn(0) lt bar_psn(3)-bar_psn(1) then begin

; Put annotation to right of bar (vertically)
; Calculate x_txt_psn = x coord for start of text
; Calculate bar_ncr = height of each bar
;x_txt_psn=bar_psn(2) + (bar_psn(2)-bar_psn(0))*.125
x_txt_psn=bar_psn(2)+chr_sz_x_nrm*0.1
bar_ncr=(bar_psn(3)-bar_psn(1))/bar_clr_nbr

; Loop for each color
for clr_crr=0,bar_clr_nbr-1 do begin

; Draw rectangle for this color
y_btm=bar_ncr*clr_crr + bar_psn(1)
y_top=y_btm + bar_ncr
polyfill, $
	[bar_psn(0),bar_psn(2),bar_psn(2),bar_psn(0)], $
	[y_btm,y_btm,y_top,y_top], $
	/normal, $
	color=bar_idx(clr_crr)

; Plot bar_lgn for this entry if one exists
if bar_lgn(clr_crr) then begin
	y_txt_psn=y_btm-bar_ncr*0.25
	xyouts, $
		x_txt_psn, $
		y_txt_psn, $
		bar_lgn(clr_crr), $
		/normal, $
		alignment=0.0, $
		font=bar_fnt, $
;		color=bar_txt_clr, $
		charsize=bar_chr_sz

	plots, $
		[bar_psn(0),bar_psn(2)], $
		[y_top,y_top], $
		linestyle=0, $
		thick=2.0, $
		color=bar_txt_clr, $
		/normal
;
endif
endfor; end loop over interior rectangles

; Plot bar_lgn for edge label if one exists
clr_crr=bar_clr_nbr
y_btm=bar_ncr*clr_crr+bar_psn(1)
if bar_lgn(clr_crr) then $ 
	xyouts, x_txt_psn, y_btm-bar_ncr*0.25, $
	bar_lgn(clr_crr), /normal, alignment=0, $
	font=bar_fnt, color=bar_txt_clr, charsize=bar_chr_sz

endif else begin; end vertical clr_bar_drw

; Put annotation below bar (horizontally)
; Calculate y_psn=y coord for start of text
; Calculate bar_ncr=width of each bar

;y_txt_psn=bar_psn(1)-(bar_psn(3)-bar_psn(1))*0.5
y_txt_psn=bar_psn(1)-chr_sz_y_nrm*1.1
bar_ncr=(bar_psn(2)-bar_psn(0))/bar_clr_nbr

; Loop for each color
for clr_crr=0,bar_clr_nbr-1 do begin

; Draw rectangle for this color
x_lft=bar_ncr*clr_crr + bar_psn(0)
x_rgt=x_lft + bar_ncr
polyfill, $
	[x_lft,x_rgt,x_rgt,x_lft], $
	[bar_psn(1),bar_psn(1),bar_psn(3),bar_psn(3)], $
	/normal, $
	color=bar_idx(clr_crr)

; Plot bar_lgn for this entry if one exists
if bar_lgn(clr_crr) then begin

	x_txt_psn=x_lft

	if clr_crr ne 0 and clr_crr ne bar_clr_nbr-1 then begin
		xyouts,x_txt_psn,y_txt_psn,bar_lgn(clr_crr),/normal,alignment=0.5,font=bar_fnt,color=bar_txt_clr,charsize=bar_chr_sz
		plots,[x_txt_psn,x_txt_psn],[bar_psn(1),bar_psn(3)],linestyle=0,thick=2.0,/normal
	endif; endif clr_crr ne 0 and clr_crr ne bar_clr_nbr-1

	if clr_crr eq 0 or clr_crr eq bar_clr_nbr-1 then begin
		xyouts,x_txt_psn,y_txt_psn,bar_lgn(clr_crr),/normal,alignment=0.5,font=bar_fnt,color=bar_txt_clr,charsize=bar_chr_sz
		plots,[x_txt_psn,x_txt_psn],[bar_psn(1),bar_psn(3)],linestyle=0,thick=2.0,/normal
	endif; endif clr_crr eq 0 or clr_crr eq bar_clr_nbr-1

endif; endif drawing clr_bar_drw label

endfor; end loop over interior rectangles

; Plot bar_lgn for edge label if one exists
clr_crr=bar_clr_nbr
x_lft=bar_ncr*clr_crr + bar_psn(0)
x_txt_psn=x_lft
if bar_lgn(clr_crr) then begin
	xyouts, $
		x_txt_psn, $
		y_txt_psn, $
		bar_lgn(clr_crr), $
		/normal, $
		alignment=0.5, $
		font=bar_fnt, $
		color=bar_txt_clr, $
		charsize=bar_chr_sz
endif; endif drawing clr_bar_drw label

endelse; end horizontal clr_bar_drw

; Draw box around colorbar
; 20000111 Is this necessary when each box is drawn individually?
plots, $
	[bar_psn(0)-0.000, $
	bar_psn(2)+0.000, $
	bar_psn(2)+0.000, $
	bar_psn(0)-0.000, $
	bar_psn(0)-0.000], $
	[bar_psn(1), $
	bar_psn(1), $
	bar_psn(3), $
	bar_psn(3), $
	bar_psn(1)], $
	linestyle=0, $
	thick=2.0, $
;	color=bar_txt_clr, $
	/normal

; Label units of colorbar
if bar_unit ne '' then $
	xyouts, $
		0.5*(bar_psn(0)+bar_psn(2)), $
		bar_psn(3)+0.013, $
		'('+bar_unit+')', $
		size=bar_lbl_sz, $
		alignment=0.5, $
		/NORMAL

end; end clr_bar_drw()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Colorbar routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Make Colors Routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro clr_mk, $
		nbr_colors, $
		clr_ord, $
		pll_idx, $
		background

common colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr

; Define a set of RGB color triples for colors 2 through nbr_colors+1.

; If by chance running on a monochrome workstation then do nothing
if !d.table_size le 2 then return

; for colors,
; clr_ord=1 colors run from blue..red
; clr_ord=0 colors run from red..blue

; for grey scale,
; clr_ord=1 colors run from black to white, except color 2 is white
; clr_ord=0 colors run from white to black

; When the background keyword is set then the continuous spectrum  
; of colors is actually displaced to begin at color index 3 because 
; color index 2 is reserved for a background color:
; 	color 2 is .8 grey for colorful colormaps
; 	color 2 is white for grey scale colormaps
; The same total number of colors is still set regardless of background.

; pll_idx=0  White
; pll_idx=1  Linear Grey
; pll_idx=2  Linear Blue
; pll_idx=3  HLS
; pll_idx=4  HSV
; pll_idx=5  grey
; pll_idx=6  delta
; pll_idx=7  Linear Grey (kludge for thesis, ok to delete)

; nbr_colors=number of colors to use starting with index 2 in the table.
; (index 0 and 1 are assumed set to black and white somewhere in the
; calling program)

if (pll_idx eq 0) then begin

; White Algorithm
; Reserve color 2 for grey background color in color plots
; That leaves nbr_colorsm1 left to specify in this white algorithm.
	nbr_unset_colors=nbr_colors
	color_idx=2
	if ((!d.table_size gt 2) and background) then begin
		color_set=color_idx
		r_curr(color_set)=0.8*255
		g_curr(color_set)=0.8*255
		b_curr(color_set)=0.8*255
		nbr_unset_colors=nbr_unset_colors-1
		color_idx=color_idx+1
	endif

	for color=0,nbr_unset_colors-1 do begin
		red = 255
		green = 255
		blue = 255
		color_idx=color+2
		if (clr_ord eq 1) then color_set=nbr_unset_colors-color_idx+3 else color_set=color_idx
		r_curr(color_set)=red
		g_curr(color_set)=green
		b_curr(color_set)=blue
		color_idx=color_idx+1
	endfor
endif 
if (pll_idx eq 1) then begin

; Greyscale algorithm
; Make sure greyscale map sets color 2 from grey to white
; so that backgrounds are in white in this palette
	nbr_unset_colors=nbr_colors
	color_idx=2
	if ((!d.table_size gt 2) and background) then begin
		color_set=color_idx
		r_curr(color_set)=255
		g_curr(color_set)=255
		b_curr(color_set)=255
		nbr_unset_colors=nbr_unset_colors-1
		color_idx=color_idx+1
	endif

	for color=0,nbr_unset_colors-1 do begin

; The 0.6 exponent looks the best on zeke
		red = 255*(color/float(nbr_unset_colors-1))^(0.6)
;		red = 255*color/float(nbr_unset_colors-1)
		green = red
		blue = red
		if (clr_ord eq 0) then color_set= $
			nbr_unset_colors-color_idx+3 else $
			color_set=color_idx
		r_curr(color_set)=red
		g_curr(color_set)=green
		b_curr(color_set)=blue
		color_idx=color_idx+1
	endfor
endif 
if (pll_idx eq 2) then begin

; Blue Algorithm
; Reserve color 2 for grey background color in color plots
; That leaves nbr_colorsm1 left to specify in this purple algorithm
	nbr_unset_colors=nbr_colors
	color_idx=2
	if ((!d.table_size gt 2) and background) then begin
		color_set=color_idx
		r_curr(color_set)=0.8*255
		g_curr(color_set)=0.8*255
		b_curr(color_set)=0.8*255
		nbr_unset_colors=nbr_unset_colors-1
		color_idx=color_idx+1
	endif

	for color=0,nbr_unset_colors-1 do begin
		red = 255*(color/float(nbr_unset_colors-1))
		green = 255-red
		blue = 255
		color_idx=color+2
		if (clr_ord eq 1) then color_set= $
			nbr_unset_colors-color_idx+3 else $
			color_set=color_idx
		r_curr(color_set)=red
		g_curr(color_set)=green
		b_curr(color_set)=blue
		color_idx=color_idx+1
	endfor
endif 
if (pll_idx eq 3) then begin

; Reserve color 2 for grey background color in color plots
; That leaves nbr_colorsm1 left to specify with this HLS algorithm
	nbr_unset_colors=nbr_colors
	color_idx=2
	special_offset=0
	if ((!d.table_size gt 2) and background) then begin
		color_set=color_idx
		r_curr(color_set)=0.8*255
		g_curr(color_set)=0.8*255
		b_curr(color_set)=0.8*255
		nbr_unset_colors=nbr_unset_colors-1
		color_idx=color_idx+1
		special_offset=1
	endif

	hue_min=324
	hue_max=360+hue_min
	for color=0,nbr_unset_colors-1 do begin
		hue_fraction=color/float(nbr_unset_colors-1)
		hue=hue_min+(hue_max - hue_min)*hue_fraction
		hue=hue mod 360
		color_convert,hue,0.60,0.75,red,green,blue,/HLS_RGB
		if (clr_ord eq 1) then color_set=nbr_unset_colors-color_idx+special_offset+4 else color_set=color_idx
;		if (clr_ord eq 1) then color_set=nbr_unset_colors-color_idx+3 else color_set=color_idx
		r_curr(color_set)=red
		g_curr(color_set)=green
		b_curr(color_set)=blue
		color_idx=color_idx+1
	endfor
endif 
if (pll_idx eq 4) then begin

; HSV algorithm
; Reserve color 2 for grey background color in color plots
; That leaves nbr_colorsm1 left to specify with this HSV algorithm
	nbr_unset_colors=nbr_colors
	color_idx=2
	special_offset=0
	if ((!d.table_size gt 2) and background) then begin
		color_set=color_idx
		r_curr(color_set)=0.8*255
		g_curr(color_set)=0.8*255
		b_curr(color_set)=0.8*255
		nbr_unset_colors=nbr_unset_colors-1
		color_idx=color_idx+1
		special_offset=1
	endif

	hue_min=-20
	hue_max=250
	for color=0,nbr_unset_colors-1 do begin
		hue_fraction=color/float(nbr_unset_colors-1)
		hue_fraction=hue_fraction+0.05*sin(4*3.1415*hue_fraction)
		hue=hue_min+(hue_max - hue_min)*hue_fraction
		color_convert,hue,1.0,1.0,red,green,blue,/HSV_RGB
		if (clr_ord eq 1) then color_set=nbr_unset_colors-color_idx+special_offset+4 else color_set=color_idx
		r_curr(color_set)=red
		g_curr(color_set)=green
		b_curr(color_set)=blue
		color_idx=color_idx+1
	endfor
endif 
if (pll_idx eq 5) then begin

; Fixed greyscale algorithm.
; This colormap is intended to provide a uniform, consistent spectrum of 
; greys for use with printed output.
	nbr_unset_colors=nbr_colors
	color_idx=2
	if ((!d.table_size gt 2) and background) then begin
		color_set=color_idx
		r_curr(color_set)=255
		g_curr(color_set)=255
		b_curr(color_set)=255
		nbr_unset_colors=nbr_unset_colors-1
		color_idx=color_idx+1
	endif

	for color=0,nbr_unset_colors-1 do begin

; The 0.6 exponent looks the best on zeke
		red = 255*(color/float(nbr_unset_colors-1))^(0.6)
;		red = 255*color/float(nbr_unset_colors-1)
		green = red
		blue = red
		if (clr_ord eq 0) then color_set= $
			nbr_unset_colors-color_idx+3 else $
			color_set=color_idx
		r_curr(color_set)=red
		g_curr(color_set)=green
		b_curr(color_set)=blue
		color_idx=color_idx+1
	endfor
endif 
if (pll_idx eq 6) then begin

; Delta colors algorithm.
; This colormap was written by SGC or MZ and produces a nice spectrum split evenly down the middle.
	nbr_unset_colors=nbr_colors
	color_idx=2
	if ((!d.table_size gt 2) and background) then begin
		color_set=color_idx
		r_curr(color_set)=255
		g_curr(color_set)=255
		b_curr(color_set)=255
		nbr_unset_colors=nbr_unset_colors-1
		color_idx=color_idx+1
	endif

	delta_colors,nbr_unset_colors,red,green,blue

	for color=0,nbr_unset_colors-1 do begin
		if (clr_ord eq 0) then color_set=nbr_unset_colors-color_idx+3 else color_set=color_idx
		r_curr(color_set)=red(color)
		g_curr(color_set)=green(color)
		b_curr(color_set)=blue(color)
		color_idx=color_idx+1
	endfor
endif 
if (pll_idx eq 7) then begin

; Greyscale algorithm.
; Make sure greyscale map sets color 2 from grey to white
; so that backgrounds are in white in this palette.
	nbr_unset_colors=nbr_colors
	color_idx=2
	if ((!d.table_size gt 2) and background) then begin
		color_set=color_idx
		r_curr(color_set)=255
		g_curr(color_set)=255
		b_curr(color_set)=255
		nbr_unset_colors=nbr_unset_colors-1
		color_idx=color_idx+1
	endif

	for color=0,nbr_unset_colors-1 do begin

; The 0.6 exponent looks the best on zeke
		red = 255*(color/float(nbr_unset_colors-1))^(0.6)
;		red = 255*color/float(nbr_unset_colors-1)
		green = red
		blue = red
		color_set=nbr_unset_colors-color_idx+3
;		color_set=color_idx
;		if (clr_ord eq 0) then color_set=nbr_unset_colors-color_idx+3 else color_set=color_idx
		r_curr(color_set)=red
		g_curr(color_set)=green
		b_curr(color_set)=blue
		color_idx=color_idx+1
	endfor
endif 
if (pll_idx) gt 7 then print,"Palette type ",pll_idx," is not supported."

;for color=0,nbr_colors+2 do begin
;	print,"color = ",color," rgb_curr(color) = ", r_curr(color),g_curr(color),b_curr(color)
;endfor

tvlct,r_curr,g_curr,b_curr

end; end clr_mk()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Make Color Routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Delta Color Routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro delta_colors,nbr_clr,red,green,blue

nc=255

r=intarr(nc)
g=intarr(nc)
b=intarr(nc)

pi=!pi
twopi   = 2.*pi
piovr2  = pi/2.
ncr     = nc*1.
ncovr2p1 = 129.
ncovr2p2 = 130.
;
;***  Negative values
;
;---  red curve 1
      ar1 = 0.
      xr1 = 0.
      wr1 = pi / 4.
;
;---  red curve 2
      ar2 = .7
      xr2 = pi
      wr2 = pi / 8.
;
;---  green curve 1
      ag1 = 1.
      xg1 = pi
      wg1 = pi / 2.
;
;---  green curve 2
      ag2 = 0.
      xg2 = pi
      wg2 = pi / 8.
;
;---  blue curve 1
      ab1 = 1.
      xb1 = pi / 3.
      wb1 = pi / 3.
;
;---  blue curve 2
      ab2 = .7
      xb2 = pi
      wb2 = pi / 8.
;
      dx = pi/ncovr2p2
      x = -dx
      for i = 0,ncovr2p1 do begin
         x = x + dx
;
         r(i) = round( ar1*ncr*exp(-((x-xr1)/wr1)^2) ) $
              + round( ar2*ncr*exp(-((x-xr2)/wr2)^2) )

         g(i) = round( ag1*ncr*exp(-((x-xg1)/wg1)^2) ) $
              + round( ag2*ncr*exp(-((x-xg2)/wg2)^2) )

         b(i) = round( ab1*ncr*exp(-((x-xb1)/wb1)^2) ) $
              + round( ab2*ncr*exp(-((x-xb2)/wb2)^2) )
;
         if (r(i) gt nc) then r(i) = nc
         if (g(i) gt nc) then g(i) = nc
         if (b(i) gt nc) then b(i) = nc
;
      endfor
;
;***  Positive values
;
;---  red curve 1
      ar1 = 1.
      xr1 = pi/2.
      wr1 = 10. * pi
;
;---  red curve 2
      ar2 = 0.
      xr2 = 5. * pi
      wr2 = pi / 5.
;
;---  green curve 1
      ag1 = 1.
      xg1 = 0.
      wg1 = pi / 2.
;
;---  green curve 2
      ag2 = 1.
      xg2 = 5. * pi
      wg2 = pi / 5.
;
;---  blue curve 1
      ab1 = 1.
      xb1 = pi
      wb1 = pi / 3.
;
;---  blue curve 2
      ab2 = .8
      xb2 = 0.
      wb2 = pi / 8.
;
      dx = pi/(ncr-ncovr2p2+1.)
      x = -dx
      for i = ncovr2p2,nc-1 do begin
         x = x + dx
;
         r(i) = round( ar1*ncr*exp(-((x-xr1)/wr1)^2) ) $
              + round( ar2*ncr*exp(-((x-xr2)/wr2)^2) )

         g(i) = round( ag1*ncr*exp(-((x-xg1)/wg1)^2) ) $
              + round( ag2*ncr*exp(-((x-xg2)/wg2)^2) )

         b(i) = round( ab1*ncr*exp(-((x-xb1)/wb1)^2) ) $
              + round( ab2*ncr*exp(-((x-xb2)/wb2)^2) )
;
         if (r(i) gt nc) then r(i) = nc
         if (g(i) gt nc) then g(i) = nc
         if (b(i) gt nc) then b(i) = nc
;
      endfor
;
	idx=round(nc*indgen(nbr_clr)/(nbr_clr-1.))
	red=r(idx)
	green=g(idx)
	blue=b(idx)
;
end; end delta_colors()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Delta Color Routine
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

