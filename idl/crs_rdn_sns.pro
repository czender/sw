;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CVS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: crs_rdn_sns.pro,v $
; $Source: /home/zender/cvs/idl/crs_rdn_sns.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; Purpose: All the IDL figures for the shape and size sensitivity paper.
;
; $Log: not supported by cvs2svn $
; Revision 1.2  2000-01-15 02:07:50  zender
; *** empty log message ***
;
; Revision 1.9  2000/01/10 19:33:32  zender
; *** empty log message ***
;
; Revision 1.8  2000/01/01 01:55:47  zender
; *** empty log message ***
;
; Revision 1.7  1999/12/31 20:12:45  zender
; *** empty log message ***
;
; Revision 1.6  1999/12/31 02:09:36  zender
; *** empty log message ***
;
; Revision 1.5  1999/12/31 00:18:14  zender
; *** empty log message ***
;
; Revision 1.4  1999/12/30 19:34:24  zender
; *** empty log message ***
;
; Revision 1.3  1999/10/04 23:37:11  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/03 16:52:03  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:46  zender
; Imported sources
;
; Revision 1.8  1993/08/26  19:10:50  zender
; this is the version sent to the reviewers, fully double spaced,
; still has glitches on some of the color plots.
;
; Revision 1.7  1993/08/06  17:37:51  zender
; last version before california. text looks nice. minor problems
; with the color figures include lousy axis numbers and some
; color spills.
;
; Revision 1.6  1993/08/02  17:14:24  zender
; fixed the contour plotting error, i was dividing by crystal_length
; not the delta_length! changed from plotting the effective_radius
; (which was bogus), to finding and plotting effective_length
;
; Revision 1.5  1993/07/28  20:09:14  zender
; added text legends to line graphs
;
; Revision 1.4  1993/07/13  21:30:02  zender
; have mad many little adjustments
;
; Revision 1.3  1993/07/12  13:40:30  zender
; another all nighter, improved figure 6 slightly, looking
; at possible systematic error in contour plotting the
; log of the crystal distributions.
;
; Revision 1.2  1993/07/09  19:17:26  zender
; Have added many more figures. still no color plots.
;
; Revision 1.1  1993/07/05  23:31:07  zender
; Initial revision
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
pro cir_rad_sns, $
	filename=filename, $
	num_layer=num_layer, $
	clr_tbl=clr_tbl, $
	pause=pause, $
	ps=ps
;
;Example usage:
;
;idl_template
;
if n_elements(filename) eq 0 then filename = 'clouds.nc'
if n_elements(num_layer) eq 0 then num_layer = 18
if n_elements(clr_tbl) eq 0 then clr_tbl = 34
if n_elements(pause) eq 0 then pause = 'y';ready to roll
if n_elements(ps) eq 0 then ps = 'n'; print to .ps file or printer?
if pause eq 'y' then !p.multi=0 else begin
        !p.multi=[0,3,1,0,0];=[0,num_cols,num_rows,0,0]
        pause = 'n'
endelse
;
sys_time = systime(1) 
close,/all
erase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Define some information useful to all the plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Define the primary filenames used in the figures
big_case="/cgd/data/zender/data/H.T.big.50.nc"
fast_case="/cgd/data/zender/data/H.T.reg.fast.50.nc"
heymsfield_case="/cgd/data/zender/data/H.M.heyms.40.nc"
high_case="/cgd/data/zender/data/H.T.high.50.nc"
init_conds="/cgd/data/zender/data/H.T.reg.init0.50.nc"
init_rad_reg="/cgd/data/zender/data/H.T.reg.init.50.nc"
init_rad_truncate="/cgd/data/zender/data/H.T.reg.init.40.nc"
init_rad_sphere="/cgd/data/zender/data/S.T.reg.init.50.nc"
low_case="/cgd/data/zender/data/H.T.low.50.nc"
radke_case="/cgd/data/zender/data/H.M.radke.50.nc"
sphere_case="/cgd/data/zender/data/S.T.reg.50.nc"
still_case="/cgd/data/zender/data/H.T.still.50.nc"
;truncate_case="/cgd/data/zender/data/H.T.reg.40.nc"
;sphere_case="/cgd/data/zender/data/S.T.reg.50.nc"
;control_case="/cgd/data/zender/data/H.T.reg.50.nc"
;control_case="/cgd/data/zender/data/H.T.reg.LWLIOU.50.nc"
control_case="/cgd/data/zender/data/H.T.final.50.nc"
truncate_case="/cgd/data/zender/data/H.T.final.35.nc"
sphere_case="/cgd/data/zender/data/S.T.final.50.nc"
sphere_truncate_case="/cgd/data/zender/data/S.T.final.35.nc"
;
; Set up some common contouring information
dist_contour_level=[ $
;  .03125,0.0625,.125, $
;  .25,.5,1.0,2.0,4., $
;  8.,16.,32.0,64.,128., $
;  256.,512.0,1024.,2048.,4096., $
;  8192.0,16384.,32768.,65536.,131072.0, $
;  262144.]
  0.0,.5,1.0,2.0,4., $
  8.,16.,32.0,64.,128., $
  256.,512.0,1024.,2048.,4096., $
  8192.0,16384.,32768.,65536.,131072.0, $
  262144.,524288.,1048576.,2097152.0,4194304.]
;
dist_contour_lbls=[ $
  '0.','.5','1','2','4', $
  '8','16','32','64','128', $
  '256','512','1k','2k','4k', $
  '8k','16k','32k','64k','128k', $
  '256k','512k','1024k','2048k','4096k'] 
;
log10_dst_contour_level=[ $
  -2.0,-1.5,-1.0,-.5,0.0, $
  .5,1.0,1.5,2.0,2.5,3.0, $
  3.5,4.,4.5,5.,5.5, $
  6.]  
;
dist_which_lbl=[ $
  0.0,.03125,0.0625,.125, $
  .25,.5,1.0,2.0,4., $
  8.,16.,32.0,64.,128., $
  256.,512.0,1024.,2048.,4096., $
  8192.0,16384.,32768.,65536.,131072.0, $
  262144.]
;
dist_ntt=[ $
  '0.','.03125','.0625','.125', $
  '.25','.5','1','2','4', $
  '8','16','32','64','128', $
  '256','512','1024','2048','4096', $
  '8192','16384','32768','65536','131072', $
  '262144']
;
neg_dst_contour_level=[ $
  -128,-64,-32,-16,-8, $
  -4,-2,-1,-.5,-.25, $
  0.0, $
  .25,.5,1,2,4, $
  8,16,32,64,128, $
  256,512,1024,2048,4096, $
  8192]
;
neg_dst_contour_level=[ $
  -128,-64,-32,-16,-8, $
  -4,-2,-1,-.5,-.25, $
  0.0, $
  .25,.5,1,2,4, $
  8,16,32,64,128, $
  256,512,1024,2048,4096, $
  8192]
;
fig7b_dst_contour_level=[ $
  -128,-64,-32,-16,-1, $
  -.5,-.4,-.3,-.2,-.1, $
  0.0, $
  .1,0.2,.3,.4,.5, $
  1,2,5,10,15, $
  20,25,30,35,40, $
  8192]
;
fig7b_which_lbl=[ $
  1,1,1,1,1, $	
  1,1,1,1,1, $	
  1, $
  1,1,1,1,1, $	
  1,1,1,1,0, $	
  1,0,1,0,1, $	
  1]
;
fig7b_ntt=[ $
  '-128','-64','-32','-16','-1', $
  '-.5','-.4','-.3','-.2','-.1', $
  '0', $
  '.1','.2','.3','.4','.5', $
  '1','2','5','10','15', $
  '20','25','30','35','40', $
  '8192']
;
neg_dstln_sty=[ $
  1,1,1,1,1, $
  1,1,1,1,1, $
  0, $
  0,0,0,0,0, $
  0,0,0,0,0, $
  0,0,0,0,0, $
  0]
;
neg_dst_thick=[ $
  2,2,2,2,2, $
  2,2,2,2,2, $
  2, $
  1,1,1,1,1, $
  1,1,1,1,1, $
  1,1,1,1,1, $
  1]
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End common information block
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 4a.
; Plot the comparison of albedo and emissivity for the three main cases.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename_1=control_case
cdfid_1=ncdf_open(filename_1)
filename_2=sphere_case
cdfid_2=ncdf_open(filename_2)
filename_3=truncate_case
cdfid_3=ncdf_open(filename_3)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
ncdf_varget1,cdfid_2,"num_step",num_step_2
ncdf_varget1,cdfid_2,"rad_step",rad_step_2
ncdf_varget1,cdfid_3,"num_step",num_step_3
ncdf_varget1,cdfid_3,"rad_step",rad_step_3
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"albedo_of_time",albedo_1
ncdf_varget,cdfid_1,"emissivity_of_time",emissivity_1
ncdf_varget,cdfid_2,"albedo_of_time",albedo_2
ncdf_varget,cdfid_2,"emissivity_of_time",emissivity_2
ncdf_varget,cdfid_3,"albedo_of_time",albedo_3
ncdf_varget,cdfid_3,"emissivity_of_time",emissivity_3
ncdf_varget,cdfid_1,"time_array",time_array
;
; say bye-bye
ncdf_close,cdfid_1
ncdf_close,cdfid_2
ncdf_close,cdfid_3
; End of NetCDF commands
;
plot, $
	emissivity_1(0:num_step_1-1), $
	albedo_1(0:num_step_1-1), $
	tit='!5 Emissivity vs. Albedo', $
	xtit='Emissivity !7e!5', $
	ytit='Albedo !8A!5', $
;	font=0, $				;uses postscript fonts
	xrange=[.85,1.01], $
	yrange=[0.1,0.6], $
	xstyle=1, $
	ystyle=1, $
	thick=2.0, $
	linestyle=0
;
oplot,	$
	emissivity_2(0:num_step_2-1), $
	albedo_2(0:num_step_2-1), $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	emissivity_3(0:num_step_3-1), $
	albedo_3(0:num_step_3-1), $
	thick=2.0, $
	linestyle=2
;
ln_lgn_x1=.22
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Control Cloud',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5Spherical',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5Truncated',size=txt_lgn_sz,/NORMAL
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 4a commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 4b.
; Plot the time-by-time evolution of the albedo difference to the control cloud.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
;albedo_2_1_compare=100.*albedo_2(0:num_step_1-1)/albedo_1(0:num_step_1-1)
;albedo_3_1_compare=100.*albedo_3(0:num_step_1-1)/albedo_1(0:num_step_1-1)
albedo_2_1_compare=100.*(albedo_1(0:num_step_1-1)-albedo_2(0:num_step_1-1))/ $
	albedo_1(0:num_step_1-1)
albedo_3_1_compare=100.*(albedo_1(0:num_step_1-1)-albedo_3(0:num_step_1-1))/ $ 
	albedo_1(0:num_step_1-1)
;
plot, $
	time_array(0:num_step_1-1)/60.0, $
	albedo_2_1_compare(0:num_step_1-1), $
	tit='!5 Time Evolution of Albedo Disparity', $
	xtit='Elapsed Time (min.)', $
;	ytit='Albedo as % of Control Case', $
	ytit='Percentage Error (%)', $
;	yrange=[50.0,100.0], $
	yrange=[0.0,50.0], $
	xstyle=1, $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	time_array(0:num_step_1-1)/60.0, $
	albedo_3_1_compare(0:num_step_1-1), $
	thick=2.0, $
	linestyle=2
;
ln_lgn_x1=.52
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=0.9
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(1),'!5Spherical',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5Truncated',size=txt_lgn_sz,/NORMAL
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 4b commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 4c.
; Plot the time-by-time evolution of the emissvity difference to the control cloud.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
emis_2_1_compare=100.*emissivity_2(0:num_step_1-1)/emissivity_1(0:num_step_1-1)
emis_3_1_compare=100.*emissivity_3(0:num_step_1-1)/emissivity_1(0:num_step_1-1)
;
plot, $
	time_array(0:num_step_1-1)/60.0, $
	emis_2_1_compare(0:num_step_1-1), $
	tit='!5 Time Evolution of Emissivity Disparity', $
	xtit='Elapsed Time (min.)', $
	ytit='Emissivity as % of Control Case', $
;	yrange=[50.0,100.0], $
	xstyle=1, $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	time_array(0:num_step_1-1)/60.0, $
	emis_3_1_compare(0:num_step_1-1), $
	thick=2.0, $
	linestyle=2
;
ln_lgn_x1=.52
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=.5
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(1),'!5Spherical',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5Truncated',size=txt_lgn_sz,/NORMAL
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 4c commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 1acolor.
; Contour plot the initial size distribution vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_sz] (in C)
; is accessed as         foo(1..num_sz,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename=init_conds
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_sz",num_sz
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"effective_radius",effective_radius
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
distribution=concentration
for layer=1,num_layer do begin
	distribution(1:num_sz,layer)=concentration(1:num_sz,layer)/ $
		delta_length(1:num_sz) ; #/m^3/m
endfor
;
min_level=min(distribution)/1.0e6
max_level=max(distribution)/1.0e6
;
num_contour_level=n_elements(dist_contour_level)
;
; enable superior dithering technique for B&W monitors
if !d.table_size eq 2 then device,/floyd	
;
; load the color table then replace the background/foreground so text is 
; always black and background is always white
common colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr
loadct,clr_tbl
if !d.name eq 'X' then begin
	clr_wht_idx=0
	color_1_value=0
	clr_blk_idx=n_elements(r_curr)-1
endif else begin
	clr_blk_idx=0
	color_1_value=255
	clr_wht_idx=n_elements(r_curr)-1
endelse
r_curr=r_orig
g_curr=g_orig
b_curr=b_orig
r_curr(clr_wht_idx)=255
g_curr(clr_wht_idx)=255
b_curr(clr_wht_idx)=255
r_curr(clr_blk_idx)=0
g_curr(clr_blk_idx)=0
b_curr(clr_blk_idx)=0
r_curr(1)=color_1_value
g_curr(1)=color_1_value
b_curr(1)=color_1_value
if(!d.table_size gt 2) then begin
	r_curr(2)=0.8*256
	g_curr(2)=0.8*256
	b_curr(2)=0.8*256
endif
;
num_cbar_colors=min([!d.table_sz-1,num_contour_level])
clr_mk,num_cbar_colors,1,3
;
tvlct,r_curr,g_curr,b_curr
;
cbar_fmt='(A5)'
cbar_idx=indgen(num_cbar_colors)+2
cbar_fnt=!p.font
cbar_txt_clr=clr_blk_idx
cbar_chr_sz=1.3
cbar_unit=""
cbar_lbl_sz=1.5
;
cbar_lgn=dist_contour_lbls
;
;cbar_lgn=strarr(num_cbar_colors)
;find the color index associated with each contour level
;for level=1,num_contour_level-2 do begin
;	color_idx=level+2
;	cbar_lgn(level)=dist_contour_lbls(level)
;endfor
;
cbar_psn=[ $
	.88, $ ; x_min
	.10, $ ; y_min
	.92, $ ; x_max
	.90] ; 7_max
;
plt_rgn_nrm=[ $
	.1, $ ; x_min
	.1, $ ; y_min
	.87, $ ; x_max
	.90] ; y_max
;
plt_rgn_dvc=convert_coord( $
	[plt_rgn_nrm(0), $
	plt_rgn_nrm(2)], $
	[plt_rgn_nrm(1), $
	plt_rgn_nrm(3)], $
	/normal, $
	/to_device)
;
plt_rgn_dvc=[ $
	plt_rgn_dvc(0), $
	plt_rgn_dvc(1), $
	plt_rgn_dvc(3), $
	plt_rgn_dvc(4)]
;	
clr_bar_drw, $
	bar_psn=cbar_psn, $
	num_bar_colors=num_cbar_colors, $
	bar_idx=cbar_idx, $
	bar_lgn=cbar_lgn, $
	bar_first=dist_contour_level(0), $
	bar_last=dist_contour_level(num_contour_level-1), $
	bar_fmt=cbar_fmt, $
	bar_fnt=cbar_fnt, $
	bar_txt_clr=cbar_txt_clr, $
	bar_chr_sz=cbar_chr_sz, $
	bar_unit=cbar_unit, $
	bar_lbl_sz=cbar_lbl_sz
;
title_sz=1.6
;
; note that contouring distribution(*,1:num_layer) always works whereas
; contouring distribution(1:num_sz,1:num_layer) sometimes crashes.
;
contour, $
;	alog10(distribution(1:num_sz,1:num_layer)/1.0e6 > 1.0e-2), $
	distribution(*,1:num_layer)/1.0e6, $
;	alog10(crystal_length(*)*1.0e6), $
	crystal_length(*)*1.0e6, $
	altitude(1:num_layer)/1000.0, $	
	tit='!5Distribution (m!E-3!N!7l!5m!E-1!N) of Control Cloud', $
	xtit='Crystal length !8L!5 (!7l!5m)', $
;	xtit='Log!I10!N(Crystal length !8L!5 (!7l!5m))', $
	ytit='Altitude !8z!5 (km)', $
;	level=log10_dst_contour_level, $      ;user supplies levels 4 all contours
	level=dist_contour_level, $            ;user supplies levels 4 all contours
;	c_labels=dist_which_lbl, $		;which levels to label
;	c_annotation=dist_ntt, $	;actual text labels
;	font=0, $				;uses postscript fonts
;	/follow, $				;labels every other contour automat.
;	/min_curve_surf, $			;makes contours smoother
	xrange=[0.0,600.0], $
;	xrange=[.5,3.], $
;	yrange=[8.,16.], $
	xstyle=1, $
	ystyle=1, $
	/closed, $				;IDL v 3.1	
	/fill, $				;IDL v 3.1	
	c_color=indgen(num_contour_level)+3, $	;IDL v 3.1	
	noclip=0, $ 				;IDL v 3.1	
	clip=[3,8,600,18], $			;IDL v 3.1	
	charsize=title_sz, $
	position=plt_rgn_nrm, $
	ycharsize=1.1, $
	xcharsize=1.1, $
	/noerase
;
contour, $
;	alog10(distribution(*,1:num_layer)/1.0e6 > 1.0e-2), $
	distribution(*,1:num_layer)/1.0e6, $
;	alog10(crystal_length(*)*1.0e6), $
	crystal_length(*)*1.0e6, $
	altitude(1:num_layer)/1000.0, $	
	tit='!5Distribution (m!E-3!N!7l!5m!E-1!N) of Control Cloud', $
	xtit='Crystal length !8L!5 (!7l!5m)', $
;	xtit='Log!I10!N(Crystal length !8L!5 (!7l!5m))', $
	ytit='Altitude !8z!5 (km)', $
;	level=log10_dst_contour_level, $      ;user supplies levels 4 all contours
	level=dist_contour_level, $            ;user supplies levels 4 all contours
;	c_labels=dist_which_lbl, $		;which levels to label
;	c_annotation=dist_ntt, $	;actual text labels
;	font=0, $				;uses postscript fonts
;	/follow, $				;labels every other contour
	xrange=[0.0,600.0], $
;	xrange=[.5,3.], $
;	yrange=[8.,16.], $
	xstyle=1, $
	ystyle=1, $
	noclip=0, $ 				;IDL v 3.1	
	clip=[3,8,600,18], $			;IDL v 3.1	
	charsize=title_sz, $
	position=plt_rgn_nrm, $
	ycharsize=1.1, $
	xcharsize=1.1, $
	/noerase
;
effective_length=effective_radius
for layer=1,num_layer do begin
	effective_length(layer)=equiv_rad_to_bullet_length(effective_radius(layer))
endfor
cloud_eff_length=equiv_rad_to_bullet_length(cloud_eff_radius(num_step))
;
oplot,	$
	effective_length(1:num_layer)*1.0e6, $
;	alog10(effective_length(1:num_layer)*1.0e6), $
	altitude(1:num_layer)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
	[cloud_eff_length,cloud_eff_length]*1.0e6, $
;	alog10([cloud_eff_length,cloud_eff_length]*1.0e6), $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 1acolor commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 1a.
; Contour plot the initial size distribution vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_sz] (in C)
; is accessed as         foo(1..num_sz,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename=init_conds
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_sz",num_sz
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"effective_radius",effective_radius
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
distribution=concentration
for layer=1,num_layer do begin
	distribution(1:num_sz,layer)=concentration(1:num_sz,layer)/ $
		delta_length(1:num_sz) ; #/m^3/m
endfor
;
contour, $
;	alog10(distribution(1:num_sz,1:num_layer)/1.0e6 > 1.0e-2), $
	distribution(1:num_sz,1:num_layer)/1.0e6, $
;	alog10(crystal_length(1:num_sz)*1.0e6), $
	crystal_length(1:num_sz)*1.0e6, $
	altitude(1:num_layer)/1000.0, $	
	tit='!5 Initial Distribution (#-m!E-3!N-!7l!5m!E-1!N) of Control Cloud', $
	xtit='Crystal length !8L!5 (!7l!5m)', $
;	xtit='Log!I10!N(Crystal length !8L!5 (!7l!5m))', $
	ytit='Altitude !8z!5 (km)', $
;	level=log10_dst_contour_level, $      ;user supplies levels 4 all contours
	level=dist_contour_level, $            ;user supplies levels 4 all contours
;	c_labels=dist_which_lbl, $		;which levels to label
;	c_annotation=dist_ntt, $	;actual text labels
;	font=0, $				;uses postscript fonts
	/follow, $				;labels every other contour
	xrange=[0.0,600.0], $
;	xrange=[.5,3.], $
;	yrange=[8.,16.], $
	xstyle=1, $
	ystyle=1, $
	/noerase
;
effective_length=effective_radius
for layer=1,num_layer do begin
	effective_length(layer)=equiv_rad_to_bullet_length(effective_radius(layer))
endfor
cloud_eff_length=equiv_rad_to_bullet_length(cloud_eff_radius(num_step))
;
oplot,	$
	effective_length(1:num_layer)*1.0e6, $
;	alog10(effective_length(1:num_layer)*1.0e6), $
	altitude(1:num_layer)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
	[cloud_eff_length,cloud_eff_length]*1.0e6, $
;	alog10([cloud_eff_length,cloud_eff_length]*1.0e6), $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 1a commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 1b.
; Contour plot the final size distribution vs. altitude
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_sz] (in C)
; is accessed as         foo(1..num_sz,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename=control_case
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_sz",num_sz
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"effective_radius",effective_radius
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
distribution=concentration
for layer=1,num_layer do begin
	distribution(1:num_sz,layer)=concentration(1:num_sz,layer)/ $
		delta_length(1:num_sz) ; #/m^3/m
endfor
;
contour, $
;	alog10(distribution(1:num_sz,1:num_layer)/1.0e6 > 1.0e-2), $
	distribution(1:num_sz,1:num_layer)/1.0e6, $
;	alog10(crystal_length(1:num_sz)*1.0e6), $
	crystal_length(1:num_sz)*1.0e6, $
	altitude(1:num_layer)/1000.0, $	
	tit='!5 Final Distribution (#-m!E-3!N-!7l!5m!E-1!N) of Control Cloud', $
	xtit='Crystal length !8L!5 (!7l!5m)', $
;	xtit='Log!I10!N(Crystal length !8L!5 (!7l!5m))', $
	ytit='Altitude !8z!5 (km)', $
;	level=log10_dst_contour_level, $      ;user supplies levels 4 all contours
	level=dist_contour_level, $            ;user supplies levels 4 all contours
;	c_labels=dist_which_lbl, $		;which levels to label
;	c_annotation=dist_ntt, $	;actual text labels
;	font=0, $				;uses postscript fonts
	/follow, $				;labels every other contour
	xrange=[0.0,600.0], $
;	xrange=[.5,3.], $
;	yrange=[8.,16.], $
	xstyle=1, $
	ystyle=1, $
;	/min_curve_surf, $
	/noerase
;
effective_length=effective_radius
for layer=1,num_layer do begin
	effective_length(layer)=equiv_rad_to_bullet_length(effective_radius(layer))
endfor
cloud_eff_length=equiv_rad_to_bullet_length(cloud_eff_radius(num_step))
;
oplot,	$
	effective_length(1:num_layer)*1.0e6, $
;	alog10(effective_length(1:num_layer)*1.0e6), $
	altitude(1:num_layer)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
	[cloud_eff_length,cloud_eff_length]*1.0e6, $
;	alog10([cloud_eff_length,cloud_eff_length]*1.0e6), $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 1b commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 5
; Plot the time evolution of radiative forcing for the three main cases.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename_1=control_case
cdfid_1=ncdf_open(filename_1)
filename_2=sphere_case
cdfid_2=ncdf_open(filename_2)
filename_3=truncate_case
cdfid_3=ncdf_open(filename_3)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
ncdf_varget1,cdfid_2,"num_step",num_step_2
ncdf_varget1,cdfid_2,"rad_step",rad_step_2
ncdf_varget1,cdfid_3,"num_step",num_step_3
ncdf_varget1,cdfid_3,"rad_step",rad_step_3
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"LW_cloud_forcing",LW_forcing_1
ncdf_varget,cdfid_1,"SW_cloud_forcing",SW_forcing_1
ncdf_varget,cdfid_1,"time_array",time_array_1
ncdf_varget,cdfid_2,"LW_cloud_forcing",LW_forcing_2
ncdf_varget,cdfid_2,"SW_cloud_forcing",SW_forcing_2
ncdf_varget,cdfid_2,"time_array",time_array_2
ncdf_varget,cdfid_3,"LW_cloud_forcing",LW_forcing_3
ncdf_varget,cdfid_3,"SW_cloud_forcing",SW_forcing_3
ncdf_varget,cdfid_3,"time_array",time_array_3
;
; say bye-bye
ncdf_close,cdfid_1
ncdf_close,cdfid_2
ncdf_close,cdfid_3
; End of NetCDF commands
;
; Convert to Diurnally Averaged forcings...
; The 1/pi factor is the diurnal average, and the 2/sqrt(3) factor converts
; the 30 degree zenith angle simulations into pseudo local-noon quantities.
SW_forcing_1=SW_forcing_1*2./(sqrt(3.)*!pi)
SW_forcing_2=SW_forcing_2*2./(sqrt(3.)*!pi)
SW_forcing_3=SW_forcing_3*2./(sqrt(3.)*!pi)
;
; ...and compute net forcings
net_forcing_1=SW_forcing_1+LW_forcing_1
net_forcing_2=SW_forcing_2+LW_forcing_2
net_forcing_3=SW_forcing_3+LW_forcing_3
;
plot, $
	time_array_1(0:num_step_1-1)/60.0, $
	net_forcing_1(0:num_step_1-1), $
	tit='!5Evolution of Diurnally Averaged Cloud Forcing', $
	xtit='Elapsed Time (min.)', $
	ytit='Forcing (W-m!E-2!N)', $
	charsize=1.2, $
	yrange=[-200.0,200.0], $
	xstyle=1, $
	thick=2.0, $
;	xmargin=[0,0], $ 
;	ymargin=[0,0], $ 
	linestyle=0
;
oplot,	$
	time_array_2(0:num_step_2-1)/60.0, $
	net_forcing_2(0:num_step_2-1), $
	thick=2.0, $
	linestyle=0
;
oplot,	$
	time_array_3(0:num_step_3-1)/60.0, $
	net_forcing_3(0:num_step_3-1), $
	thick=2.0, $
	linestyle=0
;
; Now overplot the SW forcings
oplot,	$
	time_array_1(0:num_step_1-1)/60.0, $
	SW_forcing_1(0:num_step_1-1), $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	time_array_2(0:num_step_2-1)/60.0, $
	SW_forcing_2(0:num_step_2-1), $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	time_array_3(0:num_step_3-1)/60.0, $
	SW_forcing_3(0:num_step_3-1), $
	thick=2.0, $
	linestyle=1
;
; Now overplot the LW forcings
oplot,	$
	time_array_1(0:num_step_1-1)/60.0, $
	LW_forcing_1(0:num_step_1-1), $
	thick=2.0, $
	linestyle=2
;
oplot,	$
	time_array_2(0:num_step_2-1)/60.0, $
	LW_forcing_2(0:num_step_2-1), $
	thick=2.0, $
	linestyle=2
;
oplot,	$
	time_array_3(0:num_step_3-1)/60.0, $
	LW_forcing_3(0:num_step_3-1), $
	thick=2.0, $
	linestyle=2
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5SWCF',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5LWCF',size=txt_lgn_sz,/NORMAL
;
txt_lgn_x=10
lgn_y=[65,35,-05, $
	-90,-126,-165, $
	165,140]
;
xyouts,txt_lgn_x,lgn_y(0),'!5Truncated',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-6
xyouts,txt_lgn_x,lgn_y(1),'!5Spherical',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-6
xyouts,txt_lgn_x,lgn_y(2),'!5Control',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-6
;
xyouts,txt_lgn_x,lgn_y(3),'!5Truncated',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-8
xyouts,txt_lgn_x,lgn_y(4),'!5Spherical',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-8
xyouts,txt_lgn_x,lgn_y(5),'!5Control',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-8
;
xyouts,txt_lgn_x,lgn_y(6),'!5Spherical and Control',size=txt_lgn_sz,/DATA, $
	ORIENTATION=2
xyouts,txt_lgn_x,lgn_y(7),'!5Truncated',size=txt_lgn_sz,/DATA, $
	ORIENTATION=2
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 5 commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 7a.
; Plot the vertical heating rate of the high thin cloud
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename=high_case
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"heating_rate_LW",LW_heating
ncdf_varget,cdfid,"heating_rate_SW",SW_heating
ncdf_varget,cdfid,"heating_rate_net",net_heating
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
plot, $
	net_heating(1:num_layer)*3600.0, $
	altitude(1:num_layer)/1000.0, $
	tit='!5 Integrated Heating Rates in High Thin Cirrus', $
	xtit='Heating Rate (K-hr!E-1!N)', $
	ytit='Altitude (km)', $
	xrange=[-.5,1.0], $
	yrange=[14.,18.], $
	ystyle=1, $
	thick=2.0, $
	linestyle=0
;
oplot,	$
	SW_heating(1:num_layer)*3600.0, $
	altitude(1:num_layer)/1000.0, $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	LW_heating(1:num_layer)*3600.0, $
	altitude(1:num_layer)/1000.0, $
	thick=1.0, $
	linestyle=2
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	  lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5SW',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5LW',size=txt_lgn_sz,/NORMAL
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 7a commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 7b.
; Contour plot the final integrated heating rate distribution vs. altitude
; for the high thin cloud.
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_sz] (in C)
; is accessed as         foo(1..num_sz,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename=high_case
cdfid=ncdf_open(filename)
;
;now get the scalars 
ncdf_varget1,cdfid,"num_layer",num_layer
ncdf_varget1,cdfid,"num_sz",num_sz
ncdf_varget1,cdfid,"num_step",num_step
;
;now get the one-dimensional arrays
ncdf_varget,cdfid,"altitude",altitude
ncdf_varget,cdfid,"cloud_effective_radius",cloud_eff_radius
ncdf_varget,cdfid,"crystal_length",crystal_length
ncdf_varget,cdfid,"delta_length",delta_length
ncdf_varget,cdfid,"effective_radius",effective_radius
ncdf_varget,cdfid,"env_density",env_density
;
;now get the two-dimensional arrays
ncdf_varget,cdfid,"concentration",concentration
ncdf_varget,cdfid,"crystal_heat_SW",SW_int_heat
ncdf_varget,cdfid,"crystal_heat_LW",LW_int_heat
;
; say bye-bye
ncdf_close,cdfid
; End of NetCDF commands
;
spec_heat_const_pres=1005.	;J/kg/K
seconds_per_hour=3600.		
distribution=concentration
net_int_heat=concentration
for layer=1,num_layer do begin
	distribution(1:num_sz,layer)=concentration(1:num_sz,layer)/ $
		delta_length(1:num_sz) ; #/m^3/m
	net_int_heat(1:num_sz,layer)=distribution(1:num_sz,layer)* $
		(SW_int_heat(1:num_sz,layer)+LW_int_heat(1:num_sz,layer))
	net_int_heat(1:num_sz,layer)=net_int_heat(1:num_sz,layer)/ $
		(env_density(layer)*spec_heat_const_pres)	; K/s/m
endfor
net_int_heat=net_int_heat*seconds_per_hour/1.0e6	; K/hour/micron
net_int_heat=net_int_heat*1.0e3			; mK/hour/micron
print,max(net_int_heat)
print,min(net_int_heat)
;net_int_heat=net_int_heat*1.0e6/1.0e6	; micro-Watts/m^3/micron
;
contour,net_int_heat(1:num_sz,1:num_layer), $
	alog10(crystal_length(1:num_sz)*1.0e6), $
;	crystal_length(1:num_sz)*1.0e6, $
	altitude(1:num_layer)/1000.0, $	
	tit='!5 Heating Rate Distribution (!5mK-hr!E-1!N-!7l!5m!E-1!N) in High Thin Cirrus', $
;	xtit='Crystal length !8L!5 (!7l!5m)', $
	xtit='Log!I10!N(Crystal length !8L!5 (!7l!5m))', $
	ytit='Altitude !8z!5 (km)', $
	level=fig7b_dst_contour_level, $      ;user supplies levels 4 all contours
	c_labels=fig7b_which_lbl, $		;which levels to label
	c_linestyle=neg_dstln_sty, $	;user supplies styles 4 all contours
	c_thick=neg_dst_thick, $		;user supplies thicknesses
	c_annotation=fig7b_ntt, $	;actual text labels
;	font=0, $				;uses postscript fonts
	/follow, $				;labels every other contour
;	xrange=[0.0,1000.0], $
	xrange=[.5,3.], $
	yrange=[14.,18.], $
	xstyle=1, $
	ystyle=1, $
;	/min_curve_surf, $
	/noerase
;
effective_length=effective_radius
for layer=1,num_layer do begin
	effective_length(layer)=equiv_rad_to_bullet_length(effective_radius(layer))
	if altitude(layer)/1000. gt 17.2 then effective_length(layer)=1.0e-6
endfor
cloud_eff_length=equiv_rad_to_bullet_length(cloud_eff_radius(num_step))
;
oplot,	$
;	effective_length(1:num_layer)*1.0e6, $
	alog10(effective_length(1:num_layer)*1.0e6), $
	altitude(1:num_layer)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
;	[cloud_eff_length,cloud_eff_length]*1.0e6, $
	alog10([cloud_eff_length,cloud_eff_length]*1.0e6), $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 7b commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 6a.
; Plot the vertical heating rate profile of the control cloud
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename_1=control_case
cdfid_1=ncdf_open(filename_1)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_layer",num_layer_1
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"altitude",altitude_1
ncdf_varget,cdfid_1,"heating_rate_LW",LW_heating_1
ncdf_varget,cdfid_1,"heating_rate_SW",SW_heating_1
ncdf_varget,cdfid_1,"heating_rate_net",net_heating_1
;
; say bye-bye
ncdf_close,cdfid_1
; End of NetCDF commands
;
plot, $
	net_heating_1(1:num_layer_1)*3600.0, $
	altitude_1(1:num_layer_1)/1000.0, $
	tit='!5 Heating Rates in Control Cloud', $
	xtit='Heating Rate (K-hr!E-1!N)', $
	ytit='Altitude (km)', $
	xrange=[-1.0,1.1], $
;	yrange=[8.,18.], $
	ystyle=1, $
	thick=2.0, $
	linestyle=0
;
oplot,	$
	SW_heating_1(1:num_layer_1)*3600.0, $
	altitude_1(1:num_layer_1)/1000.0, $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	LW_heating_1(1:num_layer_1)*3600.0, $
	altitude_1(1:num_layer_1)/1000.0, $
	thick=1.0, $
	linestyle=2
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5SW',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5LW',size=txt_lgn_sz,/NORMAL
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 6a commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 6b.
; Plot the (hex - sphere)/abs(hex) heating rate profile (the shape effect)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename_2=sphere_case
cdfid_2=ncdf_open(filename_2)
;
;now get the scalars 
ncdf_varget1,cdfid_2,"num_layer",num_layer_2
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_2,"altitude",altitude_2
ncdf_varget,cdfid_2,"heating_rate_LW",LW_heating_2
ncdf_varget,cdfid_2,"heating_rate_SW",SW_heating_2
ncdf_varget,cdfid_2,"heating_rate_net",net_heating_2
;
; say bye-bye
ncdf_close,cdfid_2
; End of NetCDF commands
;
; Compute the difference due to the shape effect
net_2_1_compare=altitude_2
LW_2_1_compare=altitude_2
SW_2_1_compare=altitude_2
epsilon=0.1/3600.
;
net_difference=net_heating_2-net_heating_1
net_2_1_compare=100.*net_difference*abs(net_heating_1)/ $
	(net_heating_1*net_heating_1+epsilon*epsilon)
LW_difference=LW_heating_2-LW_heating_1
LW_2_1_compare=100.*LW_difference*abs(LW_heating_1)/ $
	(LW_heating_1*LW_heating_1+epsilon*epsilon)
SW_difference=SW_heating_2-SW_heating_1
SW_2_1_compare=100.*SW_difference*abs(SW_heating_1)/ $
	(SW_heating_1*SW_heating_1+epsilon*epsilon)
;
;net_2_1_compare=net_difference*3600.
;LW_2_1_compare=LW_difference*3600.
;SW_2_1_compare=SW_difference*3600.
;
plot, $
	net_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	tit='!5 Shape Effect on Heating Rates', $
	xtit='Percentage Error = 100*(Spheres - Columns)/!10"!5Columns!10"!5', $
;	xtit='Spheres - Columns (K hr!E-1!N)', $
	ytit='Altitude (km)', $
;	xrange=[-.3,.3], $
	xrange=[-100.0,100.0], $
;	yrange=[0.1,0.6], $
	xstyle=1, $
	ystyle=1, $
	thick=2.0, $
	linestyle=0
;
oplot,	$
	SW_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	LW_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	thick=1.0, $
	linestyle=2
;
ln_lgn_x1=.22
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5SW',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5LW',size=txt_lgn_sz,/NORMAL
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 6b commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 6c.
; Plot the (control - truncated)/abs(control) heating rate profile 
; (the size effect)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename_3=truncate_case
cdfid_3=ncdf_open(filename_3)
;
;now get the scalars 
ncdf_varget1,cdfid_3,"num_layer",num_layer_3
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_3,"altitude",altitude_3
ncdf_varget,cdfid_3,"heating_rate_LW",LW_heating_3
ncdf_varget,cdfid_3,"heating_rate_SW",SW_heating_3
ncdf_varget,cdfid_3,"heating_rate_net",net_heating_3
;
; say bye-bye
ncdf_close,cdfid_3
; End of NetCDF commands
;
; Compute the difference due to the shape effect
net_3_1_compare=altitude_3
LW_3_1_compare=altitude_3
SW_3_1_compare=altitude_3
epsilon=0.1/3600.
;
net_difference=net_heating_3-net_heating_1
net_3_1_compare=100.*net_difference*abs(net_heating_1)/ $
	(net_heating_1*net_heating_1+epsilon*epsilon)
LW_difference=LW_heating_3-LW_heating_1
LW_3_1_compare=100.*LW_difference*abs(LW_heating_1)/ $
	(LW_heating_1*LW_heating_1+epsilon*epsilon)
SW_difference=SW_heating_3-SW_heating_1
SW_3_1_compare=100.*SW_difference*abs(SW_heating_1)/ $
	(SW_heating_1*SW_heating_1+epsilon*epsilon)
;
;net_3_1_compare=net_difference*3600.
;LW_3_1_compare=LW_difference*3600.
;SW_3_1_compare=SW_difference*3600.
;
plot, $
	net_3_1_compare(1:num_layer_3), $
	altitude_3(1:num_layer_3)/1000.0, $
	tit='!5 Truncation Effect on Heating Rates', $
	xtit='Percentage Error = 100*(Truncated - Control)/!10"!5Control!10"!5', $
;	xtit='Truncated - Control (K hr!E-1!N)', $
	ytit='Altitude (km)', $
;	xrange=[-.3,.3], $
	xrange=[-100.0,100.0], $
;	yrange=[0.1,0.6], $
	xstyle=1, $
	ystyle=1, $
	thick=2.0, $
	linestyle=0
;
oplot,	$
	SW_3_1_compare(1:num_layer_3), $
	altitude_3(1:num_layer_3)/1000.0, $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	LW_3_1_compare(1:num_layer_3), $
	altitude_3(1:num_layer_3)/1000.0, $
	thick=1.0, $
	linestyle=2
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5SW',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5LW',size=txt_lgn_sz,/NORMAL
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 6c commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 8.
; Plot the comparison of albedo and emissivity for the varying initial
; conditions and distributions.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename_1=still_case
cdfid_1=ncdf_open(filename_1)
filename_2=fast_case
cdfid_2=ncdf_open(filename_2)
filename_3=low_case
cdfid_3=ncdf_open(filename_3)
filename_4=heymsfield_case
cdfid_4=ncdf_open(filename_4)
filename_5=radke_case
cdfid_5=ncdf_open(filename_5)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_1,"rad_step",rad_step_1
ncdf_varget1,cdfid_2,"num_step",num_step_2
ncdf_varget1,cdfid_2,"rad_step",rad_step_2
ncdf_varget1,cdfid_3,"num_step",num_step_3
ncdf_varget1,cdfid_3,"rad_step",rad_step_3
ncdf_varget1,cdfid_4,"num_step",num_step_4
ncdf_varget1,cdfid_4,"rad_step",rad_step_4
ncdf_varget1,cdfid_5,"num_step",num_step_5
ncdf_varget1,cdfid_5,"rad_step",rad_step_5
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"albedo_of_time",albedo_1
ncdf_varget,cdfid_1,"emissivity_of_time",emissivity_1
ncdf_varget,cdfid_2,"albedo_of_time",albedo_2
ncdf_varget,cdfid_2,"emissivity_of_time",emissivity_2
ncdf_varget,cdfid_3,"albedo_of_time",albedo_3
ncdf_varget,cdfid_3,"emissivity_of_time",emissivity_3
ncdf_varget,cdfid_4,"albedo_of_time",albedo_4
ncdf_varget,cdfid_4,"emissivity_of_time",emissivity_4
ncdf_varget,cdfid_5,"albedo_of_time",albedo_5
ncdf_varget,cdfid_5,"emissivity_of_time",emissivity_5

ncdf_varget,cdfid_1,"time_array",time_array
;
; say bye-bye
ncdf_close,cdfid_1
ncdf_close,cdfid_2
ncdf_close,cdfid_3
ncdf_close,cdfid_4
ncdf_close,cdfid_5
; End of NetCDF commands
;
plot, $
	emissivity_1(0:num_step_1-1), $
	albedo_1(0:num_step_1-1), $
	tit='!5 Emissivity vs. Albedo', $
	xtit='Emissivity !7e!5', $
	ytit='Albedo !8A!5', $
;	font=0, $				;uses postscript fonts
	xrange=[0.1,1.05], $
	yrange=[0.0,0.7], $
	xstyle=1, $
	ystyle=1, $
	thick=2.0, $
	linestyle=0
;
oplot,	$
	emissivity_2(0:num_step_2-1), $
	albedo_2(0:num_step_2-1), $
	thick=2.0, $
	linestyle=1
;
oplot,	$
	emissivity_3(0:num_step_3-1), $
	albedo_3(0:num_step_3-1), $
	thick=2.0, $
	linestyle=2
;
oplot,	$
	emissivity_4(0:num_step_4-1), $
	albedo_4(0:num_step_4-1), $
	thick=2.0, $
	linestyle=3
;
oplot,	$
	emissivity_5(0:num_step_5-1), $
	albedo_5(0:num_step_5-1), $
	thick=2.0, $
	linestyle=4
;
ln_lgn_x1=.22
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.2
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	  lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=2.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(4)+0.013,linestyle=4,thick=2.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5No updraft',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5Fast updraft',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5Thin-low',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(3),'!5Heymsfield-Platt',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(4),'!5Dowling-Radke',size=txt_lgn_sz,/NORMAL
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 8 commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure 3.
; Contour plot the difference in final size distributions between 
; the hexagonal columns and the spheres.
; Recall that two-dimensional arrays are accessed in transposed 
; Fortran-like order, so foo[1..num_layer][1..num_sz] (in C)
; is accessed as         foo(1..num_sz,1..num_layer)  (in IDL)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
filename_1=control_case
cdfid_1=ncdf_open(filename_1)
filename_2=sphere_case
cdfid_2=ncdf_open(filename_2)
;
;now get the scalars 
ncdf_varget1,cdfid_1,"num_layer",num_layer_1
ncdf_varget1,cdfid_1,"num_sz",num_sz_1
ncdf_varget1,cdfid_1,"num_step",num_step_1
ncdf_varget1,cdfid_2,"num_layer",num_layer_2
ncdf_varget1,cdfid_2,"num_sz",num_sz_2
ncdf_varget1,cdfid_2,"num_step",num_step_2
;
;now get the one-dimensional arrays
ncdf_varget,cdfid_1,"altitude",altitude_1
ncdf_varget,cdfid_1,"cloud_effective_radius",cloud_eff_radius_1
ncdf_varget,cdfid_1,"crystal_length",crystal_length_1
ncdf_varget,cdfid_1,"delta_length",delta_length_1
ncdf_varget,cdfid_1,"effective_radius",effective_radius_1
ncdf_varget,cdfid_2,"altitude",altitude_2
ncdf_varget,cdfid_2,"cloud_effective_radius",cloud_eff_radius_2
ncdf_varget,cdfid_2,"crystal_length",crystal_length_2
ncdf_varget,cdfid_2,"delta_length",delta_length_2
ncdf_varget,cdfid_2,"effective_radius",effective_radius_2
;
;now get the two-dimensional arrays
ncdf_varget,cdfid_1,"concentration",concentration_1
ncdf_varget,cdfid_2,"concentration",concentration_2
;
; say bye-bye
ncdf_close,cdfid_1
ncdf_close,cdfid_2
; End of NetCDF commands
;
distribution_1=concentration_1
distribution_2=concentration_2
for layer=1,num_layer_1 do begin
	distribution_1(1:num_sz_1,layer)=concentration_1(1:num_sz_1,layer)/ $
		delta_length_1(1:num_sz_1) ; #/m^3/m
endfor
for layer=1,num_layer_2 do begin
	distribution_2(1:num_sz_2,layer)=concentration_2(1:num_sz_2,layer)/ $
		delta_length_2(1:num_sz_2) ; #/m^3/m
endfor
;
diff_dst=distribution_1-distribution_2
contour,diff_dst(1:num_sz_1,1:num_layer_1)/1.0e6, $
	alog10(crystal_length_1(1:num_sz_1)*1.0e6), $
;	crystal_length_1(1:num_sz_1)*1.0e6, $
	altitude_1(1:num_layer_1)/1000.0, $	
	tit='!5 Difference in Final Distributions (#-m!E-3!N-!7l!5m!E-1!N) Columns - Spheres # m!E-3!N !7l!5m!E-1!N', $
;	xtit='Crystal length !8L!5 (!7l!5m)', $
	xtit='Log!I10!N(Crystal length !8L!5 (!7l!5m))', $
	ytit='Altitude !8z!5 (km)', $
	level=neg_dst_contour_level, $        ;user supplies levels 4 all contours
	c_linestyle=neg_dstln_sty, $	;user supplies styles 4 all contours
	c_thick=neg_dst_thick, $		;user supplies thicknesses
;	c_labels=neg_dst_which_lbl, $	;which levels to label
;	c_annotation=dist_ntt, $	;actual text labels
;	font=0, $				;uses postscript fonts
	/follow, $				;labels every other contour
	xrange=[.5,3.], $
;	xrange=[0.0,750.0], $
;	yrange=[8.,16.], $
	xstyle=1, $
	ystyle=1, $
	/noerase
;
effective_length_1=effective_radius_1
for layer=1,num_layer_1 do begin
	effective_length_1(layer)=equiv_rad_to_bullet_length(effective_radius_1(layer))
endfor
cloud_eff_length_1=equiv_rad_to_bullet_length(cloud_eff_radius_1(num_step_1))
;
oplot,	$
;	effective_length_1(1:num_layer_1)*1.0e6, $
	alog10(effective_length_1(1:num_layer_1)*1.0e6), $
	altitude_1(1:num_layer_1)/1000.0, $
	linestyle=2, $
	thick=2.
;
oplot,	$
;	[cloud_eff_length_1,cloud_eff_length_1]*1.0e6, $
	alog10([cloud_eff_length_1,cloud_eff_length_1]*1.0e6), $
	[0.0,100.0], $ 
	linestyle=0, $
	thick=2.
;
if pause eq 'y' then begin
	print,' Hit any key to continue, p to print this graph, or q to quit ...'
	junk = get_kbrd(1)
	if junk eq 'q' then goto, exit_gracefully
	if junk eq 'p' then begin
	        halfplot
;******************** Printing routine, should be identical to above ***************
;******************** End Printing routine ******************************
	        if ps eq 'n' then lpr else lps
	endif        
endif ;endif pause
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure 3 commands
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of procedure
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;***************enforced exit********************************
goto,exit_gracefully
;
exit_gracefully: foo=1
;
end




