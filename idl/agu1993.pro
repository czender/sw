;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CVS Identification
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; $Author: zender $
; $Date$
; $Id$
; $Revision$
; $Locker:  $
; $RCSfile: agu1993.pro,v $
; $Source: /home/zender/cvs/idl/agu1993.pro,v $
; $Id$
; $State: Exp $
;
; NB: get RCS formatting in IDL files by using rcs -U -c"; " foo.pro
;
; Purpose: All the IDL figures for the shape and size sensitivity paper.
;
; $Log: not supported by cvs2svn $
; Revision 1.8  2000-01-10 23:36:26  zender
; *** empty log message ***
;
; Revision 1.6  2000/01/01 01:55:47  zender
; *** empty log message ***
;
; Revision 1.5  1999/12/31 02:09:35  zender
; *** empty log message ***
;
; Revision 1.4  1999/12/31 00:18:14  zender
; *** empty log message ***
;
; Revision 1.3  1999/10/04 23:37:11  zender
; *** empty log message ***
;
; Revision 1.2  1999/10/03 16:52:03  zender
; *** empty log message ***
;
; Revision 1.1.1.1  1999/05/11 22:20:52  zender
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
pro agu, $
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
;init_conds="/cgd/data/zender/data/H.T.reg.init0.50.nc"
;init_rad_reg="/cgd/data/zender/data/H.T.reg.init.50.nc"
;init_rad_truncate="/cgd/data/zender/data/H.T.reg.init.40.nc"
;init_rad_sphere="/cgd/data/zender/data/S.T.reg.init.50.nc"
init_conds="/cgd/data/zender/data/H.T.final.init0.50.nc"
init_rad_reg="/cgd/data/zender/data/H.T.final.init1.50.nc"
init_rad_truncate="/cgd/data/zender/data/H.T.final.init1.35.nc"
init_rad_sphere="/cgd/data/zender/data/S.T.final.init1.50.nc"
low_case="/cgd/data/zender/data/H.T.low.50.nc"
radke_case="/cgd/data/zender/data/H.M.radke.50.nc"
still_case="/cgd/data/zender/data/H.T.still.50.nc"
;truncate_case="/cgd/data/zender/data/H.T.reg.40.nc"
;sphere_case="/cgd/data/zender/data/S.T.reg.50.nc"
;control_case="/cgd/data/zender/data/H.T.reg.50.nc"
;control_case="/cgd/data/zender/data/H.T.reg.LWLIOU.50.nc"
control_case="/cgd/data/zender/data/H.T.final.50.nc"
truncate_case="/cgd/data/zender/data/H.T.final.35.nc"
sphere_case="/cgd/data/zender/data/S.T.final.50.nc"
sphere_truncate_case="/cgd/data/zender/data/S.T.final.35.nc"
;control_case="/cgd/data/zender/data/H.T.final_4hr.50.nc"
;truncate_case="/cgd/data/zender/data/H.T.final_4hr.35.nc"
;
; Set up some common contouring information
dist_contour_level=[ $
;  .03125,0.0625,.125, $
;  .25,.5,1.0,2.,4., $
;  8.,16.,32.,64.,128., $
;  256.,512.,1024.,2048.,4096., $
;  8192.,16384.,32768.,65536.,131072., $
;  262144.]
  .5,1.0,2.0,4., $
  8.,16.,32.0,64.,128., $
  256.,512.0,1024.,2048.,4096., $
  8192.0,16384.,32768.,65536.,131072.0, $
  262144.,524288.,1048576.,2097152.0,4194304.]
;
log10_dst_contour_level=[ $
  -2.0,-1.5,-1.0,-.5,0.0, $
  .5,1.0,1.5,2.0,2.5,3., $
  3.5,4.,4.5,5.,5.5, $
  6.]  
;
dist_which_lbl=[ $
  .03125,0.0625,.125, $
  .25,.5,1.0,2.0,4., $
  8.,16.,32.0,64.,128., $
  256.,512.0,1024.,2048.,4096., $
  8192.0,16384.,32768.,65536.,131072.0, $
  262144.]
;
dist_ntt=[ $
  '.03125','.0625','.125', $
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
tvlct,r_curr,g_curr,b_curr
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
filename_2=truncate_case
cdfid_2=ncdf_open(filename_2)
filename_3=sphere_case
cdfid_3=ncdf_open(filename_3)
filename_4=sphere_truncate_case
cdfid_4=ncdf_open(filename_4)
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
ncdf_varget,cdfid_1,"time_array",time_array
;
; say bye-bye
ncdf_close,cdfid_1
ncdf_close,cdfid_2
ncdf_close,cdfid_3
ncdf_close,cdfid_4
; End of NetCDF commands
;
plot, $
	emissivity_1(0:num_step_1-1), $
	albedo_1(0:num_step_1-1), $
	tit='!5Albedo vs. Emissivity', $
	xtit='Emissivity !7e!5', $
	ytit='Albedo !8A!5', $
	xrange=[.94,1.0], $
	yrange=[.20,.4], $
	xstyle=0, $
	ystyle=0, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0

oplot,	$
	emissivity_2(0:num_step_2-1), $
	albedo_2(0:num_step_2-1), $
	thick=3.0, $
	linestyle=1
;
oplot,	$
	emissivity_3(0:num_step_3-1), $
	albedo_3(0:num_step_3-1), $
	thick=3.0, $
	linestyle=2
;
oplot,	$
	emissivity_4(0:num_step_4-1), $
	albedo_4(0:num_step_4-1), $
	thick=3.0, $
	linestyle=3
;
ln_lgn_x1=.3
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(3)+0.013,linestyle=3,thick=3.0,/NORMAL
;
xyouts,txt_lgn_x,lgn_y(0),'!5Control Cloud',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'!5Truncated',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'!5Spherical Cloud',size=txt_lgn_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(3),'!5Spherical Truncated',size=txt_lgn_sz,/NORMAL
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
; Figure 6a.
; Plot the vertical heating rate profile of the control cloud
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
erase
;
; Read in binary data from the NetCDF file 
;filename_1=init_rad_reg
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
; Allow the y-axis to be within .5 degrees of the curves
data_min=3600.*min([min(LW_heating_1),min(SW_heating_1),min(net_heating_1)])
data_max=3600.*max([max(LW_heating_1),max(SW_heating_1),max(net_heating_1)])
data_min=(data_min-.5) - (data_min mod .5)
data_max=(data_max+.5) - (data_max mod .5)
;
plot, $
	net_heating_1(1:num_layer_1)*3600.0, $
	altitude_1(1:num_layer_1)/1000.0, $
	tit='!5 Heating Rates in Control Cloud', $
	xtit='Heating Rate (K-hr!E-1!N)', $
	ytit='Altitude !8z!5 (km)', $
	xrange=[data_min,data_max], $
;	yrange=[8.,18.], $
	ystyle=1, $
	xstyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	SW_heating_1(1:num_layer_1)*3600.0, $
	altitude_1(1:num_layer_1)/1000.0, $
	thick=3.0, $
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
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
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
;filename_2=init_rad_sphere
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
; Allow the y-axis to be within 10% of the curves
data_min=min([min(LW_2_1_compare),min(SW_2_1_compare),min(net_2_1_compare)])
data_max=max([max(LW_2_1_compare),max(SW_2_1_compare),max(net_2_1_compare)])
data_min=(data_min-10.) - (data_min mod 10.)
data_max=(data_max+10.) - (data_max mod 10.)
;
plot, $
	net_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	tit='!5 Shape Effect on Heating Rates', $
;	xtit='Percent Change = 100*(Spheres - Columns)/!10"!5Columns!10"!5', $
	xtit='Percent Change (%)', $
;	xtit='Spheres - Columns (K hr!E-1!N)', $
	ytit='Altitude !8z!5 (km)', $
	xrange=[data_min,data_max], $
;	yrange=[0.1,0.6], $
	xstyle=1, $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	SW_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	thick=3.0, $
	linestyle=1
;
oplot,	$
	LW_2_1_compare(1:num_layer_2), $
	altitude_2(1:num_layer_2)/1000.0, $
	thick=1.0, $
	linestyle=2
;
;ln_lgn_x1=.22
ln_lgn_x1=.3
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
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
;filename_3=init_rad_truncate
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
; Allow the y-axis to be within 10% of the curves
data_min=min([min(LW_3_1_compare),min(SW_3_1_compare),min(net_3_1_compare)])
data_max=max([max(LW_3_1_compare),max(SW_3_1_compare),max(net_3_1_compare)])
data_min=(data_min-10.) - (data_min mod 10.)
data_max=(data_max+10.) - (data_max mod 10.)
;
plot, $
	net_3_1_compare(1:num_layer_3), $
	altitude_3(1:num_layer_3)/1000.0, $
	tit='!5 Truncation Effect on Heating Rates', $
;	xtit='Percent Change = 100*(Truncated - Control)/!10"!5Control!10"!5', $
	xtit='Percent Change (%)', $
;	xtit='Truncated - Control (K hr!E-1!N)', $
	ytit='Altitude !8z!5 (km)', $
	xrange=[data_min,data_max], $
;	yrange=[0.1,0.6], $
	xstyle=1, $
	ystyle=1, $
	thick=3.0, $
	charsize=2.0, $
	linestyle=0
;
oplot,	$
	SW_3_1_compare(1:num_layer_3), $
	altitude_3(1:num_layer_3)/1000.0, $
	thick=3.0, $
	linestyle=1
;
oplot,	$
	LW_3_1_compare(1:num_layer_3), $
	altitude_3(1:num_layer_3)/1000.0, $
	thick=1.0, $
	linestyle=2
;
ln_lgn_x1=.3
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=.4
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
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
; Figure 5agu
; Plot the time evolution of radiative forcing for the two main cases.
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
	tit='!5Diurnally Averaged Cloud Forcing', $
	xtit='Elapsed Time (min.)', $
	ytit='Forcing (W-m!E-2!N)', $
	yrange=[-150.0,250.0], $
	xstyle=1, $
	ystyle=1, $
	thick=3.0, $
;	xmargin=[0,0], $ 
;	ymargin=[0,0], $ 
	charsize=2.0, $
	linestyle=0
;
;oplot,	$
;	time_array_2(0:num_step_2-1)/60.0, $
;	net_forcing_2(0:num_step_2-1), $
;	thick=3.0, $
;	linestyle=0
;
oplot,	$
	time_array_3(0:num_step_3-1)/60.0, $
	net_forcing_3(0:num_step_3-1), $
	thick=3.0, $
	linestyle=0
;
; Now overplot the SW forcings
oplot,	$
	time_array_1(0:num_step_1-1)/60.0, $
	SW_forcing_1(0:num_step_1-1), $
	thick=3.0, $
	linestyle=1
;
;oplot,	$
;	time_array_2(0:num_step_2-1)/60.0, $
;	SW_forcing_2(0:num_step_2-1), $
;	thick=3.0, $
;	linestyle=1
;
oplot,	$
	time_array_3(0:num_step_3-1)/60.0, $
	SW_forcing_3(0:num_step_3-1), $
	thick=3.0, $
	linestyle=1
;
; Now overplot the LW forcings
oplot,	$
	time_array_1(0:num_step_1-1)/60.0, $
	LW_forcing_1(0:num_step_1-1), $
	thick=3.0, $
	linestyle=2
;
;oplot,	$
;	time_array_2(0:num_step_2-1)/60.0, $
;	LW_forcing_2(0:num_step_2-1), $
;	thick=3.0, $
;	linestyle=2
;
oplot,	$
	time_array_3(0:num_step_3-1)/60.0, $
	LW_forcing_3(0:num_step_3-1), $
	thick=3.0, $
	linestyle=2
;
ln_lgn_x1=.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.5
lgn_y_top=0.8
lgn_dy=0.1
lgn_y=[lgn_y_top,lgn_y_top-lgn_dy,lgn_y_top-2.*lgn_dy, $
	lgn_y_top-3.*lgn_dy,lgn_y_top-4.*lgn_dy]
;
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=1,thick=3.0,/NORMAL
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(2)+0.013,linestyle=2,thick=3.0,/NORMAL
;
;xyouts,txt_lgn_x,lgn_y(0),'!5Net',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(1),'!5SWCF',size=txt_lgn_sz,/NORMAL
;xyouts,txt_lgn_x,lgn_y(2),'!5LWCF',size=txt_lgn_sz,/NORMAL
;
txt_lgn_x=10
lgn_y=[125,35,75, $
	-55,-126,-120, $
	200,165]
;
xyouts,txt_lgn_x,lgn_y(0),'!5Truncated NetCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
;xyouts,txt_lgn_x,lgn_y(1),'!5Spherical NetCF',size=txt_lgn_sz,/DATA, $
;	ORIENTATION=-1
xyouts,txt_lgn_x,lgn_y(2),'!5Control NetCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1 
;
xyouts,txt_lgn_x,lgn_y(3),'!5Truncated SWCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
;xyouts,txt_lgn_x,lgn_y(4),'!5Spherical SWCF',size=txt_lgn_sz,/DATA, $
;	ORIENTATION=-1
xyouts,txt_lgn_x,lgn_y(5),'!5Control SWCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
;
xyouts,txt_lgn_x,lgn_y(6),'!5Control LWCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
xyouts,txt_lgn_x,lgn_y(7),'!5Truncated LWCF',size=txt_lgn_sz,/DATA, $
	ORIENTATION=-1
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
; End of Figure 5agu commands
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

