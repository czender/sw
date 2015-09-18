; $Id$

; Purpose: Offline analysis routines to look at aerosol data

@~/idl/ibp_fnc.pro
@~/idl/ibp_clr.pro
@~/idl/aspect.pro ; Needed for polar plot aspect ratios

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure aer_gph: Graphs various aerosol optical properties from mie output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro aer_gph, $
	fl_nm=fl_nm, $
	fl_ctl=fl_ctl, $ ; [sng] Control file for differences
	fl_prn_sng=fl_prn_sng, $ ; [sng] String to add to name of graphic file
	fld_nbr=fld_nbr, $ ; [nbr] Number of fields to plot [1--3]=[abs,sct,ext]=
	mie=mie, $
	abs=abs, $ ; [flg] Scale to absorption rather than extinction
	hrz=hrz, $ ; [flg] High resolution plotting
	dff=dff, $ ; [flg] Plot difference in absorption coefficient
	ffc=ffc, $ ; [flg] Plot effective phase function
	dbg=dbg, $
	info=info, $
	abc_max=abc_max, $ ; [um] Maximum abscissa
	abc_min=abc_min, $ ; [um] Minimum abscissa
	ord_max=ord_max, $ ; [m2 kg-1] Maximum ordinate
	ord_min=ord_min, $ ; [m2 kg-1] Minimum ordinate
	xlg_flg=xlg_flg, $ ; [flg] Plot abscissa on logarithmic axis
	ylg_flg=ylg_flg, $ ; [flg] Plot ordinate on logarithmic axis
	mrg_lft=mrg_lft, $ ; [frc] Left margin, 10 is default
	prn=prn
; Usage:
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_h2o_lqd_rds_swa_10_mie16.nc',abs=1,dff=1
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/ppr_ZeT05/aer_h2o_lqd_rds_swa_10_mie16mlrz.nc',abs=1
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/mie/mie.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_h2o_lqd_rds_swa_10_lrz.nc',abc_min=0.5,abc_max=1.5
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_h2o_lqd_rds_swa_10_mie16.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/mie/mie_1000.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/mie/mie_1000.nc',ord_max=1.0e-2
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/mie/mie_hrz_1000.nc',hrz=1,prn=0
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/mie/mie_1um_cld_1000000.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/mie/aer_saharan_dust_01_foo.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/dst/mie/aer_saharan_dust_01_CCM_LW.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/dst/mie/aer_saharan_dust_01_CCM_SW.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/dst/mie/aer_saharan_dust_01_SWNB.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/dst/mie/aer_afghan_dust_01_SWNB.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/dst/mie/aer_aeronet_Bnz_01_1000.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/dst/mie/aer_saharan_dust_02_CCM_SW.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/dst/mie/aer_saharan_dust_03_CCM_SW.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/dst/mie/aer_saharan_dust_04_CCM_SW.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_mineral_dust.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_saharan_dust.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_ncsu.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_YZS99.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_H2O_lqd_SW.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_h2so4_215K.nc'
; aer_gph,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_Fe2O3.nc'
@~/idl/ibp_clr.com
!p.multi=0

if n_elements(mie) eq 0 then mie=0
if n_elements(dff) eq 0 then dff=0 ; [flg] Plot difference in absorption coefficient
if n_elements(ffc) eq 0 then ffc=0 ; [flg] Plot effective phase function
if n_elements(abc_max) eq 0 then abc_max=0.0 ; [um] Maximum abscissa
if n_elements(abc_min) eq 0 then abc_min=0.0 ; [um] Minimum abscissa
if n_elements(ord_max) eq 0 then ord_max=0.0 ; [m2 kg-1] Maximum ordinate
if n_elements(ord_min) eq 0 then ord_min=0.0 ; [m2 kg-1] Minimum ordinate
if n_elements(xlg_flg) eq 0 then xlg_flg=0 ; [flg] Plot abscissa on logarithmic axis
if n_elements(ylg_flg) eq 0 then ylg_flg=0 ; [flg] Plot ordinate on logarithmic axis
if n_elements(mrg_lft) eq 0 then mrg_lft=8 ; [frc] Left margin, 10 is default
if n_elements(dbg) eq 0 then dbg=0
if n_elements(abs) eq 0 then abs=0 ; [flg] Scale to absorption rather than extinction
if n_elements(hrz) eq 0 then hrz=0 ; [flg] High resolution plotting
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(fld_nbr) eq 0 then fld_nbr=3 ; [nbr] Number of fields to plot [1--3]=[abs,sct,ext]=
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/mie/mie.nc'
if n_elements(fl_prn_sng) eq 0 then fl_prn_sng='' ; [sng] String to add to name of graphic file
if n_elements(fl_ctl) eq 0 then fl_ctl=getenv('DATA')+'/aca/aer_h2o_lqd_rds_swa_10_lrz.nc' ; [sng] Control file for differences
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22

;if xlg_flg then x_sty=1 else x_sty=0 ; [idx] Plot box exactly circumscribes abscissa
;if ylg_flg then y_sty=1 else y_sty=0 ; [idx] Plot box exactly circumscribes ordinate

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=3.0
;	y_sz=6.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

fl_idx=0
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
dim_id=ncdf_dimid(nc_id,'wvl')
ncdf_diminq,nc_id,dim_id,foo,wvl_nbr
if dbg then print,'Current number of wavelengths = ',auto_sng(wvl_nbr,0)
; netCDF character strings are returned as byte arrays to IDL, convert with string(att_val)
ncdf_attget,nc_id,/GLOBAL,'aerosol_long_name',aer_lng_nm
ncdf_varget,nc_id,'wvl',wvl
ncdf_varget,nc_id,'abs_cff_mss',abs_cff_mss
ncdf_varget,nc_id,'abs_fsh_ffc',abs_fsh_ffc
ncdf_varget,nc_id,'ang_xpn',ang_xpn
ncdf_varget,nc_id,'asm_prm',asm_prm
ncdf_varget,nc_id,'dmt_ctr',dmt_ctr
ncdf_varget,nc_id,'dmt_nmr',dmt_nmr
ncdf_varget,nc_id,'ext_cff_mss',ext_cff_mss
ncdf_varget,nc_id,'sca_cff_mss',sca_cff_mss
ncdf_varget,nc_id,'ss_alb',ss_alb
ncdf_varget,nc_id,'ss_co_alb',ss_co_alb
; 20011111 Temporarily comment out phase function graphs while mie is broken
;ncdf_varget,nc_id,'psd_gsd_anl',psd_gsd_anl
;ncdf_varget,nc_id,'wvl_dbg',wvl_dbg
;ncdf_varget,nc_id,'sz_dbg',sz_dbg
;dmt_dbg=2.0*sz_dbg
ncdf_varget,nc_id,'ngl_dgr',ngl_dgr
ncdf_varget,nc_id,'phz_fnc',phz_fnc
ncdf_varget,nc_id,'phz_fnc_ffc',phz_fnc_ffc
;ncdf_varget,nc_id,'foo',foo
ncdf_close,nc_id

abc=wvl*1.0e6
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
foo=min((abc-abc_max)*(abc-abc_max),abc_max_idx)
foo=min((abc-abc_min)*(abc-abc_min),abc_min_idx)
abc=abc(abc_min_idx:abc_max_idx)
ord=ext_cff_mss(abc_min_idx:abc_max_idx)
ord_nm='ext_cff_mss'
if hrz or abs then begin
	ord=abs_cff_mss(abc_min_idx:abc_max_idx)
        ord_nm='abs_cff_mss'
;	ord=abs_fsh_ffc(abc_min_idx:abc_max_idx)
;        ord_nm='abs_fsh_ffc'
endif ; endif hrz or abs
foo=max(ord,ord_max_idx)
if ord_max eq 0.0 then ord_max=foo
foo=min(ord,ord_min_idx)
if ord_min eq 0.0 then ord_min=foo
print,'Ordinate = ',ord_nm
print,'Maximum ordinate = ',ord(ord_max_idx),' at ',abc(ord_max_idx),' um'
print,'Minimum ordinate = ',ord(ord_min_idx),' at ',abc(ord_min_idx),' um'
if hrz then begin
	print,'Offsetting abscissa to resonance center at ',abc(ord_max_idx),' um'
	print,'NOTE re-scaling abscissa to pm (ridiculously fine units!)'
	wvl_fst=abc(ord_max_idx)*1.0e-6
	abc=(abc-abc(ord_max_idx))*1.0e6
; Recompute abc_min and abc_max based on resonance neighborhood so rng_x comes out right
        abc_max=max(abc)
        abc_min=min(abc)
endif ; endif hrz
rng_x=[abc_min,abc_max]
;rng_x=[0.2,1.0]
;rng_x=[0.5055,0.5057] ; This range contains step jump with BoH83 data
;rng_x=[0.5000039,0.5000043] ; This range contains resonance with Wis79 data
rng_y=[ord_min,ord_max]
ord_max=0.0 ; Reset so value does not effect next plot
ord_min=0.0 ; Reset so value does not effect next plot
if prn then chr_sz=1.5 else chr_sz=1.0
ttl=string(aer_lng_nm)
;ttl='Mineral Dust'
x_ttl='!5Wavelength !7k!5 (!7l!5m)'
y_ttl='!5Specific Extinction !7w!5!Im!N (m!E2!N kg!E-1!N)'
if hrz then begin
	ttl=string(aer_lng_nm)+' Resonant Absorptance'
	x_ttl='!5Wavelength !7k!5 Offset in pm from '+auto_sng(wvl_fst*1.0e6,6)+' (!7l!5m)'
endif; endif hrz
if hrz or abs then begin
	y_ttl='!5Absorption !7w!5!Ia!N (m!E2!N kg!E-1!N)'
endif; endif hrz

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
;mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	if fl_prn_sng eq '' then fl_nm_out=getenv('DATA')+'/ps/aer_ext_cff_mss.eps' else fl_nm_out=getenv('DATA')+'/ps/aer_'+fl_prn_sng+'.eps'
	fl_psn=strpos(fl_nm,'aer_')
	nc_psn=strpos(fl_nm,'.nc')
	if fl_prn_sng eq '' and fl_psn gt -1 and nc_psn gt -1 then begin
		fl_nm_out=getenv('DATA')+'/ps/'+strmid(fl_nm,fl_psn,(nc_psn-fl_psn))+'_aer_ext_cff_mss.eps'
	endif; endif		
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot aer_gph
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=1,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
yll_idx=base_clr_idx+3
orn_idx=base_clr_idx+4
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(orn_idx)=255
g_curr(orn_idx)=128
b_curr(orn_idx)=0
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

lbl_sng=strarr(fld_nbr)
ln_sty=0.0*intarr(fld_nbr)
clr=intarr(fld_nbr)
clr(0)=clr_blk_idx
if fld_nbr ge 2 then clr(1)=grn_idx
if fld_nbr ge 3 then clr(2)=blu_idx
;clr(3)=yll_idx
;clr(4)=orn_idx
lbl_sng(0)='Absorption'
if fld_nbr ge 2 then lbl_sng(1)='Scattering'
if fld_nbr ge 3 then lbl_sng(2)='Extinction'
;lbl_sng(3)='dst_mpc_a_fix'
;lbl_sng(4)='dst_mpc_b_tphys'

if fld_nbr ge 1 then oplot,abc,abs_cff_mss(abc_min_idx:abc_max_idx),max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)
if fld_nbr ge 2 then oplot,abc,sca_cff_mss(abc_min_idx:abc_max_idx),max_value=1.0e20,thick=2.0,linestyle=0,color=clr(1)
if fld_nbr ge 3 then oplot,abc,ext_cff_mss(abc_min_idx:abc_max_idx),max_value=1.0e20,thick=2.0,linestyle=0,color=clr(2)

if info then begin
ln_lgn_x1=0.20
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

;xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5'+fl_nm,color=clr_blk_idx,size=chr_sz,/NORMAL

for fld_idx=0,fld_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fld_idx)+0.013,linestyle=ln_sty(fld_idx),color=clr(fld_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fld_idx),lbl_sng(fld_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

; Plot absorptance difference file
if dff then begin; 

nc_id=ncdf_open(fl_ctl)
dim_id=ncdf_dimid(nc_id,'wvl')
ncdf_diminq,nc_id,dim_id,foo,wvl_nbr_ctl
if wvl_nbr ne wvl_nbr_ctl then print,'ERROR: wavelength mismatch'
if dbg then print,'Control number of wavelengths = ',auto_sng(wvl_nbr,0)
ncdf_varget,nc_id,'abs_cff_mss',abs_cff_mss_ctl
ncdf_varget,nc_id,'ext_cff_mss',ext_cff_mss_ctl
ncdf_close,nc_id

abc=wvl*1.0e6
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
rng_x=[abc_min,abc_max]
foo=min((abc-abc_max)*(abc-abc_max),abc_max_idx)
foo=min((abc-abc_min)*(abc-abc_min),abc_min_idx)
abc=abc(abc_min_idx:abc_max_idx)
ord=100.0*(ext_cff_mss(abc_min_idx:abc_max_idx)-ext_cff_mss_ctl(abc_min_idx:abc_max_idx))/ext_cff_mss(abc_min_idx:abc_max_idx)
ord_nm='ext_cff_mss'
if hrz or abs then begin
        ord=100.0*(abs_cff_mss(abc_min_idx:abc_max_idx)-abs_cff_mss_ctl(abc_min_idx:abc_max_idx))/abs_cff_mss(abc_min_idx:abc_max_idx)
        ord_nm='abs_cff_mss'
endif ; endif hrz or abs
; Find RMS percent relative bias
ord_wgt=ord*1.0 ; Do not flux-weight yet
ord_rms=total(sqrt(ord*ord))/n_elements(ord)
print,'Control file: ',fl_ctl
print,'Experiment file: ',fl_nm
print,'RMS percent relative bias total(sqrt((100*(abs_cff_mss_xpt-abs_cff_mss_ctl)/abs_cff_mss_ctl))^2))/ord_nbr= ',ord_rms
foo=max(ord,ord_max_idx)
ord_max=foo
foo=min(ord,ord_min_idx)
ord_min=foo
print,'Ordinate = ',ord_nm
print,'Maximum ordinate = ',ord(ord_max_idx),' at ',abc(ord_max_idx),' um'
print,'Minimum ordinate = ',ord(ord_min_idx),' at ',abc(ord_min_idx),' um'
rng_y=[ord_min,ord_max]
;rng_y=[-1.0,2.0]
if prn then chr_sz=1.5 else chr_sz=1
ttl=string(aer_lng_nm)+' Absorption' 
x_ttl='!5Wavelength !7k!5 (!7l!5m)'
y_ttl='!5Percent Change in !7w!5!Ia!N'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
;mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/aer_'+ord_nm+'_err.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot abs_dff
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,color=clr_blk_idx

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

endif; endif dff

abc=wvl*1.0e6
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
rng_x=[abc_min,abc_max]
foo=min((abc-abc_max)*(abc-abc_max),abc_max_idx)
foo=min((abc-abc_min)*(abc-abc_min),abc_min_idx)
abc=abc(abc_min_idx:abc_max_idx)
ord=ang_xpn(abc_min_idx:abc_max_idx)
ord_nm='ang_xpn'
foo=max(ord,ord_max_idx)
ord_max=foo
foo=min(ord,ord_min_idx)
ord_min=foo
print,'Ordinate = ',ord_nm
print,'Maximum ordinate = ',ord(ord_max_idx),' at ',abc(ord_max_idx),' um'
print,'Minimum ordinate = ',ord(ord_min_idx),' at ',abc(ord_min_idx),' um'
rng_y=[-1.0,2.0]
if prn then chr_sz=2.0 else chr_sz=1
ttl=string(aer_lng_nm)
x_ttl='!5Wavelength !7k!5 (!7l!5m)'
y_ttl='!5Angstrom Exponent'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/aer_ang_xpn.eps'
	fl_psn=strpos(fl_nm,'aer_')
	nc_psn=strpos(fl_nm,'.nc')
	if fl_psn gt -1 and nc_psn gt -1 then begin
		fl_nm_out=getenv('DATA')+'/ps/'+strmid(fl_nm,fl_psn,(nc_psn-fl_psn))+'_ang_xpn.eps'
	endif; endif		
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot ang_xpn
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=1,thick=2.0,charsize=chr_sz,linestyle=0,color=clr_blk_idx

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

abc=wvl*1.0e6
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
rng_x=[abc_min,abc_max]
;rng_x=[0.2,1.0]
foo=min((abc-abc_max)*(abc-abc_max),abc_max_idx)
foo=min((abc-abc_min)*(abc-abc_min),abc_min_idx)
abc=abc(abc_min_idx:abc_max_idx)
ord=ss_alb(abc_min_idx:abc_max_idx)
ord_nm='ss_alb'
foo=max(ord,ord_max_idx)
ord_max=foo
foo=min(ord,ord_min_idx)
ord_min=foo
print,'Ordinate = ',ord_nm
print,'Maximum ordinate = ',ord(ord_max_idx),' at ',abc(ord_max_idx),' um'
print,'Minimum ordinate = ',ord(ord_min_idx),' at ',abc(ord_min_idx),' um'
;if ord_max eq 0.0 then ord_max=max(ord)
;if ord_min eq 0.0 then ord_min=min(ord)
rng_y=[ord_min,ord_max]
;rng_y=[0.0,1.0]
if prn then chr_sz=2.0 else chr_sz=1
ttl=string(aer_lng_nm)
;ttl='Saharan, Afghan, Bnz'
x_ttl='!5Wavelength !7k!5 (!7l!5m)'
y_ttl='!5Single Scatter Albedo !7x!5'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/aer_ss_alb.eps'
	fl_psn=strpos(fl_nm,'aer_')
	nc_psn=strpos(fl_nm,'.nc')
	if fl_psn gt -1 and nc_psn gt -1 then begin
		fl_nm_out=getenv('DATA')+'/ps/'+strmid(fl_nm,fl_psn,(nc_psn-fl_psn))+'_ss_alb.eps'
	endif; endif		
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot ss_alb
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=1,thick=2.0,charsize=chr_sz,linestyle=0,color=clr_blk_idx

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

abc=wvl*1.0e6
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
rng_x=[abc_min,abc_max]
;rng_x=[0.2,1.0]
foo=min((abc-abc_max)*(abc-abc_max),abc_max_idx)
foo=min((abc-abc_min)*(abc-abc_min),abc_min_idx)
abc=abc(abc_min_idx:abc_max_idx)
ord=ss_co_alb(abc_min_idx:abc_max_idx)
ord_nm='ss_co_alb'
foo=max(ord,ord_max_idx)
ord_max=foo
foo=min(ord,ord_min_idx)
ord_min=foo
print,'Ordinate = ',ord_nm
print,'Maximum ordinate = ',ord(ord_max_idx),' at ',abc(ord_max_idx),' um'
print,'Minimum ordinate = ',ord(ord_min_idx),' at ',abc(ord_min_idx),' um'
if ord_max eq 0.0 then ord_max=max(ord)
if ord_min eq 0.0 then ord_min=min(ord)
rng_y=[ord_min,ord_max]
;rng_y=[0.0,1.0]
if prn then chr_sz=2.0 else chr_sz=1
ttl=string(aer_lng_nm)
;ttl='Saharan, Afghan, Bnz'
x_ttl='!5Wavelength !7k!5 (!7l!5m)'
y_ttl='!5Single Scatter Co-albedo 1-!7x!5'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/aer_ss_co_alb.eps'
	fl_psn=strpos(fl_nm,'aer_')
	nc_psn=strpos(fl_nm,'.nc')
	if fl_psn gt -1 and nc_psn gt -1 then begin
		fl_nm_out=getenv('DATA')+'/ps/'+strmid(fl_nm,fl_psn,(nc_psn-fl_psn))+'_ss_co_alb.eps'
	endif; endif		
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot ss_co_alb
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=1,thick=2.0,charsize=chr_sz,linestyle=0,color=clr_blk_idx

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

abc=wvl*1.0e6
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
rng_x=[abc_min,abc_max]
;rng_x=[0.2,1.0]
foo=min((abc-abc_max)*(abc-abc_max),abc_max_idx)
foo=min((abc-abc_min)*(abc-abc_min),abc_min_idx)
abc=abc(abc_min_idx:abc_max_idx)
ord=asm_prm(abc_min_idx:abc_max_idx)
ord_nm='asm_prm'
foo=max(ord,ord_max_idx)
ord_max=foo
foo=min(ord,ord_min_idx)
ord_min=foo
print,'Ordinate = ',ord_nm
print,'Maximum ordinate = ',ord(ord_max_idx),' at ',abc(ord_max_idx),' um'
print,'Minimum ordinate = ',ord(ord_min_idx),' at ',abc(ord_min_idx),' um'
if ord_max eq 0.0 then ord_max=max(ord)
if ord_min eq 0.0 then ord_min=min(ord)
rng_y=[ord_min,ord_max]
;rng_y=[0.0,1.0]
if prn then chr_sz=2.0 else chr_sz=1
ttl=string(aer_lng_nm)
;ttl='Saharan, Afghan, Bnz'
x_ttl='!5Wavelength !7k!5 (!7l!5m)'
y_ttl='!5Asymmetry Parameter !8g!5'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/aer_asm_prm.eps'
	fl_psn=strpos(fl_nm,'aer_')
	nc_psn=strpos(fl_nm,'.nc')
	if fl_psn gt -1 and nc_psn gt -1 then begin
		fl_nm_out=getenv('DATA')+'/ps/'+strmid(fl_nm,fl_psn,(nc_psn-fl_psn))+'_asm_prm.eps'
	endif; endif		
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot asm_prm
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=1,thick=2.0,charsize=chr_sz,linestyle=0,color=clr_blk_idx

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

if(ffc) then abc=phz_fnc_ffc else abc=phz_fnc
abc=[reverse(abc(1:n_elements(abc)-1)),abc] ; Add symmetric data for 360 degree plots
if (max(abc)/min(abc)) gt 10 then xpn_flg=1 else xpn_flg=0; Plot data by exponent
if xpn_flg then begin
	scl_xpn=abs(floor(min(alog10(abc))))+1
	print,'scl_xpn = ',scl_xpn
	scl=10.0^scl_xpn
	abc=alog10(scl*abc)
endif; endif xpn_flg
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
rng_x=[abc_min,abc_max]
rng_x=[-abc_max,abc_max]
ord=ngl_dgr*!pi/180.0
ord=[-1.0*reverse(ord(1:n_elements(ord)-1)),ord]
if ord_max eq 0.0 then ord_max=max(ord)
if ord_min eq 0.0 then ord_min=min(ord)
rng_y=[ord_min,ord_max]
rng_y=rng_x
if prn then chr_sz=2.0 else chr_sz=1
ttl=string(aer_lng_nm)
x_ttl=''
y_ttl='!5Scattering Angle !7H!5 (radians)'

if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/phz_fnc_plr.ps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/ps,/color
endif; endif prn

; polar plot phz_fnc
print,'abc = ',abc
print,'abc_min = ',abc_min,', abc_max = ',abc_max
xlg_flg=0
ylg_flg=0
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,tit=ttl,xtit=x_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=5,ystyle=5,thick=2.0,charsize=chr_sz,linestyle=0,color=clr_blk_idx,/nodata,/polar,psym=2,position=aspect(1.0)
oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx,/polar
if xpn_flg then tck_fmt_fnc='plr_tck_fmt' else tck_fmt_fnc='pst_tck_fmt'
axis,0.0,0.0,xaxis=0,xlog=0,xtick_get=x_tck,xstyle=1,charsize=0.75*chr_sz,color=clr_blk_idx,xtickformat=tck_fmt_fnc ; IDL 3.6 UG p. 14-24
axis,0.0,0.0,yaxis=0,ylog=0,ytick_get=y_tck,ystyle=1,charsize=0.75*chr_sz,color=clr_blk_idx,ytickformat='null_fmt' ; IDL 3.6 UG p. 14-24

; Draw circles
max_dat=max(abc)
for idx=0,n_elements(x_tck)-1 do begin
	crc_val=circle(0,0,x_tck(idx))
	if idx ne 0 then plots,crc_val,color=2
endfor; end loop over tck

if info then begin
ln_lgn_x1=0.07
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.77
lgn_dy=0.05
lgn_y=lgn_y_top-indgen(6)*lgn_dy
if xpn_flg then xyouts,txt_lgn_x+0.45,0.22,'!5log!I10!N[10!E'+auto_sng(scl_xpn,0)+'!N!8p!5(!7H!5)!5]',color=clr_blk_idx,size=0.75*chr_sz,/NORMAL
xyouts,txt_lgn_x+0.45,lgn_y(0),'!7k!5 = '+auto_sng(wvl_dbg*1.0e6,2)+' !7l!5m',color=clr_blk_idx,size=0.75*chr_sz,/NORMAL
if ffc then xyouts,txt_lgn_x,lgn_y(0),'!8D!5!Inmr!N = '+auto_sng(dmt_nmr*1.0e6,2)+' !7l!5m',color=clr_blk_idx,size=0.75*chr_sz,/NORMAL
if ffc then xyouts,txt_lgn_x,lgn_y(1),'!7r!5 = '+auto_sng(psd_gsd_anl,2),color=clr_blk_idx,size=0.75*chr_sz,/NORMAL
if not ffc then xyouts,txt_lgn_x,lgn_y(0),'!8D!5 = '+auto_sng(dmt_dbg*1.0e6,2)+' !7l!5m',color=clr_blk_idx,size=0.75*chr_sz,/NORMAL
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

abc=phz_fnc
abc=[reverse(abc(1:n_elements(abc)-1)),abc]
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
rng_x=[abc_min,abc_max]
rng_x=[10^(1.0*floor(alog10(abc_min))),abc_max]
ord=ngl_dgr
ord=[-1.0*reverse(ord(1:n_elements(ord)-1)),ord]
if ord_max eq 0.0 then ord_max=max(ord)
if ord_min eq 0.0 then ord_min=min(ord)
rng_y=[ord_min,ord_max]
if prn then chr_sz=2.0 else chr_sz=1
ttl=string(aer_lng_nm)+'!5 Phase function'
x_ttl='!5Phase function !8p!5(!7H!5)'
y_ttl='!5Scattering Angle !7H!5 (degrees)'

if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/phz_fnc_xy.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot phz_fnc
xlg_flg=1
ylg_flg=0
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=0,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,color=clr_blk_idx,/nodata
oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end aer_gph()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure aer_gph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure cnc_mss_gph: Plots time series of dust mass concentrations 
; from AEROCE data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cnc_mss_gph, $
	fl_nm=fl_nm, $
	info=info, $
	prn=prn
; Usage:
; cnc_mss_gph,fl_nm=getenv('DATA')+'/aeroce/aeroce_dly.nc'
; cnc_mss_gph,fl_nm=getenv('DATA')+'/aeroce/aeroce_19890101_19961231.nc'
; cnc_mss_gph,fl_nm=getenv('DATA')+'/aeroce/aeroce_19900101_19901231.nc'
@~/idl/ibp_clr.com
!p.multi=0

if n_elements(dbg) eq 0 then dbg=0
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/aeroce/aeroce_dly.nc'
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=4.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

fl_idx=0
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
dim_id=ncdf_dimid(nc_id,'time')
ncdf_diminq,nc_id,dim_id,foo,time_nbr
if dbg then print,'Current number of time samples is '+time_nbr
; netCDF character strings are returned as byte arrays to IDL, convert with string(att_val)
cnc_mss_dst_Brm_id=ncdf_varid(nc_id,'cnc_mss_dst_Brm')
ncdf_attget,nc_id,cnc_mss_dst_Brm_id,'long_name',cnc_long_name_Brm
cnc_mss_dst_Brb_id=ncdf_varid(nc_id,'cnc_mss_dst_Brb')
ncdf_attget,nc_id,cnc_mss_dst_Brb_id,'long_name',cnc_long_name_Brb
cnc_mss_dst_Izn_id=ncdf_varid(nc_id,'cnc_mss_dst_Izn')
ncdf_attget,nc_id,cnc_mss_dst_Izn_id,'long_name',cnc_long_name_Izn
time_id=ncdf_varid(nc_id,'time')
ncdf_attget,nc_id,time_id,'long_name',time_long_name
doy_id=ncdf_varid(nc_id,'doy')
ncdf_attget,nc_id,doy_id,'long_name',doy_long_name
mth_id=ncdf_varid(nc_id,'mth')
ncdf_attget,nc_id,mth_id,'long_name',mth_long_name

ncdf_varget,nc_id,'time',time
ncdf_varget,nc_id,'doy',doy
ncdf_varget,nc_id,'mth',mth
ncdf_varget,nc_id,'cnc_mss_dst_Brm',cnc_mss_dst_Brm
ncdf_varget,nc_id,'cnc_mss_dst_Brb',cnc_mss_dst_Brb
ncdf_varget,nc_id,'cnc_mss_dst_Izn',cnc_mss_dst_Izn
;ncdf_varget,nc_id,'foo',foo
ncdf_close,nc_id

abc=mth
ord=cnc_mss_dst_Brm
xcl_mss_dat,abc,ord,max=1.0e20,min=0.0,nbr=nbr
abc_min=min(abc)
abc_max=max(abc)
rng_x=[abc_min,abc_max]
ord_min=min(ord)
ord_max=max(ord)
;rng_y=[ord_min,ord_max]
rng_y=[0,ord_max]
if prn then chr_sz=2.0 else chr_sz=1
ttl='!5Mineral Dust Mass Concentration'
print,size(ttl)
;ttl='Mineral Dust'
x_ttl='!5'+string(mth_long_name)
y_ttl='!5Concentration !8q!5 (kg m!E-3!N)'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/cnc_mss_gph.eps'
	fl_psn=strpos(fl_nm,'aer_')
	nc_psn=strpos(fl_nm,'.nc')
	if fl_psn gt -1 and nc_psn gt -1 then begin
		fl_nm_out=getenv('DATA')+'/ps/'+strmid(fl_nm,fl_psn,(nc_psn-fl_psn))+'.eps'
	endif; endif		
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot cnc_mss
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=0,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
yll_idx=base_clr_idx+3
orn_idx=base_clr_idx+4
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(orn_idx)=255
g_curr(orn_idx)=128
b_curr(orn_idx)=0
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

fld_nbr=3
lbl_sng=strarr(fld_nbr)
ln_sty=0.0*intarr(fld_nbr)
clr=intarr(fld_nbr)
clr(0)=clr_blk_idx
clr(1)=grn_idx
clr(2)=blu_idx
;clr(3)=yll_idx
;clr(4)=orn_idx
lbl_sng(0)='Bermuda'
lbl_sng(1)='Barbados'
lbl_sng(2)='Izania'
;lbl_sng(3)='dst_mpc_a_fix'
;lbl_sng(4)='dst_mpc_b_tphys'

oplot,abc,cnc_mss_dst_Brm,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)
oplot,abc,cnc_mss_dst_Brb,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(1)
oplot,abc,cnc_mss_dst_Izn,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(2)

if info then begin
ln_lgn_x1=0.50
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

;xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5'+fl_nm,color=clr_blk_idx,size=chr_sz,/NORMAL

for fl_idx=0,fld_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end cnc_mss_gph()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure cnc_mss_gph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure cnv_gph_mie: Graphs various aerosol optical convergence
; properties from mie files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cnv_gph_mie, $
	fl_nm=fl_nm, $
	mie=mie, $
	abs=abs, $ ; [flg] Scale to absorption rather than extinction
	hrz=hrz, $ ; [flg] High resolution plotting
	ffc=ffc, $ ; Plot effective phase function
	dbg=dbg, $
	info=info, $
	abc_max=abc_max, $ ; [m2 kg-1] Maximum abscissa
	abc_min=abc_min, $ ; [m2 kg-1] Minimum abscissa
	ord_max=ord_max, $ ; [m2 kg-1] Maximum ordinate
	ord_min=ord_min, $ ; [m2 kg-1] Minimum ordinate
	xlg_flg=xlg_flg, $ ; [flg] Plot abscissa on logarithmic axis
	ylg_flg=ylg_flg, $ ; [flg] Plot ordinate on logarithmic axis
	prn=prn
; Usage:
; cnv_gph_mie,dbg=1,fl_nm=getenv('DATA')+'/mie/mie.nc'
; cnv_gph_mie,dbg=1,fl_nm=getenv('DATA')+'/aca/aer_h2o_lqd_rds_swa_10_lrz.nc',abc_min=0.5,abc_max=1.5
; cnv_gph_mie,dbg=1,fl_nm=getenv('DATA')+'/ppr_ZeT05/aer_h2o_lqd_rds_swa_10_avg_mie10_mie16.nc'
@~/idl/ibp_clr.com
!p.multi=0

if n_elements(mie) eq 0 then mie=0
if n_elements(ffc) eq 0 then ffc=0 ; Plot effective phase function
if n_elements(abc_max) eq 0 then abc_max=0.0 ; [m2 kg-1] Maximum abscissa
if n_elements(abc_min) eq 0 then abc_min=0.0 ; [m2 kg-1] Minimum abscissa
if n_elements(ord_max) eq 0 then ord_max=0.0 ; [m2 kg-1] Maximum ordinate
if n_elements(ord_min) eq 0 then ord_min=0.0 ; [m2 kg-1] Minimum ordinate
if n_elements(xlg_flg) eq 0 then xlg_flg=0 ; [flg] Plot abscissa on logarithmic axis
if n_elements(ylg_flg) eq 0 then ylg_flg=0 ; [flg] Plot ordinate on logarithmic axis
if n_elements(dbg) eq 0 then dbg=0
if n_elements(abs) eq 0 then abs=0 ; [flg] Scale to absorption rather than extinction
if n_elements(hrz) eq 0 then hrz=0 ; [flg] High resolution plotting
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/mie/mie.nc'
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=3.0
;	y_sz=6.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

fl_idx=0
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
dim_id=ncdf_dimid(nc_id,'rsn')
ncdf_diminq,nc_id,dim_id,foo,rsn_nbr
if dbg then print,'Current number of resolutions = ',auto_sng(rsn_nbr,0)
; netCDF character strings are returned as byte arrays to IDL, convert with string(att_val)
ncdf_attget,nc_id,/GLOBAL,'aerosol_long_name',aer_lng_nm
ncdf_varget,nc_id,'abs_cff_mss',abs_cff_mss
ncdf_varget,nc_id,'abs_fsh_ffc',abs_fsh_ffc
ncdf_varget,nc_id,'ang_xpn',ang_xpn
ncdf_varget,nc_id,'asm_prm',asm_prm
ncdf_varget,nc_id,'dmt_ctr',dmt_ctr
ncdf_varget,nc_id,'dmt_nmr',dmt_nmr
ncdf_varget,nc_id,'ext_cff_mss',ext_cff_mss
ncdf_varget,nc_id,'sca_cff_mss',sca_cff_mss
ncdf_varget,nc_id,'ss_alb',ss_alb
ncdf_varget,nc_id,'ss_co_alb',ss_co_alb
ncdf_varget,nc_id,'sz_prm_rsn_avg',sz_prm_rsn_avg
ncdf_varget,nc_id,'wvl_ctr',wvl_ctr
ncdf_close,nc_id

abc=[1,10,100,1000,10000,100000,1000000] ; [# (0.01 um)-1] Bands per 0.01 um
abc=abc*100 ; [# (0.01 um)-1]->[# um-1] Bands per um
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
ord=abs_cff_mss
; 20041127 Get from post-processing ord_rms in aer_gph in dff mode
ord=[4.7614208,2.2497703,1.1436235,1.0565009,0.29511564,0.10467634,0.0] ; [pct] RMS percent relative bias
foo=max(ord,ord_max_idx)
if ord_max eq 0.0 then ord_max=foo
foo=min(ord,ord_min_idx)
if ord_min eq 0.0 then ord_min=foo
ord_nm='abs_cff_mss'
rng_x=[abc_min,abc_max]
rng_y=[ord_min,ord_max]
ord_max=0.0 ; Reset so value does not effect next plot
ord_min=0.0 ; Reset so value does not effect next plot
if prn then chr_sz=1.5 else chr_sz=1.0
ttl=string(aer_lng_nm)
;ttl='Mineral Dust'
;x_ttl='!5Resolution !8R!5 (# !7l!5m!E-1!N)'
x_ttl='!5Quadrature Density !8N!5 (# !7l!5m!E-1!N)'
;y_ttl='!5Specific Absorption !7w!5!Im!N (m!E2!N kg!E-1!N)'
y_ttl='!5RMS Absorption Bias (%)'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/aer_abs_cff_mss_cnv.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot cnv_gph_mie
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=1,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
yll_idx=base_clr_idx+3
orn_idx=base_clr_idx+4
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(orn_idx)=255
g_curr(orn_idx)=128
b_curr(orn_idx)=0
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

fld_nbr=1
lbl_sng=strarr(fld_nbr)
ln_sty=0.0*intarr(fld_nbr)
clr=intarr(fld_nbr)
clr(0)=clr_blk_idx
;clr(1)=grn_idx
;clr(2)=blu_idx
;clr(3)=yll_idx
;clr(4)=orn_idx
lbl_sng(0)='Absorption'
;lbl_sng(1)='Absorption'
;lbl_sng(2)='Scattering'

oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)
;oplot,abc,ext_cff_mss(abc_min_idx:abc_max_idx),max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)

if info then begin
ln_lgn_x1=0.50
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

;xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5'+fl_nm,color=clr_blk_idx,size=chr_sz,/NORMAL

for fl_idx=0,fld_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end cnv_gph_mie()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure cnv_gph_mie
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure cnv_gph_swnb2: Graphs various aerosol optical convergence
; properties from mie files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cnv_gph_swnb2, $
	fl_nm=fl_nm, $
	mie=mie, $
	abs=abs, $ ; [flg] Scale to absorption rather than extinction
	hrz=hrz, $ ; [flg] High resolution plotting
	ffc=ffc, $ ; Plot effective phase function
	dbg=dbg, $
	info=info, $
	abc_max=abc_max, $ ; [m2 kg-1] Maximum abscissa
	abc_min=abc_min, $ ; [m2 kg-1] Minimum abscissa
	ord_max=ord_max, $ ; [m2 kg-1] Maximum ordinate
	ord_min=ord_min, $ ; [m2 kg-1] Minimum ordinate
	xlg_flg=xlg_flg, $ ; [flg] Plot abscissa on logarithmic axis
	ylg_flg=ylg_flg, $ ; [flg] Plot ordinate on logarithmic axis
	mrg_lft=mrg_lft, $ ; [frc] Left margin, 10 is default
	prn=prn
; Usage:
; cnv_gph_swnb2,dbg=1,fl_nm=getenv('DATA')+'/aca/swnb.nc'
; cnv_gph_swnb2,dbg=1,fl_nm=getenv('DATA')+'/ppr_ZeT05/swnb2_aer_h2o_lqd_rds_swa_10_mie10_mie16.nc'
@~/idl/ibp_clr.com
!p.multi=0

if n_elements(mie) eq 0 then mie=0
if n_elements(ffc) eq 0 then ffc=0 ; Plot effective phase function
if n_elements(abc_max) eq 0 then abc_max=0.0 ; [m2 kg-1] Maximum abscissa
if n_elements(abc_min) eq 0 then abc_min=0.0 ; [m2 kg-1] Minimum abscissa
if n_elements(ord_max) eq 0 then ord_max=0.0 ; [m2 kg-1] Maximum ordinate
if n_elements(ord_min) eq 0 then ord_min=0.0 ; [m2 kg-1] Minimum ordinate
if n_elements(xlg_flg) eq 0 then xlg_flg=0 ; [flg] Plot abscissa on logarithmic axis
if n_elements(ylg_flg) eq 0 then ylg_flg=0 ; [flg] Plot ordinate on logarithmic axis
if n_elements(mrg_lft) eq 0 then mrg_lft=8 ; [frc] Left margin, 10 is default
if n_elements(dbg) eq 0 then dbg=0
if n_elements(abs) eq 0 then abs=0 ; [flg] Scale to absorption rather than extinction
if n_elements(hrz) eq 0 then hrz=0 ; [flg] High resolution plotting
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/aca/swnb.nc'
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=3.0
;	y_sz=6.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

fl_idx=0
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
dim_id=ncdf_dimid(nc_id,'rsn')
ncdf_diminq,nc_id,dim_id,foo,rsn_nbr
if dbg then print,'Current number of resolutions = ',auto_sng(rsn_nbr,0)
; netCDF character strings are returned as byte arrays to IDL, convert with string(att_val)
;ncdf_attget,nc_id,/GLOBAL,'aer_sng',aer_sng
ncdf_varget,nc_id,'flx_bb_abs_cld',flx_bb_abs_cld
ncdf_varget,nc_id,'flx_bb_abs_atm',flx_bb_abs_atm
ncdf_varget,nc_id,'flx_bb_abs_ttl',flx_bb_abs_ttl
ncdf_varget,nc_id,'flx_bb_abs_sfc',flx_bb_abs_sfc
ncdf_varget,nc_id,'mpc_CWP',mpc_CWP
ncdf_varget,nc_id,'abs_bb_SAS',abs_bb_SAS
ncdf_varget,nc_id,'rfl_bb_SAS',rfl_bb_SAS
ncdf_varget,nc_id,'trn_bb_atm',trn_bb_SAS
ncdf_close,nc_id

abc=[1,10,100,1000,10000,100000,1000000] ; [# (0.01 um)-1] Bands per 0.01 um
abc=abc*100 ; [# (0.01 um)-1]->[# um-1] Bands per um
if abc_max eq 0.0 then abc_max=max(abc)
if abc_min eq 0.0 then abc_min=min(abc)
;ord=flx_bb_abs_atm
ord=flx_bb_abs_cld
foo=max(ord,ord_max_idx)
if ord_max eq 0.0 then ord_max=foo
foo=min(ord,ord_min_idx)
if ord_min eq 0.0 then ord_min=foo
ord_nm='flx_bb_abs_cld'
;ord_nm='flx_bb_abs_atm'
rng_x=[abc_min,abc_max]
rng_y=[ord_min,ord_max]
ord_max=0.0 ; Reset so value does not effect next plot
ord_min=0.0 ; Reset so value does not effect next plot
if prn then chr_sz=1.5 else chr_sz=1.0
;ttl=string(aer_sng)
;ttl='Atmospheric Absorption, MLS Cloud LWP = '+auto_sng(mpc_CWP(0)*1000,1)+' g m!E-2!N'
ttl='Absorption in MLS Cloud,'+auto_sng(mpc_CWP(0)*1000,1)+' g m!E-2!N'
;x_ttl='!5Resolution !8R!5 (# !7l!5m!E-1!N)'
x_ttl='!5Quadrature Density !8N!5 (# !7l!5m!E-1!N)'
y_ttl='!5Absorption (W m!E-2!N)'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
;mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/swnb2_'+ord_nm+'_cnv.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot cnv_gph_swnb2
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=1,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
yll_idx=base_clr_idx+3
orn_idx=base_clr_idx+4
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(orn_idx)=255
g_curr(orn_idx)=128
b_curr(orn_idx)=0
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

fld_nbr=1
lbl_sng=strarr(fld_nbr)
ln_sty=0.0*intarr(fld_nbr)
clr=intarr(fld_nbr)
clr(0)=clr_blk_idx
;clr(1)=grn_idx
;clr(2)=blu_idx
;clr(3)=yll_idx
;clr(4)=orn_idx
lbl_sng(0)='Absorption'
;lbl_sng(1)='Absorption'
;lbl_sng(2)='Scattering'

oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)
;oplot,abc,ext_cff_mss(abc_min_idx:abc_max_idx),max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)

if info then begin
ln_lgn_x1=0.50
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

;xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5'+fl_nm,color=clr_blk_idx,size=chr_sz,/NORMAL

for fl_idx=0,fld_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end cnv_gph_swnb2()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure cnv_gph_swnb2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure idx_rfr_gph: Plots indices of refraction from various sources
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro idx_rfr_gph, $
	fl_nm=fl_nm, $
	info=info, $
	prn=prn
; Usage:
; idx_rfr_gph,fl_nm=getenv('DATA')+'/aca/idx_rfr_DHE02.nc'
@~/idl/ibp_clr.com
!p.multi=0

if n_elements(dbg) eq 0 then dbg=0
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/aca/idx_rfr_DHE02.nc'
if n_elements(fl_mdl) eq 0 then fl_mdl=getenv('DATA')+'/aca/idx_rfr_DHE02.nc'
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=4.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

rgn_nm=['PaG81','DKS91','SAJ93','Bhr','Bnz','Brb','CpV','Mng','Ogd','SdA']
rgn_sng=['PaG81','DKS91','SAJ93','Bhr','Bnz','Brb','CpV','Mng','Ogd','SdA']
rgn_nbr=n_elements(rgn_nm)
cmp_nm=['cor','mnt','mtx','ncl','prt']
cmp_sng=['PaG81','DKS91','SAJ93','Bhr','Bnz','Brb','CpV','Mng','Ogd','SdA']
cmp_nbr=n_elements(cmp_nm)

if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
dim_id=ncdf_dimid(nc_id,'bnd')
ncdf_diminq,nc_id,dim_id,foo,wvl_nbr
ncdf_varget,nc_id,'bnd',wvl_mcr
if dbg then print,'Current number of wavelengths is '+wvl_nbr
idx_rfr_rl=fltarr(wvl_nbr,rgn_nbr)
idx_rfr_img=fltarr(wvl_nbr,rgn_nbr)
for rgn_idx=0,rgn_nbr-1 do begin
	ncdf_varget,nc_id,'idx_rfr_aeronet_'+rgn_nm[rgn_idx]+'_rl',foo
	idx_rfr_rl(*,rgn_idx)=foo
	ncdf_varget,nc_id,'idx_rfr_aeronet_'+rgn_nm[rgn_idx]+'_img',foo
	idx_rfr_img(*,rgn_idx)=foo
	print,'idx_rfr_img(*,rgn_idx)=',idx_rfr_img(*,rgn_idx)
endfor; end loop over rgn
ncdf_close,nc_id

if dbg then print,'Processing '+fl_mdl
nc_id=ncdf_open(fl_mdl)
dim_id=ncdf_dimid(nc_id,'wvl')
ncdf_diminq,nc_id,dim_id,foo,wvl_nbr_mdl
ncdf_varget,nc_id,'wvl',wvl_mdl
if dbg then print,'Current number of wavelengths is '+wvl_nbr_mdl
idx_rfr_rl=fltarr(wvl_nbr,cmp_nbr)
idx_rfr_img=fltarr(wvl_nbr,cmp_nbr)
for cmp_idx=0,cmp_nbr-1 do begin
	ncdf_varget,nc_id,'idx_rfr_aeronet_'+cmp_nm[cmp_idx]+'_rl',foo
	idx_rfr_rl(*,cmp_idx)=foo
	ncdf_varget,nc_id,'idx_rfr_aeronet_'+cmp_nm[cmp_idx]+'_img',foo
	idx_rfr_img(*,cmp_idx)=foo
	print,'idx_rfr_img(*,cmp_idx)=',idx_rfr_img(*,cmp_idx)
endfor; end loop over cmp
ncdf_close,nc_id

abc=wvl_mcr
ord=idx_rfr_img(*,0)
;xcl_mss_dat,abc,ord,max=1.0e20,min=0.0,nbr=nbr
abc_min=min(abc)
abc_max=max(abc)
rng_x=[0.4,1.0]
ord_min=min(ord)
ord_max=max(ord)
;rng_y=[ord_min,ord_max]
rng_y=[0,0.01]
if prn then chr_sz=2.0 else chr_sz=1
ttl='Refractive Index: Lab vs. AERONET'
print,size(ttl)
x_ttl='!5Wavelength !7k!5 (!7l!5m)'
y_ttl='!5Imaginary Part'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/idx_rfr_gph.eps'
	fl_psn=strpos(fl_nm,'aer_')
	nc_psn=strpos(fl_nm,'.nc')
	if fl_psn gt -1 and nc_psn gt -1 then begin
		fl_nm_out=getenv('DATA')+'/ps/'+strmid(fl_nm,fl_psn,(nc_psn-fl_psn))+'.eps'
	endif; endif		
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot idx_rfr
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=0,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
obs_nbr=3
grn_idx=base_clr_idx+obs_nbr+0
blu_idx=base_clr_idx+obs_nbr+1
red_idx=base_clr_idx+obs_nbr+2
yll_idx=base_clr_idx+obs_nbr+3
orn_idx=base_clr_idx+obs_nbr+4
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(orn_idx)=255
g_curr(orn_idx)=128
b_curr(orn_idx)=0
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

lbl_sng=strarr(rgn_nbr)
ln_sty=indgen(rgn_nbr) mod 4
clr=intarr(rgn_nbr)
for rgn_idx=0,obs_nbr-1 do begin
	clr(rgn_idx)=clr_blk_idx
endfor; end loop over rgn
clr(obs_nbr+0)=grn_idx
clr(obs_nbr+1)=blu_idx
clr(obs_nbr+2)=yll_idx
clr(obs_nbr+3)=orn_idx
clr(obs_nbr+4)=red_idx
for rgn_idx=obs_nbr+5,rgn_nbr-1 do begin
	clr(rgn_idx)=clr(1+(rgn_idx mod 6))
endfor; end loop over rgn

for rgn_idx=0,rgn_nbr-1 do begin
	oplot,abc,idx_rfr_img(*,rgn_idx),max_value=1.0e20,thick=2.0,linestyle=ln_sty(rgn_idx),color=clr(rgn_idx)
endfor; end loop over regions

if info then begin
ln_lgn_x1=0.70
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.85
lgn_dy=0.05
lgn_y=lgn_y_top-indgen(rgn_nbr)*lgn_dy

;xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5'+fl_nm,color=clr_blk_idx,size=chr_sz,/NORMAL

for rgn_idx=0,rgn_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(rgn_idx)+0.013,linestyle=ln_sty(rgn_idx),color=clr(rgn_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(rgn_idx),rgn_sng(rgn_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end idx_rfr_gph()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure idx_rfr_gph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure odx_htg_gph: Graphs heating profiles from swnb simulations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro odx_htg_gph, $
	fl_nm=fl_nm, $
	info=info, $
	prn=prn
; Usage:
; odx_htg_gph,fl_nm=getenv('DATA')+'/tmp/mls_icrccm_clr_cln_O2O2_O2N2.nc'

@~/idl/ibp_clr.com
!p.multi=0

if n_elements(dbg) eq 0 then dbg=0
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(fl_nm) eq 0 then fl_nm=getenv('DATA')+'/tmp/mls_icrccm_clr_cln_O2O2_O2N2.nc'
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=6.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

fl_idx=0
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
dim_id=ncdf_dimid(nc_id,'lev')
ncdf_diminq,nc_id,dim_id,foo,lev_nbr
if dbg then print,'Current number of levels is '+lev_nbr
ncdf_varget,nc_id,'htg_rate_bb',htg_ttl
ncdf_varget,nc_id,'lev',lev
ncdf_varget,nc_id,'tpt',tpt
;ncdf_varget,nc_id,'foo',foo
ncdf_close,nc_id
dns_mdp=lev/(tpt*287.05) ; [kg m-3] Density

fl_nm=getenv('DATA')+'/tmp/mls_icrccm_clr_cln_O2O2_O2N2_frc.nc'
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
ncdf_varget,nc_id,'htg_rate_bb',hgt_frc_O2X_clr
ncdf_close,nc_id

fl_nm=getenv('DATA')+'/tmp/mls_icrccm_clr_cln_CO2_frc.nc'
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
ncdf_varget,nc_id,'htg_rate_bb',htg_frc_CO2
ncdf_close,nc_id

fl_nm=getenv('DATA')+'/tmp/mls_icrccm_cld_cln_O2O2_O2N2_frc.nc'
if dbg then print,'Processing '+fl_nm
nc_id=ncdf_open(fl_nm)
ncdf_varget,nc_id,'htg_rate_bb',hgt_frc_O2X_cld
ncdf_close,nc_id

abc=htg_ttl*86400.0 ; [K s-1] --> [K day-1]
abc_min=min(abc)
abc_max=max(abc)
;rng_x=[abc_min,abc_max]
rng_x=[0.0,4.0]
ord=lev/100.0 ; [Pa] --> [mb]
;ord=dns_mdp ; [kg m-3]
ord_min=min(ord)
ord_max=max(ord)
;rng_y=[ord_min,ord_max]
rng_y=[1000.0,0.0]
if prn then chr_sz=2.0 else chr_sz=1
ttl='MLS Solar Gaseous Heating'
x_ttl='!5Heating rate (!5K day!E-1!N)'
y_ttl='!5Pressure !8p!5 (mb)'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/o2x_htg.eps'
	fl_psn=strpos(fl_nm,'aer_')
	nc_psn=strpos(fl_nm,'.nc')
	if fl_psn gt -1 and nc_psn gt -1 then begin
		fl_nm_out=getenv('DATA')+'/ps/'+strmid(fl_nm,fl_psn,(nc_psn-fl_psn))+'.eps'
	endif; endif		
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot odx_htg_gph
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=0,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
yll_idx=base_clr_idx+3
orn_idx=base_clr_idx+4
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(orn_idx)=255
g_curr(orn_idx)=128
b_curr(orn_idx)=0
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

fld_nbr=3
lbl_sng=strarr(fld_nbr)
ln_sty=intarr(fld_nbr)
ln_sty(0)=0
ln_sty(1)=2
ln_sty(2)=1
clr=intarr(fld_nbr)+clr_blk_idx
lbl_sng(0)='Total'
;lbl_sng(1)='!5O!I2!N!9. !5X!N clear'
;lbl_sng(2)='!5O!I2!N!9. !5X!N cloudy'
lbl_sng(1)='100!9X!5!5O!I2!N!9. !5X!N clear'
lbl_sng(2)='100!9X!5!5O!I2!N!9. !5X!N cloud'
;lbl_sng(2)='100 !9X!5 !5CO!I2!N'

oplot,htg_ttl*86400.0,ord,max_value=1.0e20,thick=2.0,linestyle=ln_sty(0),color=clr(0)
oplot,hgt_frc_O2X_clr*8640000.0,ord,max_value=1.0e20,thick=2.0,linestyle=ln_sty(1),color=clr(1)
oplot,hgt_frc_O2X_cld*8640000.0,ord,max_value=1.0e20,thick=2.0,linestyle=ln_sty(2),color=clr(2)
;oplot,htg_frc_O2*8640000.0,ord,max_value=1.0e20,thick=2.0,linestyle=ln_sty(2),color=clr(2)
if info then begin
ln_lgn_x1=0.53
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.87
lgn_dy=0.07
lgn_y=lgn_y_top-indgen(6)*lgn_dy

;xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!5'+fl_nm,color=clr_blk_idx,size=chr_sz,/NORMAL

for fld_idx=0,fld_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fld_idx)+0.013,linestyle=ln_sty(fld_idx),color=clr(fld_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fld_idx),lbl_sng(fld_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fld
endif; info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end odx_htg_gph()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure odx_htg_gph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro gcm_vy_gph,xpt_nm=xpt_nm,mth=mth,rvr=rvr,ttl=ttl,info=info,prn=prn
; Purpose: Plot zonal average quantities from CCM runs, inverting seasons
; gcm_vy_gph,info=1,prn=0
; gcm_vy_gph,info=1,mth=['01','07'],xpt_nm='dmr04',prn=0
; gcm_vy_gph,info=1,mth=['1202','0608'],xpt_nm='dmr04',prn=0   
if n_elements(info) eq 0 then info=0
if n_elements(xpt_nm) eq 0 then xpt_nm='dmr04'
if n_elements(mth) eq 0 then mth=['1202','0608']
if n_elements(rvr) eq 0 then rvr=0
if n_elements(prn) eq 0 then prn=0
if n_elements(dbg) eq 0 then dbg=0

; Color initialization
@ibp_clr.com
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr

; Printer initialization
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=3.2
	if top then y_sz=y_sz*1.11
	if btm then y_sz=y_sz*1.22
endif; endif prn

fld_stf=[ $
;	['FSATFRC','!5Tracer Atmospheric Absorption','!5W m!E-2!N'] ]
	['NPCO2O2','!5O!I2!N!9. !5O!I2!N Column','!9X!5 10!E42!N molecule!E2!N cm!E-5!N'] ]
fld_nbr=n_elements(fld_stf(0,*))
foo={fld_sct,idx:0,nm:'',sng:'',unit:'',scl:0.0}
fld=replicate({fld_sct},fld_nbr)
fl_nm_fl_outd_sng=''
for fld_idx=0,fld_nbr-1 do begin
	fld(fld_idx).idx=fld_idx
	fld(fld_idx).nm=fld_stf(0,fld_idx)
	fld(fld_idx).sng=fld_stf(1,fld_idx)
	fld(fld_idx).unit=fld_stf(2,fld_idx)
	fld(fld_idx).scl=1.0
	if fld(fld_idx).nm eq 'NPCO2O2' then fld(fld_idx).scl=1.0e-16
	fl_nm_fl_outd_sng=fl_nm_fl_outd_sng+fld(fld_idx).nm
	if fld_idx ne fld_nbr-1 then fl_nm_fl_outd_sng=fl_nm_fl_outd_sng+'_'
endfor; end loop over fld

mth_nbr=n_elements(mth)
fl_nm_out_mth_sng=''
for mth_idx=0,mth_nbr-1 do begin
	fl_nm_out_mth_sng=fl_nm_out_mth_sng+mth(mth_idx)+'_'
endfor; end loop over mth

for fld_idx=0,fld_nbr-1 do begin
for mth_idx=0,mth_nbr-1 do begin
	fl_nm=getenv('DATA')+'/'+xpt_nm+'/'+xpt_nm+'_8589_'+mth(mth_idx)+'_x.nc'
	nc_id=ncdf_open(fl_nm)
	lat_id=ncdf_dimid(nc_id,'lat')
	ncdf_diminq,nc_id,lat_id,dim_foo,lat_nbr
	ncdf_varget,nc_id,'lat',lat
	if mth_idx eq 0 and fld_idx eq 0 then data=fltarr(lat_nbr,mth_nbr,fld_nbr)
	ncdf_varget,nc_id,fld(fld_idx).nm,data_foo
	data_foo=reform(data_foo)
	data(*,mth_idx,fld_idx)=data_foo*fld(fld_idx).scl
	ncdf_close,nc_id
endfor; end loop over mth
endfor; end loop over fld

abc=data(*,0,0)
ord=lat
rng_x=[min(data),max(data)]
;rng_x=[min(abc),max(abc)]
;rng_y=[min(ord),max(ord)]
rng_y=[-90.0,90.0]
chr_sz=1.5
if n_elements(ttl) eq 0 then if fld(0).nm eq 'NPCO2O2' then ttl='!5Summer vs. Winter !5O!I2!N!9. !5O!I2!N Abundance'
if n_elements(ttl) eq 0 then if fld(0).nm eq 'FSATFRC' then ttl='!5Summer vs. Winter Tracer Forcing'
ttl=trc_nm_sbs(ttl,xpt_nm)
x_ttl=trc_nm_sbs(fld(0).sng,xpt_nm)+' ('+fld(0).unit+')'
y_ttl='!5Latitude'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/'+xpt_nm+'_8589_'+fl_nm_out_mth_sng+'x_'+fl_nm_fl_outd_sng+'.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot npc
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=0,ystyle=13,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

ln_sty=intarr(mth_nbr)
clr=intarr(mth_nbr)
for fld_idx=0,fld_nbr-1 do begin
for mth_idx=0,mth_nbr-1 do begin
	ln_sty(mth_idx)=2*mth_idx
	clr(mth_idx)=clr_blk_idx
	abc=data(*,mth_idx,fld_idx)
	if rvr and mth_idx mod 2 eq 1 then abc=reverse(abc)
	oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=ln_sty(mth_idx),color=clr(mth_idx)
endfor; end loop over mth
endfor; end loop over fld

lat_axz_drw,lat_top=!y.crange(1),lat_btm=!y.crange(0),ntt=1,ttl=y_ttl,chr_sz=chr_sz,axz_vrt=1

if info then begin
ln_lgn_x1=0.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.6
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy
for mth_idx=0,mth_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(mth_idx)+0.013,linestyle=ln_sty(mth_idx),color=clr(mth_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(mth_idx),mth2lbl(mth(mth_idx)),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over mth
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end gcm_vy_gph()

pro trn_abs_gph,day=day,ss_alb=ss_alb,info=info,prn=prn
; Purpose: Plot transmission vs. absorption result from SWNB runs
; trn_abs_gph,day=[11],ss_alb=[0.83],info=1,prn=0
if n_elements(ss_alb) eq 0 then ss_alb=[0.83]
if n_elements(info) eq 0 then info=0
if n_elements(day) eq 0 then day=[11]
if n_elements(prn) eq 0 then prn=0

; Some initialization
sym_usr_foo=findgen(17)*(!pi*2/16.)
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=6.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

; NB: Data from CZV98 Fgr. 4
tddr_trn_obs=[0.94039,0.90052,0.89467,0.88897,0.88546,0.88520,0.88538]
tddr_abs_obs=[0.0075227,0.047542,0.051726,0.057981,0.059082,0.057251,0.057514]

fl_nbr=n_elements(ss_alb)
ss_alb_sng=strarr(fl_nbr)
fl_nm=strarr(fl_nbr)
for fl_idx=0,fl_nbr-1 do begin
	foo=fix(100.0*ss_alb(fl_idx) mod 100)
	ss_alb_sng(fl_idx)=string(format='(I3.3)',foo)
	fl_nm(fl_idx)=getenv('DATA')+'/arese/arese_951011_w_'+ss_alb_sng(fl_idx)+'_trn_abs.nc'
endfor; end loop over fl

; Get coordinate sizes to initialize multidimensional variables
nc_id=ncdf_open(fl_nm(0))
aod_id=ncdf_dimid(nc_id,'odxc_aer')
ncdf_diminq,nc_id,aod_id,dim_foo,aod_nbr
tddr_trn_mdl=fltarr(aod_nbr,fl_nbr)
tddr_abs_mdl=fltarr(aod_nbr,fl_nbr)
ncdf_close,nc_id

;fl_nm=[getenv('DATA')+'/arese/951011_arese_trn_abs.nc']
fl_nbr=n_elements(fl_nm)
for fl_idx=0,fl_nbr-1 do begin
	nc_id=ncdf_open(fl_nm(fl_idx))
	ncdf_varget,nc_id,'tddr_trn_mdl',data_foo
	tddr_trn_mdl(*,fl_idx)=data_foo
 	ncdf_varget,nc_id,'tddr_abs_mdl',data_foo
	tddr_abs_mdl(*,fl_idx)=data_foo
	ncdf_close,nc_id
endfor; end loop over fl

abc_1=tddr_trn_obs
ord_1=tddr_abs_obs
xcl_mss_dat,abc_1,ord_1,max=1.0,min=0.0,nbr=nbr_1
rng_x=[min(tddr_trn_mdl),0.95]
rng_y=[0.0,max(tddr_abs_mdl)]
;rng_x=[0.87,0.95]
;rng_y=[0.0,0.08]
chr_sz=1.5
ttl='ARESE 951011 Aircraft TDDR 0.5 !7l!5m'
x_ttl='!5Transmittance !8T!5'
y_ttl='!5Absorptance !8A!5'

if prn then begin
	mrg_top=0.6 ; 2 is default
;	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
;	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/arese_951011_w_'+ss_alb_sng(0)+'_trn_abs.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps
endif; endif prn

; scatterplot transmission vs. absorption
!p.multi=0
usersym,cos(sym_usr_foo),sin(sym_usr_foo),/fill
plot,abc_1,ord_1,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xstyle=0,ystyle=0,thick=2.0,charsize=chr_sz,xmargin=[8.5,2],ymargin=[3.5,2],psym=8

abc_2=tddr_trn_mdl(*,0)
ord_2=tddr_abs_mdl(*,0)
xcl_mss_dat,abc_2,ord_2,max=1.0,min=0.0,nbr=nbr_2
usersym,cos(sym_usr_foo),sin(sym_usr_foo)
if nbr_2 ne 0 then oplot,abc_2,ord_2,psym=8

ln_lgn_x1=0.60
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.7
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

if info then begin
xyouts,txt_lgn_x,lgn_y(0)+lgn_dy,'!7D!8A/!7D!8T!5 Slope:',size=1.0*chr_sz,/NORMAL
; Plot the least squares regressions lines on the figure
x_data=[abc_1]
y_data=[ord_1]
fit_cff=poly_fit(x_data,y_data,1,ord_fit)
;oplot,x_data,ord_fit,thick=3.0,linestyle=0	 
oplot,[min(x_data),max(x_data)],[fit_cff(0)+min(x_data)*fit_cff(1),fit_cff(0)+max(x_data)*fit_cff(1)],thick=3.0,linestyle=0	 
m_sng=auto_sng(fit_cff(1),2)
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(0)+0.013,linestyle=0,thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(0),'Obs = '+m_sng,size=1.0*chr_sz,/NORMAL

x_data=[abc_2]
y_data=[ord_2]
fit_cff=poly_fit(x_data,y_data,1,ord_fit)
;oplot,x_data,ord_fit,thick=3.0,linestyle=2	 
oplot,[min(x_data),max(x_data)],[fit_cff(0)+min(x_data)*fit_cff(1),fit_cff(0)+max(x_data)*fit_cff(1)],thick=2.0,linestyle=2	 
m_sng=auto_sng(fit_cff(1),2)
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(1)+0.013,linestyle=2,thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'Mdl = '+m_sng,size=1.0*chr_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(2),'(!7x!5 = '+auto_sng(ss_alb(0),2)+')',size=1.0*chr_sz,/NORMAL
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end trn_abs_gph()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure alb_ocn
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro alb_ocn, $
	info=info, $
	prn=prn, $ 
	scl=scl, $ 
	slr_cst=slr_cst, $ 
	odxc_Ray=odxc_Ray, $ 
	odxc_aer=odxc_aer, $ 
	ss_alb=ss_alb, $ 
	asm_prm=asm_prm, $ 
	ngl_nbr=ngl_nbr, $
	dbg=dbg
; Purpose: Compute and display ocean surface albedo as specified by CCM SW code 

; alb_ocn,dbg=1,ngl_nbr=91,prn=0
; alb_ocn,odxc_aer=0.27,prn=0

@~/idl/ibp_clr.com
!p.multi=0

if n_elements(ngl_nbr) eq 0 then ngl_nbr=91
if n_elements(odxc_Ray) eq 0 then odxc_Ray=0.15
if n_elements(odxc_aer) eq 0 then odxc_aer=0.14
if n_elements(ss_alb) eq 0 then ss_alb=0.95
if n_elements(asm_prm) eq 0 then asm_prm=0.7
if n_elements(slr_cst) eq 0 then slr_cst=1367.0
if n_elements(scl) eq 0 then scl=0

if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22
if n_elements(dbg) eq 0 then dbg=0

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=6.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

; ngl goes from 0 to 0.5*pi
ngl=0.5*!pi*findgen(ngl_nbr)/(ngl_nbr-1.0)
ngl_dgr=180.0*ngl/!pi
ngl_dlt=0.5*!pi*(fltarr(ngl_nbr)+1.0)/(ngl_nbr-1.0)
cos_ngl=cos(ngl)
alb_ocn_drc=0.026/(0.065+cos_ngl^1.7) + 0.15*(cos_ngl-0.10)*(cos_ngl-0.50)*(cos_ngl-1.0)
alb_ocn_dff=fltarr(ngl_nbr)+0.06

bsf_Ray=0.5 ; Backscattered fraction for Rayleigh scattering
bsf_aer=(1.0-asm_prm)/2.0 ; Backscattered fraction for aerosol
odxc_ttl=odxc_Ray+odxc_aer
odxc_drc=odxc_ttl

if scl then begin
	odxc_Ray=odxc_Ray*(1.0-1.0*0.1)
	odxc_aer=odxc_aer*(1.0-ss_alb*asm_prm*asm_prm)
	odxc_drc=odxc_ttl+(odxc_ttl-(odxc_Ray+odxc_aer))
endif; not scl

;ss_alb_ttl=(1.0*odxc_Ray+ss_alb*odxc_aer)/odxc_ttl
;asm_prm_ttl=(0.0*1.0*odxc_Ray+asm_prm*ss_alb*odxc_aer)/(ss_alb_ttl*odxc_ttl)
;fsf_ttl=(0.1*1.0*odxc_Ray+fsf_aer*ss_alb*odxc_aer)/(ss_alb_ttl*odxc_ttl) ; ???

flx_dwn_drc=cos_ngl*slr_cst*exp(-odxc_drc/cos_ngl)

flx_dff_aer=cos_ngl*slr_cst*(1.0-exp(-odxc_aer/cos_ngl))*ss_alb ; Aerosol scattered diffuse flux
flx_up_dff_aer=flx_dff_aer*bsf_aer
flx_dwn_dff_aer=flx_dff_aer*(1.0-bsf_aer)

flx_dff_Ray=cos_ngl*slr_cst*(1.0-exp(-odxc_Ray/cos_ngl))*1.0 ; Rayleigh scattered diffuse flux
flx_up_dff_Ray=flx_dff_Ray*bsf_Ray
flx_dwn_dff_Ray=flx_dff_Ray*(1.0-bsf_Ray)

flx_dwn_dff=flx_dwn_dff_Ray+flx_dwn_dff_aer
flx_up_dff_ocn=flx_dwn_dff*alb_ocn_dff+flx_dwn_drc*alb_ocn_drc
flx_up_dff=flx_up_dff_aer+flx_up_dff_Ray+flx_up_dff_ocn

bad_idx=where(abs(cos_ngl) lt 1.0e-5,cnt)
if cnt gt 0 then flx_dwn_drc(bad_idx)=0.0
if cnt gt 0 then flx_dwn_dff_aer(bad_idx)=0.0
if cnt gt 0 then flx_dwn_dff_Ray(bad_idx)=0.0
if cnt gt 0 then flx_dwn_dff(bad_idx)=0.0
if cnt gt 0 then flx_up_dff_aer(bad_idx)=0.0
if cnt gt 0 then flx_up_dff_Ray(bad_idx)=0.0
if cnt gt 0 then flx_up_dff_ocn(bad_idx)=0.0
if cnt gt 0 then flx_dwn_dff(bad_idx)=0.0

; Compute TOA energy balance
flx_nrg_bln= $ ; Energy balance
	+cos_ngl*slr_cst $ ; Input flux
	-flx_up_dff $ ; Reflected flux
	-flx_dwn_dff*(1.0-alb_ocn_dff) $ ; Sfc abs dff flx
	-flx_dwn_drc*(1.0-alb_ocn_drc) $ ; Sfc abs drc flx
	-cos_ngl*slr_cst*(1.0-exp(-odxc_aer/cos_ngl))*(1.0-ss_alb) ; Aer abs flx

if dbg then print,'flx_nrg_bln - ',flx_nrg_bln

abc=ngl_dgr
abc_min=min(abc)
abc_max=max(abc)
rng_x=[abc_min,abc_max]
ord=alb_ocn_drc
ord_min=min(ord)
ord_max=max(ord)
rng_y=[ord_min,ord_max]
chr_sz=2.0
ttl='CCM Ocean Surface Albedo'
x_ttl='!5Solar Zenith Angle !7h!5 (degrees)'
y_ttl='!5Ocean Surface Albedo !8A!5'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/ocn_alb.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot alb
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

oplot,abc,alb_ocn_drc,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx
oplot,abc,alb_ocn_dff,max_value=1.0e20,thick=2.0,linestyle=0,color=grn_idx

fl_nbr=4
lbl_sng=strarr(fl_nbr)
ln_sty=0.0*intarr(fl_nbr)
clr=intarr(fl_nbr)
clr(0)=clr_blk_idx
clr(1)=grn_idx
clr(2)=blu_idx
clr(3)=red_idx
lbl_sng(0)='Direct'
lbl_sng(1)='Diffuse'

if info then begin
ln_lgn_x1=0.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

;xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!7s!5!Iaer = '+auto_sng(odxc_aer,2),color=clr_blk_idx,size=chr_sz,/NORMAL

for fl_idx=0,1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; not info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

ord=flx_dwn_dff
ord_min=min(ord)
;ord_max=max(ord)
ord_max=500.0
rng_y=[ord_min,ord_max]
ttl='Approx. CCM Maritime Fluxes'
x_ttl='!5Solar Zenith Angle !7h!5 (degrees)'
y_ttl='!5Upwelling Flux !8F!5 (W m!E-2!N)'

if prn then begin
	fl_nm_out=getenv('DATA')+'/ps/ocn_flx.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot flx_up
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

oplot,abc,flx_dwn_drc,thick=2.0,linestyle=0,color=clr_blk_idx
oplot,abc,flx_up_dff,thick=2.0,linestyle=0,color=grn_idx
oplot,abc,flx_dwn_dff,thick=2.0,linestyle=0,color=blu_idx
oplot,abc,flx_nrg_bln,thick=2.0,linestyle=0,color=red_idx

lbl_sng(0)='flx_dwn_drc'
lbl_sng(1)='flx_up_dff'
lbl_sng(2)='flx_dwn_dff'
lbl_sng(3)='flx_nrg_bln'

if info then begin
xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!7s!5!Iaer!N = '+auto_sng(odxc_aer,2)+', !8g!5!Iaer!N = '+auto_sng(asm_prm,2)+', !7x!5!Iaer!N = '+auto_sng(ss_alb,2),color=clr_blk_idx,size=chr_sz,/NORMAL
for fl_idx=0,fl_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),lbl_sng(fl_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end alb_ocn()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End Figure alb_ocn
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Begin Figure phz_fnc
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro phz_fnc, $
	fov_dgr=fov_dgr, $
	ngl_nbr=ngl_nbr, $
	asm_prm=asm_prm, $
	dbg=dbg
; Computes phase function forward scattering corrections

; phz_fnc,dbg=1,ngl_nbr=361,fov_dgr=2,asm_prm=0.70
if n_elements(fov_dgr) eq 0 then fov_dgr=2.0 ; observational full-width field-of-view in degrees
if n_elements(ngl_nbr) eq 0 then ngl_nbr=181
if n_elements(asm_prm) eq 0 then asm_prm=0.7
if n_elements(dbg) eq 0 then dbg=0

; Convert fov_dgr to halfwidth (using symmetry of phase function) in radians
fov=0.5*!pi*fov_dgr/180.0

; ngl goes from 0 to pi
ngl=!pi*findgen(ngl_nbr)/(ngl_nbr-1.0)
ngl_dgr=180.0*ngl/!pi
ngl_dlt=!pi*(fltarr(ngl_nbr)+1.0)/(ngl_nbr-1.0)
cos_ngl=cos(ngl)
sin_ngl=sin(ngl)
phg=(1.0-asm_prm*asm_prm)/((1.0+asm_prm*asm_prm-2.0*asm_prm*cos_ngl)^1.5)
; Normalize and sum
phg_ttl=total(0.5*phg*sin_ngl*ngl_dlt)

; Compute partial sums 
phg_psm=phg
phg_psm(0)=0.5*phg(0)*sin_ngl(0)*ngl_dlt(0)
for idx=1,ngl_nbr-1 do begin
	phg_psm(idx)=phg_psm(idx-1)+0.5*phg(idx)*sin_ngl(idx)*ngl_dlt(idx)
endfor; end loop over ngl
phg_frc=phg_psm

; Compute forward scattered fraction for given detector

fl_nm_out=getenv('DATA')+'/mie/phz_fnc_HG.txt'
prn_unit=1 
close,prn_unit
openw,prn_unit,fl_nm_out
printf,prn_unit,'This file generated automatically by ~/idl/mie.pro on '+systime(0)
printf,prn_unit,'Columns 1--4 are Henyey Greenstein phase function'
printf,prn_unit,'idx	ngl_dgr	ngl	cos_ngl	phg	phg_psm	phg_frc'
for idx=0,ngl_nbr-1 do begin
printf,prn_unit,string(format='(7(F7.3,1x))',idx,ngl_dgr(idx),ngl(idx),cos_ngl(idx),phg(idx),phg_psm(idx),phg_frc(idx))
endfor; end loop over ngl
close,prn_unit

end; end phz_fnc()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure psd
; Purpose: Plots number, surface area, and volume distributions for 
; arbitrary particle size distributions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro psd_bch, $
	ln_flg=ln_flg, $ ; [flg] Compute distributions as functions of natural logarithm of the size
	nrm_flg=nrm_flg, $ ; [flg] Plot vertical axis of distributions in normalized coords [0.0,1.0]
	scl_ord_flg=scl_ord_flg, $ ; [flg] Scale ordinate to microns and cubic centimeters
	anl_flg=anl_flg, $ ; [flg] Compute distributions interactively rather than reading from disk
	log_flg=log_flg, $ ; [flg] Compute distributions as functions of base 10 logarithm of the size
	xlg_flg=xlg_flg, $ ; [flg] Plot abscissa on logarithmic axis
	ylg_flg=ylg_flg, $ ; [flg] Plot ordinate on logarithmic axis
	stt_flg=stt_flg, $ ; [flg] Plot statistics of diameter, e.g., dmt_nmr, dmt_nwr
	dmt_flg=dmt_flg, $ ; [flg] Plot distributions as function of diameter not radius
	fmt_nsm=fmt_nsm, $ ; [flg] Format for number, surface area, mass column of graphs
	mmr_rsl=mmr_rsl, $ ; [kg kg-1] Mixing ratio
	fl_lbl=fl_lbl, $ ; [sng] Label to append to file
	dns_mdp=dns_mdp, $ ; [kg m-3] Density
	gsd_anl=gsd_anl, $ ; [frc] Geometric standard deviation
	dns_aer=dns_aer, $ ; [kg m-3] Aerosol density
	psd_nbr=psd_nbr, $ ; [nbr] Number of size bins
	sgs_nbr=sgs_nbr, $ ; [nbr] Number of sub-gridscale bins
	sz_mnm=sz_mnm, $ ; [m] Minimum particle size
	sz_mxm=sz_mxm, $ ; [m] Maximum particle size
	dsd_idx_dbg=dsd_idx_dbg, $ ; [idx] Debugging index for droplet size distribution
	sz_idx_dbg=sz_idx_dbg, $ ; [idx] Debugging index for particle size distribution
	dmt_pcp_dbg=dmt_pcp_dbg, $ ; [m] Debugging diameter for droplet size distribution
	foo=foo, $ ; []
	dbg=dbg, $
	info=info, $
	prn=prn

; Usage:
; psd_bch,dbg,log_flg=0
; psd_bch,mmr_rsl=[3.64e-08,1.63e-07,1.28e-07,3.93e-08] ; dstmch41 Annual mean DSTQXX at Lake Chad 
; psd_bch,mmr_rsl=[5.26e-08,2.28e-07,8.77e-08,1.32e-08] ; dstmch41 Annual mean DSTQXX at Lake Chad 
; psd_bch,mmr_rsl=[1.17e-09,1.50e-08,1.21e-08,1.49e-09] ; dstmch41 July mean DSTQXX at Barbados
; psd_bch,mmr_rsl=[1.54e-09,6.68e-09,5.51e-09,1.05e-09],fl_lbl='_dstmch90_199807_Brb',stt_flg=0,prn=0 ; dstmch90 July mean DSTQXX at Barbados
; psd_bch,mmr_rsl=[2.09e-09,9.13e-09,8.18e-09,1.58e-09],fl_lbl='_dstmch90_clm07_Brb',stt_flg=0,prn=0 ; dstmch90 climatological July mean DSTQXX at Barbados
; psd_bch,mmr_rsl=[2.46e-08,1.08e-07,7.95e-08,1.79e-08] ; dstmch41 Annual mean DSTQXX at Adkins' ODP658C site
; psd_bch,anl_flg=1 ; 
; psd_bch,mmr_rsl=[1.17e-09,1.50e-08,1.21e-08,1.49e-09],anl_flg=1 ; 
; psd_bch,mmr_rsl=[1.17e-09,1.50e-08,1.21e-08,1.49e-09],log_flg=0 ; 
; psd_bch,mmr_rsl=[1.17e-09,1.50e-08,1.21e-08,1.49e-09],dmt_flg=0 ; 
; psd_bch,mmr_rsl=[1.17e-09,1.50e-08,1.21e-08,1.49e-09],fmt_nsm=1,prn=0 ; 
; Find relative mass amounts predicted by dust model for a point or a region, time-averaged or singular event, using, e.g.,
; ncks -H -C -F -u -s "%9.2e, " -d lev,1000.0 -d lat,15.0,15.0 -d lon,15.0,15.0 -v DSTQ01,DSTQ02,DSTQ03,DSTQ04 ${DATA}/dstmch46/dstmch46_1998.nc ; dstmch46 Annual mean DSTQXX at Lake Chad
; ncks -H -C -F -u -s "%9.2e, " -d lev,1000.0 -d lat,15.0,15.0 -d lon,305.0,305.0 -v DSTQ01,DSTQ02,DSTQ03,DSTQ04 ${DATA}/dstmch46/dstmch46_199807.nc ; dstmch46 July mean DSTQXX at Barbados
; ncks -H -C -F -u -s "%9.2e, " -d lev,1000.0 -d lat,15.0,15.0 -d lon,305.0,305.0 -v DSTQ01,DSTQ02,DSTQ03,DSTQ04 ${DATA}/dstmch84/dstmch84_clm07.nc ; dstmch84 July mean DSTQXX at Barbados
; ncks -H -C -F -u -s "%9.2e, " -d lev,1000.0 -d lat,15.0,15.0 -d lon,305.0,305.0 -v DSTQ01,DSTQ02,DSTQ03,DSTQ04 ${DATA}/dstmch90/dstmch90_199807.nc ; dstmch90 July mean DSTQXX at Barbados
; ncks -H -C -F -u -s "%9.2e, " -d lev,1000.0 -d lat,22.0 -d lon,345.0 -v DSTQ01,DSTQ02,DSTQ03,DSTQ04 ${DATA}/dstmch46/dstmch46_1998.nc ; dstmch46 Annual mean DSTQXX at Adkins' site
; ncks -H -C -F -u -v DSTQ01,DSTQ02,DSTQ03,DSTQ04 ${DATA}/dstmch39/dstmch39_19940901_19940930_LkC.nc
; ncks -H -C -F -u -v DSTQ01,DSTQ02,DSTQ03,DSTQ04 ${DATA}/dstmch39/dstmch39_19940901_19940930_Brb.nc

if n_elements(info) eq 0 then info=0
if n_elements(dbg) eq 0 then dbg=0
if n_elements(ln_flg) eq 0 then ln_flg=0 ; [flg] Compute distributions as functions of natural logarithm of the size
if n_elements(nrm_flg) eq 0 then nrm_flg=1 ; [flg] Plot vertical axis of distributions in normalized coords [0.0,1.0]
if n_elements(scl_ord_flg) eq 0 then scl_ord_flg=1 ; [flg] Scale ordinate to microns and cubic centimeters
if n_elements(anl_flg) eq 0 then anl_flg=0 ; [flg] Compute distributions interactively rather than reading from disk
if n_elements(log_flg) eq 0 then log_flg=1 ; [flg] Compute distributions as functions of base 10 logarithm of the size
if n_elements(xlg_flg) eq 0 then xlg_flg=1 ; [flg] Plot abscissa on logarithmic axis
if n_elements(ylg_flg) eq 0 then ylg_flg=0 ; [flg] Plot ordinate on logarithmic axis
if n_elements(stt_flg) eq 0 then stt_flg=1 ; [flg] Plot statistics of diameter, e.g., dmt_nmr, dmt_nwr
if n_elements(dmt_flg) eq 0 then dmt_flg=1 ; [flg] Plot distributions as function of diameter not radius
if n_elements(fmt_nsm) eq 0 then fmt_nsm=1 ; [flg] Format for number, surface area, mass column of graphs
if n_elements(fl_lbl) eq 0 then fl_lbl='' ; [sng] Label to append to file
if n_elements(dns_mdp) eq 0 then dns_mdp=1.15 ; [kg m-3] Density
if n_elements(dns_aer) eq 0 then dns_aer=2.5e3 ; [kg m-3] Aerosol density
if n_elements(psd_nbr) eq 0 then psd_nbr=4 ; [nbr] Number of size bins
if n_elements(sgs_nbr) eq 0 then sgs_nbr=200 ; [nbr] Number of sub-gridscale bins
if n_elements(sz_mnm) eq 0 then sz_mnm=0.1e-6 ; [m] Minimum particle size
if n_elements(sz_mxm) eq 0 then sz_mxm=10.0e-6 ; [m] Maximum particle size
if n_elements(dsd_idx_dbg) eq 0 then dsd_idx_dbg=0 ; [idx] Debugging index for droplet size distribution
if n_elements(sz_idx_dbg) eq 0 then sz_idx_dbg=0 ; [idx] Debugging index for particle size distribution
if n_elements(dmt_pcp_dbg) eq 0 then dmt_pcp_dbg=1.0e-3 ; [m] Debugging diameter for droplet size distribution
if n_elements(foo) eq 0 then foo=0 ; []
if n_elements(prn) eq 0 then prn=0 ; [flg] Print to postscript file

; Parameters that depend on command line input
if n_elements(gsd_anl) eq 0 then gsd_anl=replicate(2.0,psd_nbr) ; [frc] Geometric standard deviation

; Vet input
if ln_flg and log_flg then ln_flg=0

psd, $
	ln_flg=ln_flg, $ ; [flg] Compute distributions as functions of natural logarithm of the size
	scl_ord_flg=scl_ord_flg, $ ; [flg] Scale ordinate to microns and cubic centimeters
	anl_flg=anl_flg, $ ; [flg] Compute distributions interactively rather than reading from disk
	log_flg=log_flg, $ ; [flg] Compute distributions as functions of base 10 logarithm of the size
	xlg_flg=xlg_flg, $ ; [flg] Plot abscissa on logarithmic axis
	ylg_flg=ylg_flg, $ ; [flg] Plot ordinate on logarithmic axis
	stt_flg=stt_flg, $ ; [flg] Plot statistics of diameter, e.g., dmt_nmr, dmt_nwr
	dmt_flg=dmt_flg, $ ; [flg] Plot distributions as function of diameter not radius
	fmt_nsm=fmt_nsm, $ ; [flg] Format for number, surface area, mass column of graphs
	mmr_rsl=mmr_rsl, $ ; [kg kg-1] Mixing ratio
	fl_lbl=fl_lbl, $ ; [sng] Label to append to file
	dns_mdp=dns_mdp, $ ; [kg m-3] Density
	gsd_anl=gsd_anl, $ ; [frc] Geometric standard deviation
	dns_aer=dns_aer, $ ; [kg m-3] Aerosol density
	psd_nbr=psd_nbr, $ ; [nbr] Number of size bins
	sgs_nbr=sgs_nbr, $ ; [nbr] Number of sub-gridscale bins
	sz_mnm=sz_mnm, $ ; [m] Minimum particle size
	sz_mxm=sz_mxm, $ ; [m] Maximum particle size
	dsd_idx_dbg=dsd_idx_dbg, $ ; [idx] Debugging index for droplet size distribution
	sz_idx_dbg=sz_idx_dbg, $ ; [idx] Debugging index for particle size distribution
	dmt_pcp_dbg=dmt_pcp_dbg, $ ; [m] Debugging diameter for droplet size distribution
	foo=foo, $ ; []
	info=info, $
	dbg=dbg, $
	prn=prn

end; end psd_bch()

pro psd, $
	ln_flg=ln_flg, $ ; [flg] Compute distributions as functions of natural logarithm of the size
	nrm_flg=nrm_flg, $ ; [flg] Plot vertical axis of distributions in normalized coords [0.0,1.0]
	scl_ord_flg=scl_ord_flg, $ ; [flg] Scale ordinate to microns and cubic centimeters
	anl_flg=anl_flg, $ ; [flg] Compute distributions interactively rather than reading from disk
	log_flg=log_flg, $ ; [flg] Compute distributions as functions of base 10 logarithm of the size
	xlg_flg=xlg_flg, $ ; [flg] Plot abscissa on logarithmic axis
	ylg_flg=ylg_flg, $ ; [flg] Plot ordinate on logarithmic axis
	stt_flg=stt_flg, $ ; [flg] Plot statistics of diameter, e.g., dmt_nmr, dmt_nwr
	dmt_flg=dmt_flg, $ ; [flg] Plot distributions as function of diameter not radius
	fmt_nsm=fmt_nsm, $ ; [flg] Format for number, surface area, mass column of graphs
	mmr_rsl=mmr_rsl, $ ; [kg kg-1] Mixing ratio
	fl_lbl=fl_lbl, $ ; [sng] Label to append to file
	dns_mdp=dns_mdp, $ ; [kg m-3] Density
	gsd_anl=gsd_anl, $ ; [frc] Geometric standard deviation
	dns_aer=dns_aer, $ ; [kg m-3] Aerosol density
	psd_nbr=psd_nbr, $ ; [nbr] Number of size bins
	sgs_nbr=sgs_nbr, $ ; [nbr] Number of sub-gridscale bins
	sz_mnm=sz_mnm, $ ; [m] Minimum particle size
	sz_mxm=sz_mxm, $ ; [m] Maximum particle size
	dsd_idx_dbg=dsd_idx_dbg, $ ; [idx] Debugging index for droplet size distribution
	sz_idx_dbg=sz_idx_dbg, $ ; [idx] Debugging index for particle size distribution
	dmt_pcp_dbg=dmt_pcp_dbg, $ ; [m] Debugging diameter for droplet size distribution
	foo=foo, $ ; []
	dbg=dbg, $
	info=info, $
	prn=prn

@~/idl/ibp_clr.com
!p.multi=0

if n_elements(dbg) eq 0 then dbg=0
if n_elements(prn) eq 0 then prn=0
if n_elements(mss_val) eq 0 then mss_val=1.0e20
if n_elements(info) eq 0 then info=1
if n_elements(ln_flg) eq 0 then ln_flg=0 ; [flg] Compute distributions as functions of natural logarithm of the size
if n_elements(nrm_flg) eq 0 then nrm_flg=1 ; [flg] Plot vertical axis of distributions in normalized coords [0.0,1.0]
if n_elements(scl_ord_flg) eq 0 then scl_ord_flg=1 ; [flg] Scale ordinate to microns and cubic centimeters
if n_elements(anl_flg) eq 0 then anl_flg=0 ; [flg] Compute distributions interactively rather than reading from disk
if n_elements(log_flg) eq 0 then log_flg=1 ; [flg] Compute distributions as functions of base 10 logarithm of the size
if n_elements(xlg_flg) eq 0 then xlg_flg=1 ; [flg] Plot abscissa on logarithmic axis
if n_elements(ylg_flg) eq 0 then ylg_flg=0 ; [flg] Plot ordinate on logarithmic axis
if n_elements(stt_flg) eq 0 then stt_flg=1 ; [flg] Plot statistics of diameter, e.g., dmt_nmr, dmt_nwr
if n_elements(dmt_flg) eq 0 then dmt_flg=1 ; [flg] Plot distributions as function of diameter not radius
if n_elements(fmt_nsm) eq 0 then fmt_nsm=1 ; [flg] Format for number, surface area, mass column of graphs
if n_elements(fl_lbl) eq 0 then fl_lbl='' ; [sng] Label to append to file
if n_elements(dns_mdp) eq 0 then dns_mdp=1.15 ; [kg m-3] Density
if n_elements(dns_aer) eq 0 then dns_aer=2.5e3 ; [kg m-3] Aerosol density
if n_elements(psd_nbr) eq 0 then psd_nbr=4 ; [nbr] Number of size bins
if n_elements(sgs_nbr) eq 0 then sgs_nbr=200 ; [nbr] Number of sub-gridscale bins
if n_elements(sz_mnm) eq 0 then sz_mnm=0.1e-6 ; [m] Minimum particle size
if n_elements(sz_mxm) eq 0 then sz_mxm=10.0e-6 ; [m] Maximum particle size
if n_elements(dsd_idx_dbg) eq 0 then dsd_idx_dbg=0 ; [idx] Debugging index for droplet size distribution
if n_elements(sz_idx_dbg) eq 0 then sz_idx_dbg=0 ; [idx] Debugging index for particle size distribution
if n_elements(dmt_pcp_dbg) eq 0 then dmt_pcp_dbg=1.0e-3 ; [m] Debugging diameter for droplet size distribution
if n_elements(foo) eq 0 then foo=0 ; []
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22

; Parameters that depend on command line input
if n_elements(gsd_anl) eq 0 then gsd_anl=replicate(2.0,psd_nbr) ; [frc] Geometric standard deviation

; Vet input
if ln_flg and log_flg then ln_flg=0
if xlg_flg then x_sty=1 else x_sty=0 ; [idx] Plot box exactly circumscribes abscissa
if ylg_flg then y_sty=1 else y_sty=0 ; [idx] Plot box exactly circumscribes ordinate

; Initialize constants
ln10=alog(10.0)

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=3.0
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

;fl_nm=[getenv('DATA')+'/mie/mie.nc']
;fl_nm=[getenv('DATA')+'/aca/aer_dust_like.nc']
fl_nm=[ $
	getenv('DATA')+'/dst/mie/aer_saharan_dust_01_CCM_SW.nc', $
	getenv('DATA')+'/dst/mie/aer_saharan_dust_02_CCM_SW.nc', $
	getenv('DATA')+'/dst/mie/aer_saharan_dust_03_CCM_SW.nc', $
	getenv('DATA')+'/dst/mie/aer_saharan_dust_04_CCM_SW.nc' $
	] ; end fl_nm[]
fl_nbr=n_elements(fl_nm)
if not anl_flg then psd_nbr=fl_nbr
if n_elements(mmr_rsl) ne 0 and n_elements(mmr_rsl) ne psd_nbr then print,'psd: n_elements(mmr_rsl) = ',auto_sng(n_elements(mmr_rsl),0),' but psd_nbr = ',auto_sng(psd_nbr,0)

; Declare and initialize arrays depending on number of size distributions (bins)
psd_sz=intarr(psd_nbr) ; [nbr] Number of sub-gridscale particle size bins
dsd_sz=intarr(psd_nbr) ; [nbr] Number of sub-gridscale raindrop size bins
sz_idx_srt=intarr(psd_nbr) ; [idx] Starting index of sub-gridscale bins for each particle size distribution
sz_idx_end=intarr(psd_nbr) ; [idx] Ending index of sub-gridscale bins for each particle size distribution
dsd_idx_srt=intarr(psd_nbr) ; [idx] Starting index of sub-gridscale bins for each raindrop size distribution
dsd_idx_end=intarr(psd_nbr) ; [idx] Ending index of sub-gridscale bins for each raindrop size distribution

; Initialize array properly and determine sz_nbr
if anl_flg then begin ; 
; Code to generate a series of size distributions ab initio

dmt_grd=fltarr(psd_nbr+1); [m] Particle diameter grid
if n_elements(gsd_anl) ne psd_nbr then gsd_anl=replicate(gsd_anl,psd_nbr) ; [frc] Geometric standard deviation
if n_elements(dns_aer) ne psd_nbr then dns_aer=replicate(dns_aer,psd_nbr) ; [kg m-3] Aerosol density
if psd_nbr eq 4 then begin
; These parameters should match those specified in psd.pl and dstszdstbd.F
; Standard dust model
	gsd_anl=[ 2.0    ,  2.0   ,  2.0   ,  2.0   ] ; [frc] Geometric standard deviation SBG98 p. 75 Table 1
;	gsd_anl=[ 2.19987    ,  2.19987   ,  2.19987   ,  2.19987   ] ; [frc] Geometric standard deviation PaG77
	dmt_grd=[ 0.1e-6 ,  1.0e-6,  2.5e-6,  5.0e-6, 10.0e-6 ] ; [m] Particle diameter grid
	dns_aer=[ 2.5e+3 ,  2.5e+3, 2.5e+3 ,  2.5e+3 ] ; [kg m-3] Aerosol density
	psd_sz_usr=[200,25,25,25]
endif else begin ; endif psd_nbr eq 4
	szgrd_mnm=0.1e-6; [m] Minimum particle diameter
	szgrd_mxm=10.0e-6; [m] Maximum particle diameter
	szgrd_mk,szgrd_mnm,szgrd_mxm,psd_nbr,grd_sng='logarithmic',sz_grd=dmt_grd
endelse ; endelse psd_nbr ne 4
if dbg then print,'psd_nbr = ',auto_sng(psd_nbr,0)
if dbg then print,'gsd_anl = ',gsd_anl
if dbg then print,'dns_aer = ',dns_aer
if dbg then print,'dmt_grd = ',dmt_grd

dmt_nma=fltarr(psd_nbr) ; [m] Number median diameter analytic
dmt_vma=fltarr(psd_nbr) ; [m] Volume median diameter analytic
for psd_idx=0,psd_nbr-1 do begin
	if n_elements(psd_sz_usr) ne 0 then psd_sz(psd_idx)=psd_sz_usr(psd_idx) else psd_sz(psd_idx)=sgs_nbr ; [nbr] Number of sub-gridscale bins
	if psd_idx eq 0 then sz_idx_srt(psd_idx)=0 else sz_idx_srt(psd_idx)=sz_idx_end(psd_idx-1)+1
	if psd_idx eq 0 then sz_idx_end(psd_idx)=psd_sz(psd_idx)-1 else sz_idx_end(psd_idx)=sz_idx_end(psd_idx-1)+psd_sz(psd_idx)
	if n_elements(dsd_sz_usr) ne 0 then dsd_sz(psd_idx)=dsd_sz_usr(psd_idx) else dsd_sz(psd_idx)=sgs_nbr ; [nbr] Number of sub-gridscale bins
	if psd_idx eq 0 then dsd_idx_srt(psd_idx)=0 else dsd_idx_srt(psd_idx)=dsd_idx_end(psd_idx-1)+1
	if psd_idx eq 0 then dsd_idx_end(psd_idx)=dsd_sz(psd_idx)-1 else dsd_idx_end(psd_idx)=dsd_idx_end(psd_idx-1)+dsd_sz(psd_idx)

	; Set a fundamental statistic equal to the bin center
	; Choosing "right" statistic is key to success of method
	; Number median diameter is a simple, but not necessarily good, choice, used for awhile in dust version 2
	; dmt_nma(psd_idx)=0.5*(dmt_grd(psd_idx)+dmt_grd(psd_idx+1)) ; [m] Number median diameter analytic
endfor; end loop over psd
; Dust version 2 prescribes mass median diameter
dmt_vma=replicate(2.524e-6,psd_nbr) ; [m] Mass median diameter analytic She84 p. 75 Table 1
gsd_anl=replicate(2.0,psd_nbr) ; [frc] Geometric standard deviation SBG98 p. 75 Table 1
; dmt_vma=replicate(4.82e-6,psd_nbr) ; [m] Mass median diameter analytic BSM96 p. 73 Table 2
; gsd_anl=replicate(1.9,psd_nbr) ; [frc] Geometric standard deviation BSM96 p. 73 Table 2
; dmt_vma=replicate(5.6e-6,psd_nbr) ; [m] Mass median diameter analytic PaG77 p. 2080 Table 1
; gsd_anl=replicate(2.2,psd_nbr) ; [frc] Geometric standard deviation PaG77 p. 2080 Table 1
dmt_nma=dmt_vma*exp(-3.0*alog(gsd_anl)*alog(gsd_anl)) ; [m] Number median diameter analytic

endif else begin ; endif anl_flg
; Code to read a series of input files generated by the mie() program

; Get coordinate sizes from all files
for fl_idx=0,fl_nbr-1 do begin
	print,'Ingesting ',fl_nm(fl_idx)
	nc_id=ncdf_open(fl_nm(fl_idx))
	dim_id=ncdf_dimid(nc_id,'sz')
	ncdf_diminq,nc_id,dim_id,foo,sgs_nbr
	psd_sz(fl_idx)=sgs_nbr
	dim_id=ncdf_dimid(nc_id,'dsd_sz')
	ncdf_diminq,nc_id,dim_id,foo,dsd_nbr
	dsd_sz(fl_idx)=dsd_nbr
	ncdf_close,nc_id
	if fl_idx eq 0 then sz_idx_srt(fl_idx)=0 else sz_idx_srt(fl_idx)=sz_idx_end(fl_idx-1)+1
	if fl_idx eq 0 then sz_idx_end(fl_idx)=sgs_nbr-1 else sz_idx_end(fl_idx)=sz_idx_end(fl_idx-1)+sgs_nbr
	if fl_idx eq 0 then dsd_idx_srt(fl_idx)=0 else dsd_idx_srt(fl_idx)=dsd_idx_end(fl_idx-1)+1
	if fl_idx eq 0 then dsd_idx_end(fl_idx)=dsd_nbr-1 else dsd_idx_end(fl_idx)=dsd_idx_end(fl_idx-1)+dsd_nbr
endfor; end loop over fl

endelse ; endif not anl_flg

; Total number of bins is number of size distributions times number of sub-gridscale bins per size distribution
sz_nbr=total(psd_sz)
if dbg then print,'sz_nbr = ',auto_sng(sz_nbr,0)
if dbg then print,'psd: psd_sz = ',psd_sz
if dbg then print,'sz_idx_srt = ',sz_idx_srt
if dbg then print,'sz_idx_end = ',sz_idx_end

dsd_nbr=total(dsd_sz)
if dbg then print,'dsd_nbr = ',auto_sng(dsd_nbr,0)
if dbg then print,'dsd: dsd_sz = ',dsd_sz
if dbg then print,'dsd_idx_srt = ',dsd_idx_srt
if dbg then print,'dsd_idx_end = ',dsd_idx_end

; Initialize size-dependent variables
sfc=fltarr(sz_nbr)
vlm=fltarr(sz_nbr)
mss=fltarr(sz_nbr)
dst=fltarr(sz_nbr)
cnc=fltarr(sz_nbr)
sfc_spc_rsl=fltarr(sz_nbr)
vlm_spc_rsl=fltarr(sz_nbr)
cnc_spc_rsl=fltarr(sz_nbr)
stk_nbr=fltarr(sz_nbr) ; [frc] Stokes number
ryn_nbr_grv=fltarr(sz_nbr) ; [frc] Reynolds number at terminal velocity
shm_nbr=fltarr(sz_nbr) ; [frc] Schmidt number
scv_cff=fltarr(sz_nbr) ; [s-1] Scavenging coefficient
vlc_grv=fltarr(sz_nbr) ; [m s-1] Gravitational settling speed
rds_ctr=fltarr(sz_nbr) ; [m] Radius
dmt_ctr=fltarr(sz_nbr) ; [m] Diameter
dmt_pcp=fltarr(dsd_nbr) ; [m] Raindrop diameter
rds_pcp=fltarr(dsd_nbr) ; [m] Raindrop radius
stk_nbr_crt=fltarr(dsd_nbr) ; [frc] Critical Stokes number
ryn_nbr_pcp=fltarr(dsd_nbr) ; [frc] Reynolds number at terminal velocity of raindrop
; Multi-dimensional arrays
cll_fsh=fltarr(sz_nbr) ; [frc] Collision efficiency
cll_fsh_brn_dff=fltarr(sz_nbr) ; [frc] Collision efficiency of Brownian diffusion
cll_fsh_ntc=fltarr(sz_nbr) ; [frc] Collision efficiency of interception
cll_fsh_mpc=fltarr(sz_nbr) ; [frc] Collision efficiency of impaction
stk_nbr_rlt=fltarr(sz_nbr) ; [frc] Stokes number of relative flow

; Initialize file-dependent variables
flx_vlm_pcp_rsl=fltarr(psd_nbr)
dmt_pcp_nma=fltarr(psd_nbr)
vlc_grv_nwr=fltarr(psd_nbr)
vlc_grv_vwr=fltarr(psd_nbr)
cnc_nbr_anl=fltarr(psd_nbr)
sfc_anl=fltarr(psd_nbr)
vlm_anl=fltarr(psd_nbr)
mss_anl=fltarr(psd_nbr)
cnc_nbr_rsl=fltarr(psd_nbr)
mss_rsl=fltarr(psd_nbr)
rds_nwr=fltarr(psd_nbr)
rds_swr=fltarr(psd_nbr)
rds_vwr=fltarr(psd_nbr)
rds_nmr=fltarr(psd_nbr)
rds_smr=fltarr(psd_nbr)
rds_vmr=fltarr(psd_nbr)
dmt_nwr=fltarr(psd_nbr)
dmt_swr=fltarr(psd_nbr)
dmt_vwr=fltarr(psd_nbr)
dmt_nmr=fltarr(psd_nbr)
dmt_min_min=fltarr(psd_nbr)
dmt_max_max=fltarr(psd_nbr)

if anl_flg then begin ; 
; Code to generate a series of size distributions ab initio

dmt_dlt=fltarr(sz_nbr) ; [m] Diameter bin size
for psd_idx=0,psd_nbr-1 do begin
	sz_ctr=fltarr(psd_sz(psd_idx)) ; [m] Diameter
	sz_dlt=fltarr(psd_sz(psd_idx)) ; [m] Diameter bin size
	szgrd_mk,dmt_grd(psd_idx),dmt_grd(psd_idx+1),psd_sz(psd_idx),grd_sng='logarithmic',sz_dlt=sz_dlt,sz_ctr=sz_ctr
	dmt_ctr(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))=sz_ctr ; [m] Diameter
	dmt_dlt(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))=sz_dlt ; [m] Diameter bin size
	dst(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))=lgn_evl(dmt_ctr(sz_idx_srt(psd_idx):sz_idx_end(psd_idx)),dmt_nma(psd_idx),gsd_anl(psd_idx)) ; [# m-3 m-1 (diameter)]
endfor; end loop over psd
; We have just computed distribution in terms of per unit diameter, as is done in aerosol model
; However, mie() processed files against which we shall compare contain all distributions in terms of per unit radius
; Thus we convert distribution to per unit radius here by multiplying by 2.0
; After conversion, dst is in same units as variable dst in mie() input files
; Note that dst may well be converted back to per unit diameter for plotting purposes by scl_abc later in the procedure
dst=2.0*dst ; [# m-3 m-1 (per unit diameter)] --> [# m-3 m-1 (per unit radius)]

; Derive radial quantities used in rest of routine
rds_nma=dmt_nma/2.0 ; [m] Number median radius analytic
rds_ctr=dmt_ctr/2.0 ; [m] Radius
rds_dlt=dmt_dlt/2.0 ; [m] Radius bin size

; Define concentration and moments of diameter
cnc=dst*rds_dlt ; [# m-3 m-1 (per unit radius)]
xsa=!pi*0.25*dmt_ctr^2 ; [m2] Cross-sectional area
sfc=!pi*dmt_ctr^2 ; [m2] Surface area
vlm=!pi*(1.0/6.0)*dmt_ctr^3 ; [m3] Volume

; Compute resolved quantities
for psd_idx=0,psd_nbr-1 do begin
	cnc_nbr_rsl(psd_idx)=total(cnc(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))) ; [# m-3] Number concentration resolved
;	sfc_rsl(psd_idx)=total(sfc(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))*cnc(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))) ; [m2 m-3] Surface area resolved
;	vlm_rsl(psd_idx)=total(vlm(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))*cnc(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))) ; [m3 m-3] Volume resolved

	mss(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))=vlm(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))*dns_aer(psd_idx) ; [kg] Mass
	mss_rsl(psd_idx)=total(mss(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))*cnc(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))) ; [kg m-3] Mass concentration resolved

	cnc_nbr_anl(psd_idx)=(1.0/2.0)* $ ; 1.0 represents N_0, the total number concentration in all bins
	(errorf(alog(dmt_grd(psd_idx+1)/dmt_nma(psd_idx))/(sqrt(2.0)*alog(gsd_anl(psd_idx))))- $
	 errorf(alog(dmt_grd(psd_idx  )/dmt_nma(psd_idx))/(sqrt(2.0)*alog(gsd_anl(psd_idx)))))
 ; [# m-3] Number concentration analytic
	sfc_anl(psd_idx)=!pi*cnc_nbr_anl(psd_idx)*dmt_nma(psd_idx)*dmt_nma(psd_idx)*exp(2.0*alog(gsd_anl(psd_idx))*alog(gsd_anl(psd_idx))) ; [m2 m-3] Surface area concentration analytic
	vlm_anl(psd_idx)=!pi*(1.0/6.0)*cnc_nbr_anl(psd_idx)*dmt_nma(psd_idx)*dmt_nma(psd_idx)*dmt_nma(psd_idx)*exp(4.5*alog(gsd_anl(psd_idx))*alog(gsd_anl(psd_idx))); [m3 m-3] Volume concentration analytic
	mss_anl(psd_idx)=vlm_anl(psd_idx)*dns_aer(psd_idx); [kg m-3] Mass concentration analytic

endfor; end loop over psd
cnc_nbr_anl_ttl=total(cnc_nbr_anl)
sfc_anl_ttl=total(sfc_anl)
vlm_anl_ttl=total(vlm_anl)
mss_anl_ttl=total(mss_anl)

endif else begin ; endif anl_flg

; Read a series of psd_nbr input files generated by the mie() program
for psd_idx=0,psd_nbr-1 do begin
	if dbg then print,'Processing '+fl_nm(psd_idx)
	nc_id=ncdf_open(fl_nm(psd_idx))
	ncdf_varget_1D,nc_id,'sfc',sfc,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'vlm',vlm,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'mss',mss,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'dst',dst,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'cnc',cnc,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'stk_nbr',stk_nbr,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'ryn_nbr_grv',ryn_nbr_grv,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'shm_nbr',shm_nbr,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'scv_cff',scv_cff,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'vlc_grv',vlc_grv,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'rds_ctr',rds_ctr,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'dmt_ctr',dmt_ctr,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'cnc_nbr_rsl',cnc_nbr_rsl,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'cnc_nbr_anl',cnc_nbr_anl,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'mss_rsl',mss_rsl,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'mss_anl',mss_anl,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'flx_vlm_pcp_rsl',flx_vlm_pcp_rsl,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'dmt_pcp_nma',dmt_pcp_nma,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'vlc_grv_nwr',vlc_grv_nwr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'vlc_grv_vwr',vlc_grv_vwr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'rds_nwr',rds_nwr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'rds_swr',rds_swr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'rds_vwr',rds_vwr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'rds_nmr',rds_nmr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'rds_smr',rds_smr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'rds_vmr',rds_vmr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'dmt_nwr',dmt_nwr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'dmt_swr',dmt_swr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'dmt_vwr',dmt_vwr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'dmt_nmr',dmt_nmr,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'cnc_spc_rsl',cnc_spc_rsl,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'sfc_spc_rsl',sfc_spc_rsl,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'vlm_spc_rsl',vlm_spc_rsl,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'dmt_min_min',dmt_min_min,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'dmt_max_max',dmt_max_max,[psd_idx,psd_idx]
	ncdf_varget_1D,nc_id,'dmt_pcp',dmt_pcp,[dsd_idx_srt(psd_idx),dsd_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'rds_pcp',rds_pcp,[dsd_idx_srt(psd_idx),dsd_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'stk_nbr_crt',stk_nbr_crt,[dsd_idx_srt(psd_idx),dsd_idx_end(psd_idx)]
	ncdf_varget_1D,nc_id,'ryn_nbr_pcp',ryn_nbr_pcp,[dsd_idx_srt(psd_idx),dsd_idx_end(psd_idx)]
;	Multi-dimensional arrays
	if dmt_pcp_dbg ne 0.0 then begin
;	User-specified particular value of dmt_pcp to examine
		foo=min(abs(dmt_pcp[dsd_idx_srt(psd_idx):dsd_idx_end(psd_idx)]-dmt_pcp_dbg),dsd_idx_dbg)
		if dbg gt 0 then print,'dmt_pcp_dbg = ', dmt_pcp_dbg,', dsd_idx_dbg = ',dsd_idx_dbg
	endif; endif dmt_pcp_dbg
	ncdf_varget_1D,nc_id,'cll_fsh',cll_fsh,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)],dmn0_idx=dsd_idx_dbg
	ncdf_varget_1D,nc_id,'cll_fsh_brn_dff',cll_fsh_brn_dff,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)],dmn0_idx=dsd_idx_dbg
	ncdf_varget_1D,nc_id,'cll_fsh_ntc',cll_fsh_ntc,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)],dmn0_idx=dsd_idx_dbg
	ncdf_varget_1D,nc_id,'cll_fsh_mpc',cll_fsh_mpc,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)],dmn0_idx=dsd_idx_dbg
	ncdf_varget_1D,nc_id,'stk_nbr_rlt',stk_nbr_rlt,[sz_idx_srt(psd_idx),sz_idx_end(psd_idx)],dmn0_idx=dsd_idx_dbg
	ncdf_close,nc_id
endfor; end loop over fl
dmt_grd=[dmt_min_min,dmt_max_max(psd_nbr-1)]; [m] Particle diameter grid

endelse ; endif not anl_flg

; rds_grd is only used to place boundaries on plots
rds_grd=dmt_grd/2.0 ; [m] Radius grid

; Prepare scales required to convert from input distribution to plotted ordinate
; scl_fnc is used to convert from linear scale to logarithmic scale
; As discussed in SeP97 p. 418, distributions are subtle 
; It is most intuitive to plot distributions so area under curve is physically meaningful (preserves, e.g., number concentration)
; Plotting linear diameters vs linear distributions does this 
; Plotting logarithmic diameters vs functions of the logarithmic diameter also does this
; Since horizontal axis is usually logarithmic (linear in log_10 D_p), distribution on vertical axis
; should be what SeP97 terms n^0_N, n^0_S, n^0_V, n^0_M so that area under curve preserves meaningful quantity
; Relations between distributions as functions of linear size and logarithm of size are presented on SeP97 p. 418
; scl_fnc being non-unity implements these relations
; Regarding units:
; As discussed in SeP97, dimensions of n^0_N, n^0_S, n^0_V, n^0_M do not end in per unit particle length
; This can be seen from dimensional analysis 
; For example total number concentration does NOT equal \int_0_infty n^0_N d (log_10 D_p)
; This is because taking logarithm of a quantity with physical dimensions is meaningless
; As SeP97 point out, log_10 D_p stands for log_10 (D_p / D_rfr) where D_rfr = 1 m (or 1 um ...)
; Horizontal scale may be labeled, e.g., [10^-1 um..10 um] but one horizontal inch of graph is a constant in log_10 D_p
; For example, if horizontal scale represents 0.1 < D_p < 10.0 um, then -1 < log_10 D_p < 1 is the actual scale no matter what labels say
; Obviously then horizontal scale is physically dimensionless (since it represents log of a quantity)
; Thus when log_flg or ln_flg is set, we follow SeP97 convention and drop m-1 or um-1 suffix from distribution units on y-axis
; ylg_flg determines whether ordinate is plotted on logarithmic axis (regardless of whether ordinate if function of D_p or log D_p)
; Usually it is best to plot distributions so that curve area preserves the integral of the distribution
; In this case usually plot log D_p (log axis) against n^0_X (linear axis) (i.e., log_flg=1, ylg_flg=0)
; If distributions are very wide, plotting D_p (linear) axis against n_X (linear axis) also preserves area (i.e., log_flg=0, ylg_flg=0)
; If distributions are narrow and span many orders of magnitude, then it may be impossible to preserve area and see details
; In this case plot log D_p (log axis) against n^0_X or n_X (logarithmic axis) and area does NOT preserve integral (i.e., log_flg=0 or 1, ylg_flg=1) 
scl_fnc=1.0
if log_flg then scl_fnc=ln10*rds_ctr ; Compute distributions as functions of base 10 logarithm of size SeP97 p. 418
if ln_flg then scl_fnc=rds_ctr ; Compute distributions as functions of natural logarithm of size SeP97 p. 418

; scl_abc is purely a physical dimension scaling to present the independent variable (abscissa) in traditional units of (microns) instead of SI units (meters)
; scl_abc is applied to both abscissa and ordinate so that when abscissa is in microns, ordinate is in per micron
; If scl_abc = 1.0 then abscissa is plotted in default units (meters) and ordinate is in per meter (SI units)
; If scl_abc = 1.0e6 then abscissa is plotted in microns and ordinate is in per micron (traditional units)
; This consistency is required so that area under curve equals plotted distribution 
scl_abc=1.0e6 ; [m] --> [um] Convert from per meter of particle size to per micron of particle size
if dmt_flg then scl_abc=scl_abc*2.0 ; [frc] Convert abscissa from radius to diameter and ordinate from per unit radius to per unit diameter

; scl_ord scales ordinate from SI units to traditional units in microns and cubic centimeters
; For example, SI units for surface area distribution are m2 m-3 m-1 = m-2
; Setting scl_ord_flg converts surface area distribution to um2 cm-3 m-1
scl_ord_nbr_dst=1.0; [# m-3 m-1] --> [# cm-3 m-1] Ordinate scale factor for number distribution
scl_ord_sfc_dst=1.0; [m2 m-3 m-1] --> [um2 cm-3 m-1] Ordinate scale factor for surface area distribution
scl_ord_vlm_dst=1.0; [m3 m-3 m-1] --> [um3 cm-3 m-1] Ordinate scale factor for volume distribution
scl_ord_mss_dst=1.0; [kg m-3 m-1] --> [ug m-3 m-1] Ordinate scale factor for mass distribution
if scl_ord_flg then begin
	scl_ord_nbr_dst=1.0e-6 ; [# m-3 m-1] --> [# cm-3 m-1] Ordinate scale factor for number distribution
	scl_ord_sfc_dst=1.0e12*1.0e-6 ; [m2 m-3 m-1] --> [um2 cm-3 m-1] Ordinate scale factor for surface area distribution
	scl_ord_vlm_dst=1.0e18*1.0e-6 ; [m3 m-3 m-1] --> [um3 cm-3 m-1] Ordinate scale factor for volume distribution
	scl_ord_mss_dst=1.0e9 ; [kg m-3 m-1] --> [ug m-3 m-1] Ordinate scale factor for mass distribution
endif; endif scl_ord_flg

; scl_ord and scl_abc are independent of eachother (because scl_abc is always used, not always the case with scl_ord)
; scl_ord and scl_abc must both be set in order to plot data in traditional units of, e.g., um2 cm-3 um-1

scl_mss_usr=replicate(1.0,sz_nbr) ; [frc] Scale factor to user-specified mass concentrations
mss_rsl_usr=0.0*fltarr(sz_nbr) ; [kg m-3] Specified mass concentration resolved
if n_elements(mmr_rsl) ne 0 then begin
; User specified absolute mass amounts in each bin
; Convert specified aerosol mass mixing ratio to mass concentration
	mss_rsl_usr=mmr_rsl*dns_mdp ; [kg m-3] Specified mass concentration resolved
; Factor which scales input distribution to specified mass concentration
; Input distribution is usually, but not always, based on 1 particle m-3 number concentration (default in mie())
; mss_rsl is actual mass concentration of input distribution (computed by mie())
; mss_rsl_usr is user-specified mass concentration (derived from observations or global model DSTQ fields)
; Ratio of mss_rsl_usr to mss_rsl is factor by which mass-intensive quantities must be multiplied to obtain absolute quantitites
; All size distribution fields plotted by psd() are mass-intensive, and are defined in terms of dst, the number distribution
; Scale dst so distributions defined in terms of dst, e.g., surface area, will also be absolute distributions
	for psd_idx=0,psd_nbr-1 do begin
		scl_mss_usr(sz_idx_srt(psd_idx):sz_idx_end(psd_idx))=mss_rsl_usr(psd_idx)/mss_rsl(psd_idx) ; [frc] Scale factor to user-specified mass concentrations
	endfor; end loop over fl
endif; endif mss_rsl

if dbg ne 0 then begin
; 	fxm: Using reasonable mass concentrations from dust simulations, e.g., at Barbados,
; 	this model shows ridiculously low number concentrations, e.g., < 1.0 m-3
	if n_elements(mmr_rsl) ne 0 then print,'psd: mmr_rsl = ',mmr_rsl
	print,'psd: cnc_nbr_rsl_ttl = ',total(cnc_nbr_rsl)
	print,'psd: cnc_nbr_anl_ttl = ',total(cnc_nbr_anl)
	print,'psd: mss_rsl_ttl = ',total(mss_rsl)
	print,'psd: mss_anl_ttl = ',total(mss_anl)
	print,'psd: cnc_nbr_rsl = ',cnc_nbr_rsl
	print,'psd: cnc_nbr_anl = ',cnc_nbr_anl
	print,'psd: mss_rsl = ',mss_rsl
	print,'psd: mss_anl = ',mss_anl
	print,'psd: mss_rsl_usr = ',mss_rsl_usr
	print,'psd: rds_ctr = ',rds_ctr
	print,'psd: dst = ',dst
	print,'psd: cnc = ',cnc
	print,'psd: scl_mss_usr = ',scl_mss_usr
	print,'ord=dst*scl_ord_nbr_dst*scl_mss_usr*scl_fnc/scl_abc for sz_idx_dbg = ',sz_idx_dbg
endif; endif dbg

abc=rds_ctr*scl_abc
abc_min=min(abc)
abc_max=max(abc)
rng_x=[abc_min,abc_max]
ord=dst*scl_ord_nbr_dst*scl_mss_usr*scl_fnc/scl_abc
if nrm_flg then ord=ord/max(ord)
if dbg then print,auto_sng(ord(sz_idx_dbg),2),' = ',auto_sng(dst(sz_idx_dbg),2),'*',auto_sng(scl_ord_nbr_dst,2),'*',auto_sng(scl_mss_usr(sz_idx_dbg),2),'*',auto_sng(scl_fnc(sz_idx_dbg),2),'/',auto_sng(scl_abc,2)
ord_min=min(ord)
ord_max=max(ord)
rng_y=[ord_min,ord_max]
chr_sz=2.0
ttl='Number, Surface Area, Mass'
if dmt_flg then x_ttl='!5Diameter !8D!5 (!7l!5m)' else x_ttl='!5Radius !8r!5 (!7l!5m)'
if log_flg then dst_sym='!8n!S!5!IN!R!E0!N' else if ln_flg then dst_sym='!8n!S!5!IN!R!Ee!N' else dst_sym='!8n!5!IN!N'
if scl_ord_flg then begin
	y_ttl=dst_sym+' (# !5cm!E-3!N'
endif else begin
	y_ttl=dst_sym+' (# !5m!E-3!N'
endelse; endelse scl_ord_flg
if not prn then y_ttl='!5Number '+y_ttl
if log_flg or ln_flg then y_ttl=y_ttl+')' else y_ttl=y_ttl+' !7l!5m!E-1!N)'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	if fmt_nsm then top=1
	if fmt_nsm then btm=0
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	x_sz=6.5
	y_sz=3.0
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz+0.73
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/psd_dst'+fl_lbl+'.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot psd
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=12+x_sty,ystyle=y_sty,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
tvlct,r_curr,g_curr,b_curr

if btm then ntt=1 else ntt=0
gnr_axz_drw,ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,tck_lng=0.02,dbg_lvl=0,axz_vrt=0,axz_clr=clr_blk_idx

oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx
for psd_idx=0,psd_nbr-1 do begin
	if stt_flg then plots,[rds_nwr(psd_idx)*scl_abc,rds_nwr(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=0,thick=0.5,color=clr_blk_idx,/DATA
	if stt_flg then plots,[rds_nmr(psd_idx)*scl_abc,rds_nmr(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=0,thick=0.5,color=clr_blk_idx,/DATA
	if psd_idx ne 0 then plots,[rds_grd(psd_idx)*scl_abc,rds_grd(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=1,thick=0.5,color=clr_blk_idx,/DATA
endfor; end loop over fl

lbl_sng=strarr(psd_nbr)
ln_sty=0.0*intarr(psd_nbr)
clr=intarr(psd_nbr)
clr(0)=clr_blk_idx
if psd_nbr ge 2 then clr(1)=grn_idx
if psd_nbr ge 3 then clr(2)=blu_idx
lbl_sng(0)='Size 01'
if psd_nbr ge 2 then lbl_sng(1)='Size 02'
if psd_nbr ge 3 then lbl_sng(2)='Size 03'

if info then begin
ln_lgn_x1=0.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

xyouts,ln_lgn_x1,lgn_y(0)+lgn_dy,'!7s!5 = 1.0, !7h!5 = 30!E!12_!5!N',color=clr_blk_idx,size=chr_sz,/NORMAL

for psd_idx=0,psd_nbr-1 do begin
	plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(psd_idx)+0.013,linestyle=ln_sty(psd_idx),color=clr(psd_idx),thick=3.0,/NORMAL
	xyouts,txt_lgn_x,lgn_y(psd_idx),lbl_sng(psd_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

ord=sfc*dst*scl_ord_sfc_dst*scl_mss_usr*scl_fnc/scl_abc
if nrm_flg then ord=ord/max(ord)
rng_y=[min(ord),max(ord)]
ttl=''
if dmt_flg then x_ttl='!5Diameter !8D!5 (!7l!5m)' else x_ttl='!5Radius !8r!5 (!7l!5m)'
if log_flg then dst_sym='!8n!S!5!IS!R!E0!N' else if ln_flg then dst_sym='!8n!S!5!IS!R!Ee!N' else dst_sym='!8n!5!IS!N'
if scl_ord_flg then begin
	y_ttl=dst_sym+' (!7l!5m!E2!N !5cm!E-3!N'
endif else begin
	y_ttl=dst_sym+' (!5m!E2!N !5m!E-3!N'
endelse; endelse scl_ord_flg
if not prn then y_ttl='!5Surface Area '+y_ttl
if log_flg or ln_flg then y_ttl=y_ttl+')' else y_ttl=y_ttl+' !7l!5m!E-1!N)'

if prn then begin
	if fmt_nsm then top=0
	if fmt_nsm then btm=0
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	x_sz=6.5
	y_sz=3.0
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz+0.73
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/psd_sfc'+fl_lbl+'.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot sfc
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=12+x_sty,ystyle=y_sty,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

if btm then ntt=1 else ntt=0
gnr_axz_drw,ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,tck_lng=0.02,dbg_lvl=0,axz_vrt=0,axz_clr=clr_blk_idx

oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx
for psd_idx=0,psd_nbr-1 do begin
	if stt_flg then plots,[rds_swr(psd_idx)*scl_abc,rds_swr(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=0,thick=0.5,color=clr_blk_idx,/DATA
	if stt_flg then plots,[rds_smr(psd_idx)*scl_abc,rds_smr(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=0,thick=0.5,color=clr_blk_idx,/DATA
	if psd_idx ne 0 then plots,[rds_grd(psd_idx)*scl_abc,rds_grd(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=1,thick=0.5,color=clr_blk_idx,/DATA
endfor; end loop over fl

if info then begin
for psd_idx=1,psd_nbr-1 do begin
	plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(psd_idx-1)+0.013,linestyle=ln_sty(psd_idx),color=clr(psd_idx),thick=3.0,/NORMAL
	xyouts,txt_lgn_x,lgn_y(psd_idx-1),lbl_sng(psd_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

ord=vlm*dst*scl_ord_vlm_dst*scl_mss_usr*scl_fnc/scl_abc
if nrm_flg then ord=ord/max(ord)
;rng_x=[0.1,10.0]
rng_y=[min(ord),max(ord)]
ttl=''
if dmt_flg then x_ttl='!5Diameter !8D!5 (!7l!5m)' else x_ttl='!5Radius !8r!5 (!7l!5m)'
if log_flg then dst_sym='!8n!S!5!IV!R!E0!N' else if ln_flg then dst_sym='!8n!S!5!IV!R!Ee!N' else dst_sym='!8n!5!IV!N'
if scl_ord_flg then begin
	y_ttl=dst_sym+' (!7l!5m!E3!N !5cm!E-3!N'
endif else begin
	y_ttl=dst_sym+' (!5m!E3!N !5m!E-3!N'
endelse; endelse scl_ord_flg
if not prn then y_ttl='!5Volume '+y_ttl
if log_flg or ln_flg then y_ttl=y_ttl+')' else y_ttl=y_ttl+' !7l!5m!E-1!N)'

if prn then begin
	if fmt_nsm then top=0
	if fmt_nsm then btm=0
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	x_sz=6.5
	y_sz=3.0
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz+0.73
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/psd_vlm'+fl_lbl+'.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot volume
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=12+x_sty,ystyle=y_sty,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

if btm then ntt=1 else ntt=0
gnr_axz_drw,ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,tck_lng=0.02,dbg_lvl=0,axz_vrt=0,axz_clr=clr_blk_idx

oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx
for psd_idx=0,psd_nbr-1 do begin
	if stt_flg then plots,[rds_vwr(psd_idx)*scl_abc,rds_vwr(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=0,thick=0.5,color=clr_blk_idx,/DATA
	if stt_flg then plots,[rds_vmr(psd_idx)*scl_abc,rds_vmr(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=0,thick=0.5,color=clr_blk_idx,/DATA
	if psd_idx ne 0 then plots,[rds_grd(psd_idx)*scl_abc,rds_grd(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=1,thick=0.5,color=clr_blk_idx,/DATA
endfor; end loop over fl

lbl_sng(0)='Volume 01'
if psd_nbr ge 2 then lbl_sng(1)='Volume 02'
if psd_nbr ge 3 then lbl_sng(2)='Volume 03'

if info then begin
for psd_idx=1,psd_nbr-1 do begin
	plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(psd_idx-1)+0.013,linestyle=ln_sty(psd_idx),color=clr(psd_idx),thick=3.0,/NORMAL
	xyouts,txt_lgn_x,lgn_y(psd_idx-1),lbl_sng(psd_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

ord=mss*dst*scl_ord_mss_dst*scl_mss_usr*scl_fnc/scl_abc
if nrm_flg then ord=ord/max(ord)
;rng_x=[0.1,10.0]
rng_y=[min(ord),max(ord)]
ttl=''
if dmt_flg then x_ttl='!5Diameter !8D!5 (!7l!5m)' else x_ttl='!5Radius !8r!5 (!7l!5m)'
if log_flg then dst_sym='!8n!S!5!IM!R!E0!N' else if ln_flg then dst_sym='!8n!S!5!IM!R!Ee!N' else dst_sym='!8n!5!IM!N'
if scl_ord_flg then begin
	y_ttl=dst_sym+' (!7l!5g cm!E-3!N'
endif else begin
	y_ttl=dst_sym+' (!5kg m!E-3!N'
endelse; endelse scl_ord_flg
if not prn then y_ttl='!5Mass '+y_ttl
if log_flg or ln_flg then y_ttl=y_ttl+')' else y_ttl=y_ttl+' !7l!5m!E-1!N)'

if prn then begin
	if fmt_nsm then top=0
	if fmt_nsm then btm=1
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	x_sz=6.5
	y_sz=3.0
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz+0.73
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/psd_mss'+fl_lbl+'.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot mass
plot,abc,ord,xlog=xlg_flg,ylog=ylg_flg,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=12+x_sty,ystyle=y_sty,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

if btm then ntt=1 else ntt=0
gnr_axz_drw,ntt=ntt,ttl=x_ttl,chr_sz=chr_sz,tck_lng=0.02,dbg_lvl=0,axz_vrt=0,axz_clr=clr_blk_idx

oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx
for psd_idx=0,psd_nbr-1 do begin
; Plot volume weighted size since they equal mass-weighted size
	if stt_flg then plots,[rds_vwr(psd_idx)*scl_abc,rds_vwr(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=0,thick=0.5,color=clr_blk_idx,/DATA
	if stt_flg then plots,[rds_vmr(psd_idx)*scl_abc,rds_vmr(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=0,thick=0.5,color=clr_blk_idx,/DATA
	if psd_idx ne 0 then plots,[rds_grd(psd_idx)*scl_abc,rds_grd(psd_idx)*scl_abc],[!y.crange(0),!y.crange(1)],linestyle=1,thick=0.5,color=clr_blk_idx,/DATA
endfor; end loop over fl

lbl_sng(0)='Mass 01'
if psd_nbr ge 2 then lbl_sng(1)='Mass 02'
if psd_nbr ge 3 then lbl_sng(2)='Mass 03'

if info then begin
for psd_idx=1,psd_nbr-1 do begin
	plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(psd_idx-1)+0.013,linestyle=ln_sty(psd_idx),color=clr(psd_idx),thick=3.0,/NORMAL
	xyouts,txt_lgn_x,lgn_y(psd_idx-1),lbl_sng(psd_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

ord=vlc_grv
abc=rds_ctr*scl_abc
rng_x=[0.01,100.0]
rng_y=[min(ord),max(ord)]
ttl='Stokes Settling Velocity'
if dmt_flg then x_ttl='!5Diameter !8D!5 (!7l!5m)' else x_ttl='!5Radius !8r!5 (!7l!5m)'
y_ttl='!5Velocity !8V!5!Ig!N (!5m s!E-1!N)'

if prn then begin
	fl_nm_out=getenv('DATA')+'/ps/psd_vlc_grv.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot settling velocity
plot,abc,ord,xlog=xlg_flg,ylog=1,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx,yticklen=1,xticklen=1

oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx
for psd_idx=0,psd_nbr-1 do begin
endfor; end loop over fl

lbl_sng(0)='Velocity 01'
if psd_nbr ge 2 then lbl_sng(1)='Velocity 02'
if psd_nbr ge 3 then lbl_sng(2)='Velocity 03'

if info then begin
for psd_idx=1,psd_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(psd_idx-1)+0.013,linestyle=ln_sty(psd_idx),color=clr(psd_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(psd_idx-1),lbl_sng(psd_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

ord=scv_cff*3600.0
abc=rds_ctr*scl_abc
rng_x=[0.001,20.0]
;rng_x=[min(abc),max(abc)]
rng_y=[min(ord),max(ord)]
ttl='Scavenging Coefficient'
if dmt_flg then x_ttl='!5Diameter !8D!5 (!7l!5m)' else x_ttl='!5Radius !8r!5 (!7l!5m)'
y_ttl='!5Scav. Coeff !7K!5 (!5hr!E-1!N)'

if prn then begin
	fl_nm_out=getenv('DATA')+'/ps/psd_scv_cff.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot scavenging coefficient
plot,abc,ord,xlog=xlg_flg,ylog=1,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx,yticklen=1,xticklen=1

oplot,abc,ord,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx
for psd_idx=0,psd_nbr-1 do begin
endfor; end loop over fl

lbl_sng(0)='scv_cff 01'
if psd_nbr ge 2 then lbl_sng(1)='scv_cff 02'
if psd_nbr ge 3 then lbl_sng(2)='scv_cff 03'

if info then begin
for psd_idx=1,psd_nbr-1 do begin
;plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(psd_idx-1)+0.013,linestyle=ln_sty(psd_idx),color=clr(psd_idx),thick=3.0,/NORMAL
;xyouts,txt_lgn_x,lgn_y(psd_idx-1),lbl_sng(psd_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
dmt_pcp_nma_sng=auto_sng(dmt_pcp_nma(0)*1.0e6,0)
dmt_pcp_dbg_sng=auto_sng(dmt_pcp_dbg(0)*1.0e6,0)
flx_vlm_pcp_rsl_sng=auto_sng(flx_vlm_pcp_rsl(0)*1000.0*3600.0,1)
xyouts,txt_lgn_x,lgn_y(0),'Raindrop !8D!5 = '+dmt_pcp_dbg_sng+' !7l!5m',color=clr_blk_idx,size=chr_sz,/NORMAL
xyouts,txt_lgn_x,lgn_y(1),'Precip.  !8P!5 = '+flx_vlm_pcp_rsl_sng+' !5mm hr!E-1!N',color=clr_blk_idx,size=chr_sz,/NORMAL
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

ord=cll_fsh
abc=rds_ctr*scl_abc
rng_x=[0.001,20.0]
;rng_x=[min(abc),max(abc)]
rng_y=[min(ord),max(ord)]
ttl='Collision Efficiency'
if dmt_flg then x_ttl='!5Particle Diameter !8D!5 (!7l!5m)' else x_ttl='!5Particle Radius !8r!5 (!7l!5m)'
y_ttl='!5Efficiency !8E!5'

if prn then begin
	fl_nm_out=getenv('DATA')+'/ps/psd_cll_fsh.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot scavenging coefficient
plot,abc,ord,xlog=xlg_flg,ylog=1,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx,yticklen=0.02,xticklen=0.02

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

cll_nbr=4
lbl_sng=strarr(cll_nbr)
ln_sty=0.0*intarr(cll_nbr)
clr=intarr(cll_nbr)
clr(0)=clr_blk_idx
clr(1)=grn_idx
clr(2)=blu_idx
clr(3)=red_idx

oplot,rds_ctr*scl_abc,cll_fsh_brn_dff,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(1)
oplot,rds_ctr*scl_abc,cll_fsh_ntc,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(2)
oplot,rds_ctr*scl_abc,cll_fsh_mpc,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(3)
oplot,rds_ctr*scl_abc,cll_fsh,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)

lbl_sng(0)='Total'
if cll_nbr ge 2 then lbl_sng(1)='Brownian Diffusion'
if cll_nbr ge 3 then lbl_sng(2)='Interception'
if cll_nbr ge 4 then lbl_sng(3)='Impaction'

if info then begin
for cll_idx=0,cll_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(cll_idx)+0.013,linestyle=ln_sty(cll_idx),color=clr(cll_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(cll_idx),lbl_sng(cll_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
dmt_pcp_nma_sng=auto_sng(dmt_pcp_nma(0)*1.0e6,0)
dmt_pcp_dbg_sng=auto_sng(dmt_pcp_dbg(0)*1.0e6,0)
flx_vlm_pcp_rsl_sng=auto_sng(flx_vlm_pcp_rsl(0)*1000.0*3600.0,1)
xyouts,txt_lgn_x,lgn_y(0)+lgn_dy,'Raindrop !8D!5 = '+dmt_pcp_dbg_sng+' !7l!5m',color=clr_blk_idx,size=chr_sz,/NORMAL
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

ord=stk_nbr
abc=rds_ctr*scl_abc
rng_x=[0.001,20.0]
;rng_x=[min(abc),max(abc)]
rng_y=[min(ord),max(ord)]
ttl='Aerodynamic Parameters'
if dmt_flg then x_ttl='!5Particle Diameter !8D!5 (!7l!5m)' else x_ttl='!5Particle Radius !8r!5 (!7l!5m)'
y_ttl='!5Number !5St, St!E*!N, Re, Sh'

if prn then begin
	fl_nm_out=getenv('DATA')+'/ps/psd_stk_nbr.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot flow numbers
plot,abc,ord,xlog=xlg_flg,ylog=1,max_value=1.0e20,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx,yticklen=1,xticklen=1

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

cll_nbr=4
lbl_sng=strarr(cll_nbr)
ln_sty=0.0*intarr(cll_nbr)
clr=intarr(cll_nbr)
clr(0)=clr_blk_idx
clr(1)=grn_idx
clr(2)=blu_idx
clr(3)=red_idx

oplot,rds_ctr*scl_abc,stk_nbr,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(0)
oplot,rds_ctr*scl_abc,stk_nbr_rlt,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(1)
oplot,rds_ctr*scl_abc,shm_nbr,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(2)
oplot,rds_ctr*scl_abc,ryn_nbr_grv,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(3)
oplot,rds_pcp*scl_abc,stk_nbr_crt,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(3)
oplot,rds_pcp*scl_abc,ryn_nbr_pcp,max_value=1.0e20,thick=2.0,linestyle=0,color=clr(3)
;print,'stk_nbr_crt = ',stk_nbr_crt

lbl_sng(0)='Stokes number'
if cll_nbr ge 2 then lbl_sng(1)='Relative Stokes number'
if cll_nbr ge 3 then lbl_sng(2)='Schmidt number'
if cll_nbr ge 4 then lbl_sng(3)='Reynolds number'

if info then begin
for cll_idx=0,cll_nbr-1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(cll_idx)+0.013,linestyle=ln_sty(cll_idx),color=clr(cll_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(cll_idx),lbl_sng(cll_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
dmt_pcp_nma_sng=auto_sng(dmt_pcp_nma(0)*1.0e6,0)
dmt_pcp_dbg_sng=auto_sng(dmt_pcp_dbg(0)*1.0e6,0)
flx_vlm_pcp_rsl_sng=auto_sng(flx_vlm_pcp_rsl(0)*1000.0*3600.0,1)
xyouts,txt_lgn_x,lgn_y(0)+lgn_dy,'Raindrop !8D!5 = '+dmt_pcp_dbg_sng+' !7l!5m',color=clr_blk_idx,size=chr_sz,/NORMAL
endif; endif info

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end psd()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure psd
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro prn_tst, $
	info=info, $
	prn=prn, $ 
	scl=scl, $ 
	ngl_nbr=ngl_nbr, $
	dbg=dbg
; Purpose: Test the printing and display of IDL graphics

; prn_tst,dbg=1,ngl_nbr=91,prn=0
; prn_tst,prn=0

@~/idl/ibp_clr.com
!p.multi=0

if n_elements(ngl_nbr) eq 0 then ngl_nbr=91
if n_elements(scl) eq 0 then scl=0

if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22
if n_elements(dbg) eq 0 then dbg=0

; Some initialization
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=6.5
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

; ngl goes from 0 to 0.5*pi
ngl=0.5*!pi*findgen(ngl_nbr)/(ngl_nbr-1.0)
ngl_dgr=180.0*ngl/!pi
ngl_dlt=0.5*!pi*(fltarr(ngl_nbr)+1.0)/(ngl_nbr-1.0)
cos_ngl=cos(ngl)
alb_ocn_drc=0.026/(0.065+cos_ngl^1.7) + 0.15*(cos_ngl-0.10)*(cos_ngl-0.50)*(cos_ngl-1.0)
alb_ocn_dff=fltarr(ngl_nbr)+0.06

abc=ngl_dgr
abc_min=min(abc)
abc_max=max(abc)
rng_x=[abc_min,abc_max]
ord=alb_ocn_drc
ord_min=min(ord)
ord_max=max(ord)
rng_y=[ord_min,ord_max]
chr_sz=2.0
ttl='CCM Ocean Surface Albedo'
x_ttl='!5Solar Zenith Angle !7h!5 (degrees)'
y_ttl='!5Ocean Surface Albedo !8A!5'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/prn_tst.eps'
	open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot alb
plot,abc,ord,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=1,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

oplot,abc,alb_ocn_drc,max_value=1.0e20,thick=2.0,linestyle=0,color=clr_blk_idx
oplot,abc,alb_ocn_dff,max_value=1.0e20,thick=2.0,linestyle=0,color=grn_idx

psd_nbr=4
lbl_sng=strarr(psd_nbr)
ln_sty=0.0*intarr(psd_nbr)
clr=intarr(psd_nbr)
clr(0)=clr_blk_idx
clr(1)=grn_idx
clr(2)=blu_idx
clr(3)=red_idx
lbl_sng(0)='Direct'
lbl_sng(1)='Diffuse'

if info then begin
ln_lgn_x1=0.30
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.75
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

for psd_idx=0,1 do begin
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(psd_idx)+0.013,linestyle=ln_sty(psd_idx),color=clr(psd_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(psd_idx),lbl_sng(psd_idx),color=clr_blk_idx,size=chr_sz,/NORMAL
endfor; end loop over fl
endif; not info

print,'!d.name = ',!d.name
print,'!d.n_colors = ',!d.n_colors
print,'!d.table_size = ',!d.table_size
print,'clr_blk_idx = ',clr_blk_idx
print,'clr_wht_idx = ',clr_wht_idx
;print,'cbar_clr_nbr = ',cbar_clr_nbr
;print,'cbar_txt_clr = ',cbar_txt_clr
print,'n_elements(r_curr) = ',n_elements(r_curr)
print,'Current color map: '
print,'idx red grn blu'
for idx=0,4 do begin
	print,string(format='(4(I3,1x))',idx,r_curr(idx),g_curr(idx),b_curr(idx))
endfor; end loop over idx
for idx=!d.table_size-4,!d.table_size-1 do begin
	print,string(format='(4(I3,1x))',idx,r_curr(idx),g_curr(idx),b_curr(idx))
endfor; end loop over idx
print,'Original color map: '
print,'idx red grn blu'
for idx=0,4 do begin
	print,string(format='(4(I3,1x))',idx,r_orig(idx),g_orig(idx),b_orig(idx))
endfor; end loop over idx
for idx=!d.table_size-4,!d.table_size-1 do begin
	print,string(format='(4(I3,1x))',idx,r_orig(idx),g_orig(idx),b_orig(idx))
endfor; end loop over idx

if dbg then begin
fl_nm_out='~/idl/dbg.tbl'
prn_unit=1 
openw,prn_unit,fl_nm_out
for idx=0,n_elements(r_curr)-1 do begin
printf,prn_unit,string(format='(3(I3,1x))',r_curr(idx),g_curr(idx),b_curr(idx))
endfor; end loop over idx
close,prn_unit
endif; endif dbg

if prn then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end prn_tst()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure scv_gph: Graphs various aerosol scavenging properties from mie output
; Procedure reads in an ensemble of files, each containing a different
; scavenging/deposition experiment but all on the same grid.
; Procedure then plots sensitivity experiments to tested parameters
; Procedure works for both wet scavenging and dry deposition, but not at same time
; Since parameters of interest are usually different for these processes,
; must specify which process with dps=dry or dps=wet (default) on command line
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro scv_gph, $
	fl_nm=fl_nm, $
	mie=mie, $
	dbg=dbg, $
	dps=dps, $
	info=info, $
	prn=prn
; Usage:
; scv_gph,fl_nm=getenv('DATA')+'/mie/mie.nc'
; scv_gph,dps='wet',prn=1 
; scv_gph,dps='dry',prn=1
@~/idl/ibp_clr.com
!p.multi=0

if n_elements(mie) eq 0 then mie=0
if n_elements(dbg) eq 0 then dbg=0
if n_elements(prn) eq 0 then prn=0
if n_elements(info) eq 0 then info=1
if n_elements(palette) eq 0 then palette=4
if n_elements(clr_tbl) eq 0 then clr_tbl=22
if n_elements(dps) eq 0 then dps='wet'

;fl_nm=[getenv('DATA')+'/mie/mie.nc']
fl_nm=[ $
	getenv('DATA')+'/mie/scv_cnv.nc', $
	getenv('DATA')+'/mie/scv_str.nc', $
	getenv('DATA')+'/mie/scv_cnv_mono.nc', $
	getenv('DATA')+'/mie/scv_str_mono.nc' $
;	getenv('DATA')+'/mie/scv_tst.nc' $
;	getenv('DATA')+'/mie/dps_01.nc', $
;	getenv('DATA')+'/mie/dps_02.nc', $
;	getenv('DATA')+'/mie/dps_03.nc', $
;	getenv('DATA')+'/mie/dps_04.nc' $
;	getenv('DATA')+'/mie/dps_tst.nc' $
	] ; end fl_nm[]
fl_nbr=n_elements(fl_nm)

; Some initialization
sym_usr_foo=findgen(17)*(!pi*2/16.)
color_order=0
clr_rfr,clr_tbl,"",0
clr_mk,10,0,palette,0
tvlct,r_curr,g_curr,b_curr
top=1
btm=1
if prn then begin
	x_sz=6.5
	y_sz=4.0
	if top then y_sz=y_sz*1.1
	if btm then y_sz=y_sz*1.1
endif; endif prn

; Get coordinate sizes from all files
sz_sz=intarr(fl_nbr) ; [nbr] Number of sub-gridscale bins
print,'Scanning all files for dimension sizes...'
for fl_idx=0,fl_nbr-1 do begin
	nc_id=ncdf_open(fl_nm(fl_idx))
	dim_id=ncdf_dimid(nc_id,'sz')
	ncdf_diminq,nc_id,dim_id,foo,sz_nbr
	sz_sz(fl_idx)=sz_nbr
	ncdf_close,nc_id
endfor; end loop over fl

sz_sz_max=max(sz_sz)
; 20000905 Forgetting to dimension by both dimensions may crash X server!
scv_cff_pcp_nrm=fltarr(sz_sz_max,fl_nbr) ; [m2 kg-1] Scavenging coefficient, precipitation-normalized
vlc_dry=fltarr(sz_sz_max,fl_nbr) ; [m s-1] Total dry deposition velocity
vlc_grv=fltarr(sz_sz_max,fl_nbr) ; [m s-1] Settling velocity
vlc_trb=fltarr(sz_sz_max,fl_nbr) ; [m s-1] Turbulent deposition velocity
dmt_vmr=fltarr(fl_nbr) ; [m] Diameter volume median resolved
rgh_mmn_dps=fltarr(fl_nbr) ; [m] Roughness length momentum
dns_aer=fltarr(fl_nbr) ; [kg m-3] Aerosol density
wnd_frc_dps=fltarr(fl_nbr) ; [m s-1] Surface friction velocity
dmt_pcp_nma=fltarr(fl_nbr) ; [m] Diameter number median analytic, raindrop
dmt_nmr=fltarr(fl_nbr) ; [m] Diameter number median resolved
scv_cff_mss_avg_pcp_nrm=fltarr(fl_nbr) ; [m2 kg-1] Mass mean scavenging coefficient of aerosol, precipitation-normalized
gsd_pcp_anl=fltarr(fl_nbr) ; [frc] Geometric standard deviation, raindrop
flx_vlm_pcp_rsl=fltarr(fl_nbr) ; [m s-1]=[m3 m-2 s-1] Precipitation volume flux, resolved
for fl_idx=0,fl_nbr-1 do begin
	print,'Ingesting ',fl_nm(fl_idx)
	nc_id=ncdf_open(fl_nm(fl_idx))
	ncdf_varget_1D,nc_id,'dmt_nmr',dmt_nmr,[fl_idx,fl_idx]
	ncdf_varget_1D,nc_id,'wnd_frc_dps',wnd_frc_dps,[fl_idx,fl_idx]
	ncdf_varget_1D,nc_id,'rgh_mmn_dps',rgh_mmn_dps,[fl_idx,fl_idx]
	ncdf_varget_1D,nc_id,'dns_aer',dns_aer,[fl_idx,fl_idx]
	ncdf_varget_1D,nc_id,'scv_cff_mss_avg_pcp_nrm',scv_cff_mss_avg_pcp_nrm,[fl_idx,fl_idx]
	ncdf_varget_1D,nc_id,'dmt_vmr',dmt_vmr,[fl_idx,fl_idx]
	ncdf_varget_1D,nc_id,'dmt_pcp_nma',dmt_pcp_nma,[fl_idx,fl_idx]
	ncdf_varget_1D,nc_id,'gsd_pcp_anl',gsd_pcp_anl,[fl_idx,fl_idx]
	ncdf_varget_1D,nc_id,'flx_vlm_pcp_rsl',flx_vlm_pcp_rsl,[fl_idx,fl_idx]
	ncdf_varget_2D,nc_id,'scv_cff_pcp_nrm',scv_cff_pcp_nrm,fl_idx
	ncdf_varget_2D,nc_id,'vlc_dry',vlc_dry,fl_idx
	ncdf_varget_2D,nc_id,'vlc_trb',vlc_trb,fl_idx
	ncdf_varget_2D,nc_id,'vlc_grv',vlc_grv,fl_idx
	if fl_idx eq 0 then ncdf_varget,nc_id,'dmt_ctr',dmt_ctr
	ncdf_close,nc_id
endfor; end loop over fl

abc=dmt_ctr*1.0e6
abc_min=min(abc)
abc_max=max(abc)
rng_x=[abc_min,abc_max]
ord=scv_cff_pcp_nrm
ord_min=min(ord)
ord_max=max(ord)
;rng_y=[ord_min,ord_max]
rng_y=[1.0e-4,1.0e0]
if prn then chr_sz=2.0 else chr_sz=1
;ttl=string(aer_lng_nm)
ttl='Normalized Scavenging'
x_ttl='!5Diameter !8D!5 (!7l!5m)'
y_ttl='!5Scav. Coeff. !7K!5 (!5mm!E-1!N)'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/scv_gph.eps'
	if dps eq 'wet' then open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot scv_gph
plot,abc,ord(*,0),xlog=1,ylog=1,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=0,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx

base_clr_idx=2
grn_idx=base_clr_idx+0
blu_idx=base_clr_idx+1
red_idx=base_clr_idx+2
yll_idx=base_clr_idx+3
orn_idx=base_clr_idx+4
r_curr(grn_idx)=0
g_curr(grn_idx)=255
b_curr(grn_idx)=0
r_curr(blu_idx)=0
g_curr(blu_idx)=0
b_curr(blu_idx)=255
r_curr(yll_idx)=255
g_curr(yll_idx)=255
b_curr(yll_idx)=0
r_curr(orn_idx)=255
g_curr(orn_idx)=128
b_curr(orn_idx)=0
r_curr(red_idx)=255
g_curr(red_idx)=0
b_curr(red_idx)=0
tvlct,r_curr,g_curr,b_curr

lbl_sng=strarr(fl_nbr)
ln_sty=0.0*intarr(fl_nbr)
clr=intarr(fl_nbr)
if fl_nbr gt 0 then clr(0)=clr_blk_idx
if fl_nbr gt 1 then clr(1)=grn_idx
if fl_nbr gt 2 then clr(2)=blu_idx
if fl_nbr gt 3 then clr(3)=yll_idx
if fl_nbr gt 4 then clr(4)=red_idx
if fl_nbr gt 5 then clr(5)=orn_idx
if fl_nbr gt 0 then lbl_sng(0)='Convective '
if fl_nbr gt 1 then lbl_sng(1)='Stratiform '
if fl_nbr gt 2 then lbl_sng(2)='Convective '
if fl_nbr gt 3 then lbl_sng(3)='Stratiform '
if fl_nbr gt 4 then lbl_sng(4)='Test'

usersym,cos(sym_usr_foo),sin(sym_usr_foo),/fill
for fl_idx=0,fl_nbr-1 do begin
	oplot,abc,ord(*,fl_idx),max_value=1.0e20,thick=2.0,linestyle=0,color=clr(fl_idx)
;	print,'ord(*,fl_idx) = ',ord(*,fl_idx)
	plots,dmt_nmr(fl_idx)*1.0e6,scv_cff_mss_avg_pcp_nrm(fl_idx),color=clr(fl_idx),psym=8
;	print,'scv_cff_mss_avg_pcp_nrm(',fl_idx,') = ',scv_cff_mss_avg_pcp_nrm(fl_idx)
endfor; end loop over fl

if info then begin
ln_lgn_x1=0.23
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.85
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

for fl_idx=0,fl_nbr-1 do begin
dsd_sng='!5D!In!N = '+auto_sng(dmt_pcp_nma(fl_idx)*1.0e6,0)+' !7l!5m, !7r!5!Ig!N = '+auto_sng(gsd_pcp_anl(fl_idx),2); +', !5P = '+auto_sng(flx_vlm_pcp_rsl(fl_idx)*3.6e6,1)+' !5mm hr!E-1!N'
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),dsd_sng,color=clr_blk_idx,size=chr_sz*0.5,/NORMAL
endfor; end loop over fl
endif; info

if prn and dps eq 'wet' then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

abc=dmt_ctr*1.0e6 ; [m] -> [um]
abc_min=min(abc)
abc_max=max(abc)
rng_x=[abc_min,abc_max]
rng_x=[abc_min,10.0]
ord=vlc_dry*100.0 ; [m s-1] -> [cm s-1]
ord_min=min(ord)
ord_max=max(ord)
;rng_y=[ord_min,ord_max]
rng_y=[ord_min,1.0e1]
if prn then chr_sz=2.0 else chr_sz=1
;ttl=string(aer_lng_nm)
ttl='Dry Deposition'
x_ttl='!5Diameter !8D!5 (!7l!5m)'
y_ttl='!5Dep. Vel. !8v!Id!N!5 (!5cm s!E-1!N)'

mrg_top=2.5 ; 2 is default
mrg_btm=1 ; 4 is default
if btm then mrg_btm=mrg_btm+2
mrg_lft=8 ; 10 is default
mrg_rgt=1.1 ; 3 is default
if prn then begin
	mrg_top=0.6 ; 2 is default
	if top then mrg_top=mrg_top+1.1
	mrg_btm=0.5 ; 4 is default
	if btm then mrg_btm=mrg_btm+2.7
	if not top then ttl=''
	if not btm then x_ttl=''
	fl_nm_out=getenv('DATA')+'/ps/vlc_dps_aer.eps'
	if dps eq 'dry' then open_ps,fl_nm=fl_nm_out,x_sz=x_sz,y_sz=y_sz,/eps,/color
endif; endif prn

; xyplot scv_gph
plot,abc,ord(*,0),xlog=1,ylog=1,tit=ttl,xtit=x_ttl,ytit=y_ttl,xrange=rng_x,yrange=rng_y,xmargin=[mrg_lft,mrg_rgt],ymargin=[mrg_btm,mrg_top],xstyle=0,ystyle=0,thick=2.0,charsize=chr_sz,linestyle=0,/nodata,color=clr_blk_idx,yticklen=1,xticklen=1

lbl_sng=strarr(fl_nbr)
ln_sty=0.0*intarr(fl_nbr)
clr=intarr(fl_nbr)
if fl_nbr gt 0 then clr(0)=clr_blk_idx
if fl_nbr gt 1 then clr(1)=grn_idx
if fl_nbr gt 2 then clr(2)=blu_idx
if fl_nbr gt 3 then clr(3)=yll_idx
if fl_nbr gt 4 then clr(4)=red_idx
if fl_nbr gt 5 then clr(5)=orn_idx
if fl_nbr gt 0 then lbl_sng(0)='Hi '
if fl_nbr gt 1 then lbl_sng(1)='Hi '
if fl_nbr gt 2 then lbl_sng(2)='Hi '
if fl_nbr gt 3 then lbl_sng(3)='Hi '
if fl_nbr gt 4 then lbl_sng(4)='Test'

usersym,cos(sym_usr_foo),sin(sym_usr_foo),/fill
for fl_idx=0,fl_nbr-1 do begin
	oplot,abc,ord(*,fl_idx),max_value=1.0e20,thick=2.0,linestyle=0,color=clr(fl_idx)
	oplot,abc,vlc_grv(*,fl_idx)*100.0,max_value=1.0e20,thick=2.0,linestyle=2,color=clr(fl_idx)
	oplot,abc,vlc_trb(*,fl_idx)*100.0,max_value=1.0e20,thick=2.0,linestyle=1,color=clr(fl_idx)
endfor; end loop over fl

if info then begin
ln_lgn_x1=0.23
ln_lgn_dx=0.07
ln_lgn_x2=ln_lgn_x1+ln_lgn_dx
txt_lgn_x=ln_lgn_x2+0.01
txt_lgn_sz=1.75
lgn_y_top=0.85
lgn_dy=0.1
lgn_y=lgn_y_top-indgen(6)*lgn_dy

for fl_idx=0,fl_nbr-1 do begin
dps_sng='!7q!5 = '+auto_sng(dns_aer(fl_idx)*0.001,1)+' !5g cm!e-3!N, !5z!I0!N = '+auto_sng(rgh_mmn_dps(fl_idx)*100.0,2)+' cm, !5u!I*!N = '+auto_sng(wnd_frc_dps(fl_idx)*100.0,0)+' !5cm s!E-1!N'
plots,[ln_lgn_x1,ln_lgn_x2],lgn_y(fl_idx)+0.013,linestyle=ln_sty(fl_idx),color=clr(fl_idx),thick=3.0,/NORMAL
xyouts,txt_lgn_x,lgn_y(fl_idx),dps_sng,color=clr_blk_idx,size=chr_sz*0.5,/NORMAL
endfor; end loop over fl
endif; info

if prn and dps eq 'dry' then close_ps,fl_nm=fl_nm_out else begin
print,'Hit any key to continue, or q to quit ...'
foo=get_kbrd(1)
if foo eq 'q' then goto,end_of_procedure
endelse; endif not prn

end_of_procedure: foo=1

end; end scv_gph()
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of Figure scv_gph
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
